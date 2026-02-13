"""Command-line interface for residue-scan."""

from __future__ import annotations

import argparse
import logging
import multiprocessing as mp
import os
from itertools import repeat
from pathlib import Path

import numpy as np
import pandas as pd

logger = logging.getLogger("residue_scanning")


def _parse_args(argv=None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        prog="residue-scan",
        description="Protein interface residue scanning using PyRosetta",
    )
    p.add_argument("-f", "--file", required=True, help="Path to PDB file with complex")
    p.add_argument(
        "-m",
        "--mode",
        required=True,
        choices=[
            "FR",
            "FastRelax",
            "RS",
            "Residue_Scanning",
            "DM",
            "Double_Mut_Searching",
            "CM",
            "Custom_Mutations",
        ],
        help=(
            "FR/FastRelax — minimize structure only; "
            "RS/Residue_Scanning — compute ΔΔG for all interface mutations; "
            "DM/Double_Mut_Searching — RS then evaluate double-mutation combinations; "
            "CM/Custom_Mutations — compute ΔΔG for user-supplied mutations/positions"
        ),
    )
    p.add_argument(
        "-r",
        "--receptor",
        required=True,
        help="Receptor chain IDs without spaces (e.g. HL for chains H and L)",
    )
    p.add_argument("-l", "--ligand", required=True, help="Ligand chain IDs without spaces (e.g. A)")
    p.add_argument("--cpu", required=True, type=int, help="Number of CPU cores")
    p.add_argument(
        "-o",
        "--output",
        default="results",
        help="Output directory (default: 'results' next to the input file)",
    )
    p.add_argument(
        "--not_relax",
        action="store_true",
        help="Skip FastRelax (uses PackRotamers instead; faster but less accurate)",
    )
    p.add_argument(
        "--replics",
        type=int,
        default=5,
        help="Max number of stochastic replicas per mutation (default: 5)",
    )
    p.add_argument(
        "--radius",
        type=float,
        default=8.0,
        help="Neighbour shell radius in Å for repacking (default: 8.0)",
    )
    p.add_argument(
        "--condition",
        default="ddG_complex < 0.5 and ddG_interface < 0.5",
        help="Pandas query condition for selecting single mutations in DM mode",
    )
    p.add_argument(
        "--debug",
        action="store_true",
        help="Debug mode: limit mutation set to 4 positions for quick testing",
    )
    p.add_argument(
        "--mutations_csv",
        type=str,
        default=None,
        help="Path to CSV with custom mutations (AH98D) or positions (AH98). "
        "Required for CM mode.",
    )
    args = p.parse_args(argv)
    if args.mode in ("CM", "Custom_Mutations") and not args.mutations_csv:
        p.error("--mutations_csv is required when mode is CM/Custom_Mutations")
    return args


def main(argv=None) -> None:
    args = _parse_args(argv)

    # Import core first to get PYROSETTA_FLAGS, then init before any Rosetta objects.
    import pyrosetta  # noqa: PLC0415

    from .core import (
        AA,
        PYROSETTA_FLAGS,
        _worker_init,
        df_ddG_postprocessing,
        df_FastRelax,
        mutate_pose_FR,
        prepare_custom_mut_df,
        prepare_mut_df,
        score_fxn_checker,
    )
    from .preprocessing import prepare_pdb

    # ------------------------------------------------------------------
    # PyRosetta init — must happen before any Rosetta object is created
    # ------------------------------------------------------------------

    pyrosetta.init(PYROSETTA_FLAGS)

    # ------------------------------------------------------------------
    # Paths
    # ------------------------------------------------------------------
    input_file = Path(args.file).resolve()
    if Path(args.output).is_absolute():
        output_dir = Path(args.output)
    else:
        output_dir = input_file.parent / args.output
    output_dir.mkdir(parents=True, exist_ok=True)
    os.chdir(output_dir)

    # ------------------------------------------------------------------
    # Logging setup — console + file in output directory
    # ------------------------------------------------------------------
    logger.setLevel(logging.DEBUG if args.debug else logging.INFO)

    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_handler.setFormatter(logging.Formatter("%(message)s"))
    logger.addHandler(console_handler)

    file_handler = logging.FileHandler(output_dir / "residue_scan.log")
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(
        logging.Formatter(
            "%(asctime)s | %(levelname)-8s | %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
        )
    )
    logger.addHandler(file_handler)

    mode_names = {
        "FR": "FastRelax",
        "RS": "Residue_Scanning",
        "DM": "Double_Mut_Searching",
        "CM": "Custom_Mutations",
    }
    mode_display = mode_names.get(args.mode, args.mode)

    logger.info("=" * 60)
    logger.info("Residue Scanning Pipeline")
    logger.info("=" * 60)
    logger.info("Input PDB:       %s", input_file)
    logger.info("Output dir:      %s", output_dir)
    logger.info("Mode:            %s", mode_display)
    logger.info("Receptor chains: %s", args.receptor)
    logger.info("Ligand chains:   %s", args.ligand)
    logger.info("CPUs:            %d", args.cpu)
    logger.info("Replicas:        %d", args.replics)
    logger.info("Radius:          %.1f A", args.radius)
    logger.info(
        "Scoring:         %s",
        "ref2015 (--not_relax)" if args.not_relax else "ref2015_cart (FastRelax)",
    )
    logger.info("Debug mode:      %s", "ON" if args.debug else "OFF")
    if args.mutations_csv:
        logger.info("Mutations CSV:   %s", Path(args.mutations_csv).resolve())
    if args.mode in ("DM", "Double_Mut_Searching"):
        logger.info("DM condition:    %s", args.condition)
    logger.info("=" * 60)

    # ------------------------------------------------------------------
    # PDB preprocessing
    # ------------------------------------------------------------------
    pdb = prepare_pdb(input_file, output_dir)

    pose = pyrosetta.pose_from_pdb(str(pdb))
    original_pose = pose.clone()
    sfxn = score_fxn_checker(args.not_relax)
    initial_score = sfxn.score(pose)
    logger.info("Input structure:  %.2f REU", initial_score)

    receptor_chains: str = args.receptor
    ligand_chains: str = args.ligand
    replics: int = args.replics
    not_relax: bool = args.not_relax
    debug: bool = args.debug
    radius: float = args.radius
    ddg_condition: str = args.condition

    # ------------------------------------------------------------------
    # FastRelax block  (runs for ALL modes)
    # ------------------------------------------------------------------
    if args.mode in [
        "FR",
        "FastRelax",
        "RS",
        "Residue_Scanning",
        "DM",
        "Double_Mut_Searching",
        "CM",
        "Custom_Mutations",
    ]:
        relaxed_scores: dict[str, float] = {}
        lowest_energy = 0.0
        lowest_pose_path = ""

        with mp.Pool(replics, initializer=_worker_init) as pool:
            pool.starmap(
                mutate_pose_FR,
                [
                    (pose, i + 1, False, receptor_chains, ligand_chains, not_relax, debug, radius)
                    for i in range(replics)
                ],
            )

        for i in range(1, replics + 1):
            candidate = str(pdb).replace(".pdb", f"_{i}.pdb")
            if Path(candidate).exists():
                rpose = pyrosetta.pose_from_pdb(candidate)
                relaxed_scores[candidate] = sfxn.score(rpose)

        for path, energy in relaxed_scores.items():
            if energy < lowest_energy:
                lowest_energy = energy
                lowest_pose_path = path

        if not lowest_pose_path:
            # fallback: all scores were 0 or equal — take the first replica
            lowest_pose_path = list(relaxed_scores.keys())[0]
            lowest_energy = relaxed_scores[lowest_pose_path]

        logger.info("")
        logger.info("FastRelax: %s has the smallest REU", Path(lowest_pose_path).name)
        logger.info("  Original pose:  %.2f REU", sfxn.score(original_pose))
        logger.info("  Relaxed pose:   %.2f REU", lowest_energy)

        pose = pyrosetta.pose_from_pdb(lowest_pose_path)
        original_pose = pose.clone()
        sfxn.score(pose)

        fr_dir = output_dir / "FastRelax"
        fr_dir.mkdir(exist_ok=True)
        for f in output_dir.glob(f"{pdb.stem}_*.pdb"):
            f.rename(fr_dir / f.name)

        pose.dump_pdb("FastRelax.pdb")

    # ------------------------------------------------------------------
    # Residue Scanning block
    # ------------------------------------------------------------------
    if args.mode in ["RS", "Residue_Scanning", "DM", "Double_Mut_Searching"]:
        pose = pyrosetta.pose_from_pdb("FastRelax.pdb")
        original_pose = pose.clone()

        logger.info("")
        logger.info("Residue Scanning")
        logger.info("-" * 40)

        df = prepare_mut_df(pose, receptor_chains, ligand_chains, replics, radius, debug)

        # AA set already excludes C/P, but guard against from-AA being C or P
        df = df.loc[~(df["Name"].str.startswith("P") | df["Name"].str.startswith("C")), :]

        logger.info(
            "After Cys/Pro filter: %d mutations x %d replicas = %d structures",
            len(df),
            replics,
            len(df) * replics,
        )

        if debug:
            df = df.head(len(AA) * 4)

        n_proc = min(args.cpu, len(df))

        with mp.Pool(n_proc, initializer=_worker_init) as pool:
            df_split = np.array_split(df, n_proc)
            df_concat = pd.concat(
                pool.starmap(
                    df_FastRelax,
                    zip(
                        df_split,
                        repeat(pose),
                        repeat(replics),
                        repeat(receptor_chains),
                        repeat(ligand_chains),
                        repeat(not_relax),
                        repeat(debug),
                        repeat(radius),
                    ),
                )
            )

        df_concat.to_csv("Rosetta_results_REU.csv", index=False)
        df3 = df_ddG_postprocessing(df_concat, pose, replics)

    # ------------------------------------------------------------------
    # Custom Mutations block
    # ------------------------------------------------------------------
    if args.mode in ["CM", "Custom_Mutations"]:
        pose = pyrosetta.pose_from_pdb("FastRelax.pdb")
        original_pose = pose.clone()

        logger.info("")
        logger.info("Custom Mutations")
        logger.info("-" * 40)

        mutations_csv = Path(args.mutations_csv).resolve()
        df = prepare_custom_mut_df(
            pose,
            mutations_csv,
            receptor_chains,
            ligand_chains,
            replics,
            debug,
        )

        logger.info(
            "%d mutations x %d replicas = %d structures",
            len(df),
            replics,
            len(df) * replics,
        )

        if debug:
            df = df.head(len(AA) * 4)

        n_proc = min(args.cpu, len(df))

        with mp.Pool(n_proc, initializer=_worker_init) as pool:
            df_split = np.array_split(df, n_proc)
            df_concat = pd.concat(
                pool.starmap(
                    df_FastRelax,
                    zip(
                        df_split,
                        repeat(pose),
                        repeat(replics),
                        repeat(receptor_chains),
                        repeat(ligand_chains),
                        repeat(not_relax),
                        repeat(debug),
                        repeat(radius),
                    ),
                )
            )

        df_concat.to_csv("Rosetta_results_REU.csv", index=False)
        df3 = df_ddG_postprocessing(df_concat, pose, replics)

    # ------------------------------------------------------------------
    # Double Mutation block
    # ------------------------------------------------------------------
    if args.mode in ["DM", "Double_Mut_Searching"]:
        if debug:
            selected_positions = df3.Name.apply(lambda x: x[1:-1]).unique().tolist()  # type: ignore[union-attr]
        else:
            selected_positions = (
                df3.query(ddg_condition).Name.apply(lambda x: x[1:-1]).unique().tolist()  # type: ignore[union-attr]
            )

        chains_map: dict[str, list] = {}
        for i in range(1, pose.num_chains() + 1):
            cs, ce = pose.chain_begin(i), pose.chain_end(i)
            cname = pose.pdb_info().chain(cs)
            chains_map[cname] = [i, cs, ce]

        neighbors: dict[str, list[str]] = {}
        for pos in selected_positions:
            chain = pos[0]
            if pos[-1].isalpha():
                insertion_code = pos[-1]
                pdb_num = int(pos[1:-1])
                pyrosetta_pos = pose.pdb_info().pdb2pose(chain, pdb_num, insertion_code)
            else:
                pdb_num = int(pos[1:])
                pyrosetta_pos = pose.pdb_info().pdb2pose(chain, pdb_num)

            res_sel = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
            res_sel.set_index(pyrosetta_pos)
            nbr_sel = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
            nbr_sel.set_focus_selector(res_sel)
            nbr_sel.set_include_focus_in_subset(False)

            for nbr in nbr_sel.selection_positions(pose):
                nbr_pdb = pose.pdb_info().pose2pdb(nbr).split()
                nbr_chain = nbr_pdb[1]
                if nbr_chain in receptor_chains:
                    neighbors.setdefault(pos, []).append(f"{nbr_chain}{nbr_pdb[0]}")

        # Keep only mutual neighbours
        mutual: dict[str, list[str]] = {}
        for key, values in neighbors.items():
            for val in values:
                if val in neighbors and key not in mutual.get(val, []):
                    mutual.setdefault(key, []).append(val)

        if debug:
            df4 = df3.copy()  # type: ignore[union-attr]
        else:
            df4 = df3.query(ddg_condition).copy()  # type: ignore[union-attr]

        df4.loc[:, "Position"] = df4.Name.apply(lambda x: x[1:-1])
        df4.loc[:, "Mutate_from"] = df4.Name.apply(lambda x: x[0])
        df4.loc[:, "Mutate_to"] = df4.Name.apply(lambda x: x[-1])

        df5 = (
            df4.groupby("Position")
            .agg(lambda x: x.tolist())
            .loc[:, ["Mutate_to", "Mutate_from"]]
            .reset_index()
        )

        mut_combinations: list[list] = []
        for pos1, nbrs in mutual.items():
            from1 = df5.loc[df5["Position"] == pos1, "Mutate_from"].iat[0][0]
            for pos2 in nbrs:
                from2 = df5.loc[df5["Position"] == pos2, "Mutate_from"].iat[0][0]
                for to1 in df5.loc[df5["Position"] == pos1, "Mutate_to"].iat[0]:
                    for to2 in df5.loc[df5["Position"] == pos2, "Mutate_to"].iat[0]:
                        row = [f"{from1}{pos1}{to1}_{from2}{pos2}{to2}"]
                        row.extend(["NaN"] * (replics * 2))
                        mut_combinations.append(row)
                # Wildtype pair
                wt_row = [f"{from1}{pos1}{from1}_{from2}{pos2}{from2}"]
                wt_row.extend(["NaN"] * (replics * 2))
                mut_combinations.append(wt_row)

        logger.info("%d double mutation combinations", len(mut_combinations))

        col_names = ["Name"]
        col_names += [f"Replica_{i+1}_Complex" for i in range(replics)]
        col_names += [f"Replica_{i+1}_dG" for i in range(replics)]
        df6 = pd.DataFrame(mut_combinations, columns=col_names)

        n_proc = min(args.cpu, len(df6))

        with mp.Pool(n_proc, initializer=_worker_init) as pool:
            df_split = np.array_split(df6, n_proc)
            df_concat = pd.concat(
                pool.starmap(
                    df_FastRelax,
                    zip(
                        df_split,
                        repeat(pose),
                        repeat(replics),
                        repeat(receptor_chains),
                        repeat(ligand_chains),
                        repeat(not_relax),
                        repeat(debug),
                        repeat(radius),
                    ),
                )
            )

        df_concat.to_csv("Rosetta_results_REU_double_mut.csv", index=False)
        df_ddG_postprocessing(df_concat, pose, replics, postfix="_double")

    logger.info("")
    logger.info("DONE")


if __name__ == "__main__":
    main()
