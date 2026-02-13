"""Core PyRosetta computation: FastRelax, residue scanning, ddG calculation."""

from __future__ import annotations

import logging
import os
from random import randint

import numpy as np
import pandas as pd
from pyrosetta import *  # type: ignore[misc]
from pyrosetta.rosetta import *  # type: ignore[misc]
from scipy.spatial import cKDTree

logger = logging.getLogger("residue_scanning")

# Flags used for every pyrosetta.init() call — main process and pool workers.
PYROSETTA_FLAGS = (
    "-use_input_sc -flip_HNQ -ex1 -ex2 -relax:cartesian "
    "-ignore_unrecognized_res -mute all -ignore_zero_occupancy false"
)


def _worker_init() -> None:
    """Pool initializer: call pyrosetta.init() once per spawned worker process.

    On macOS (and any platform using the 'spawn' start method) each worker
    process starts with a fresh Python interpreter that has NOT called
    pyrosetta.init().  Without this initializer the workers hang waiting
    to unpickle a Pose that requires Rosetta to be set up.
    """
    pyrosetta.init(PYROSETTA_FLAGS)  # noqa: F405 — pyrosetta imported via star-import


AA = "ADEFGHIKLMNQRSTYVW"  # 18 amino acids — Cys and Pro excluded

THREE_TO_ONE: dict[str, str] = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLU": "E",
    "GLN": "Q",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}
ONE_TO_THREE: dict[str, str] = {v: k for k, v in THREE_TO_ONE.items()}


# ---------------------------------------------------------------------------
# Utilities
# ---------------------------------------------------------------------------


def random_seed_pyrosetta(debug: bool = False) -> None:
    seed = randint(0, 9_999_999)
    pyrosetta.rosetta.basic.random.init_random_generators(seed, "mt19937")
    if debug:
        print(f"Seed is {seed}")


def wildtyper(mut: str) -> str:
    """Return the wildtype name for a mutation string.

    e.g. 'GH122AL' -> 'GH122AG'  (last char replaced by first char)
    Double mutations: 'GH122AL_SH130T' -> 'GH122AG_SH130SS'
    """
    return "_".join(f"{m[:-1]}{m[0]}" for m in mut.replace(",", "_").split("_"))


def score_fxn_checker(not_relax: bool = False):
    """Return ref2015_cart (with minimization) or ref2015 (without)."""
    name = "ref2015" if not_relax else "ref2015_cart"
    return pyrosetta.rosetta.core.scoring.ScoreFunctionFactory.create_score_function(name)


# ---------------------------------------------------------------------------
# Residue contact detection
# ---------------------------------------------------------------------------


def _extract_ligand_coords(pose, ligand_resids: list[int]) -> np.ndarray:
    """Pre-compute an (N, 3) array of all ligand atom coordinates."""
    coords: list[list[float]] = []
    for resid in ligand_resids:
        res = pose.residue(resid)
        for atom_idx in range(1, res.natoms() + 1):
            xyz = res.xyz(atom_idx)
            coords.append([xyz.x, xyz.y, xyz.z])
    return np.array(coords)


def is_interface_residue(
    pose,
    receptor_resid: int,
    ligand_tree: cKDTree,
    cutoff: float = 8.0,
) -> bool:
    """Two-tier hybrid filter for interface residue detection.

    Uses a pre-built ``cKDTree`` of all ligand atom coordinates for fast
    nearest-neighbour queries.

    **Tier 1 — direct contact** (SC_min ≤ ``DIRECT_CONTACT``):
        Any sidechain atom within 4.0 Å of the nearest ligand atom.
        Accepted unconditionally — proven interaction regardless of
        backbone orientation.

    **Tier 2 — oriented contact** (``DIRECT_CONTACT`` < SC_min ≤ ``ORIENTED_MAX``):
        Sidechain atom within 5.0 Å, but only if:
        (a) CB is oriented toward the ligand (CA→CB vs CA→nearest-ligand
            angle < 90°), AND
        (b) backbone CA is within ``cutoff + 1`` Å of the ligand (prevents
            long sidechains reaching across from distant backbone).

    **Glycine** — CA within ``DIRECT_CONTACT`` Å (mutation from Gly adds a
        sidechain; only truly embedded glycines are included).

    Thresholds calibrated so that the default cutoff of 8.0 Å captures all
    26 manually curated interface positions on the 5JXE test case with zero
    false positives.
    """
    DIRECT_CONTACT = 4.0  # Å — unconditional sidechain contact
    ORIENTED_MAX = min(DIRECT_CONTACT + 1.0, cutoff)  # 5.0 Å at default
    CA_MAX = cutoff + 1.0  # 9.0 Å at default

    rec_res = pose.residue(receptor_resid)

    # --- Gather sidechain coords ---
    sc_coords: list[list[float]] = []
    for atom_idx in range(1, rec_res.natoms() + 1):
        if not rec_res.atom_is_backbone(atom_idx):
            xyz = rec_res.xyz(atom_idx)
            sc_coords.append([xyz.x, xyz.y, xyz.z])

    ca_xyz = rec_res.xyz(rec_res.atom_index("CA"))
    ca = np.array([ca_xyz.x, ca_xyz.y, ca_xyz.z])

    # --- Glycine proxy ---
    if not sc_coords:
        dist, _ = ligand_tree.query(ca)
        return dist <= DIRECT_CONTACT

    # --- Tier 1: very close SC contact — definite interaction ---
    sc_arr = np.array(sc_coords)
    sc_dists, _ = ligand_tree.query(sc_arr)
    sc_min = float(np.min(sc_dists))
    if sc_min <= DIRECT_CONTACT:
        return True

    # --- Tier 2: medium SC contact — require orientation + backbone proximity ---
    if sc_min <= ORIENTED_MAX:
        ca_dist, _ = ligand_tree.query(ca)
        if ca_dist > CA_MAX:
            return False
        try:
            cb_idx = rec_res.atom_index("CB")
            cb_xyz = rec_res.xyz(cb_idx)
            cb = np.array([cb_xyz.x, cb_xyz.y, cb_xyz.z])
            _, cb_nearest_idx = ligand_tree.query(cb)
            nearest_lig = ligand_tree.data[cb_nearest_idx]
            v_cb = cb - ca
            v_lig = nearest_lig - ca
            norm_cb = np.linalg.norm(v_cb)
            norm_lig = np.linalg.norm(v_lig)
            if norm_cb > 0 and norm_lig > 0:
                cos_angle = np.dot(v_cb, v_lig) / (norm_cb * norm_lig)
                return cos_angle > 0  # angle < 90°
        except Exception:
            pass

    return False


# ---------------------------------------------------------------------------
# Mutation and minimization
# ---------------------------------------------------------------------------


def mutate(pose, mutation: str, debug: bool = False):
    """Apply a single point mutation to *pose* in-place.

    Returns *(pose, rosetta_resid)* on success or *(KeyError, rosetta_resid)*
    if the residue identity in the name does not match the pose.
    """
    if debug:
        print(f"Step mutate: {pose.pdb_info().name()} {mutation}")

    chain = mutation[1]
    from_aa = mutation[0]
    to_aa = mutation[-1]

    if mutation[-2].isalpha():
        insertion_code = mutation[-2]
        position = int(mutation[2:-2])
        rosetta_resid = pose.pdb_info().pdb2pose(chain, position, insertion_code)
    else:
        position = int(mutation[2:-1])
        rosetta_resid = pose.pdb_info().pdb2pose(chain, position)

    if from_aa != pose.sequence()[rosetta_resid - 1]:
        print(f"{mutation} -- Incorrect mutation: AA in name ≠ AA in pose")
        return KeyError, rosetta_resid

    mutate_residue = pyrosetta.rosetta.protocols.simple_moves.MutateResidue()
    sel = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
    sel.set_index(rosetta_resid)
    mutate_residue.set_selector(sel)
    mutate_residue.set_res_name(ONE_TO_THREE[to_aa])
    mutate_residue.apply(pose)

    return pose, rosetta_resid


def interface_analyzer(pose, receptor_chains: str, ligand_chains: str, not_relax: bool = False):
    """Compute and return interface ΔG for the given pose."""
    iam = pyrosetta.rosetta.protocols.analysis.InterfaceAnalyzerMover()
    sfxn = score_fxn_checker(not_relax)
    iam.set_interface(f"{''.join(ligand_chains)}_{''.join(receptor_chains)}")
    iam.set_compute_packstat(False)
    iam.set_pack_input(False)
    iam.set_pack_separated(True)
    iam.set_scorefunction(sfxn)
    iam.apply(pose)
    return iam.get_interface_dG()


def repack_and_minimize(
    pose,
    replica: int,
    receptor_chains: str,
    ligand_chains: str,
    mutation=False,
    radius: float = 8.0,
    adjacent_aa_range: int = 1,
    not_relax: bool = False,
    debug: bool = False,
):
    """Repack (and optionally FastRelax) the pose around the mutated position."""
    if debug:
        print(f"Step repack_and_minimize: {pose.pdb_info().name()}")

    sfxn = score_fxn_checker(not_relax)

    all_sel = pyrosetta.rosetta.core.select.residue_selector.TrueResidueSelector()
    not_sel = pyrosetta.rosetta.core.select.residue_selector.NotResidueSelector()
    not_sel.set_residue_selector(all_sel)

    tf = pyrosetta.rosetta.core.pack.task.TaskFactory()
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.InitializeFromCommandline())
    tf.push_back(pyrosetta.rosetta.core.pack.task.operation.NoRepackDisulfides())
    tf.push_back(
        pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
            pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), not_sel
        )
    )
    tf.push_back(
        pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
            pyrosetta.rosetta.core.pack.task.operation.RestrictToRepackingRLT(), all_sel
        )
    )

    mmf = pyrosetta.rosetta.core.select.movemap.MoveMapFactory()
    mmf.set_cartesian(True)
    mmf.all_bb(False)
    mmf.all_bondangles(True)
    mmf.all_bondlengths(True)
    mmf.all_chi(True)
    mmf.all_jumps(False)

    if mutation is not False:
        mut_sel = pyrosetta.rosetta.core.select.residue_selector.ResidueIndexSelector()
        for m in mutation.split("_"):
            pose, rosetta_resid = mutate(pose=pose, mutation=m, debug=debug)
            mut_sel.append_index(rosetta_resid)

        adj_sel = (
            pyrosetta.rosetta.core.select.residue_selector.PrimarySequenceNeighborhoodSelector()
        )
        adj_sel.set_selector(mut_sel)
        adj_sel.set_lower_residues(adjacent_aa_range)
        adj_sel.set_upper_residues(adjacent_aa_range)
        mmf.add_bb_action(pyrosetta.rosetta.core.select.movemap.mm_enable, adj_sel)

        nbr_sel = pyrosetta.rosetta.core.select.residue_selector.NeighborhoodResidueSelector()
        nbr_sel.set_distance(radius)
        nbr_sel.set_focus_selector(mut_sel)
        nbr_sel.set_include_focus_in_subset(True)
        tf.push_back(
            pyrosetta.rosetta.core.pack.task.operation.OperateOnResidueSubset(
                pyrosetta.rosetta.core.pack.task.operation.PreventRepackingRLT(), nbr_sel, True
            )
        )

    if not_relax or debug:
        packer = pyrosetta.rosetta.protocols.minimization_packing.PackRotamersMover()
        packer.score_function(sfxn)
        packer.task_factory(tf)
        packer.apply(pose)
    else:
        fr = pyrosetta.rosetta.protocols.relax.FastRelax()
        fr.set_scorefxn(sfxn)
        fr.set_task_factory(tf)
        fr.set_movemap_factory(mmf)
        fr.max_iter(200)
        fr.cartesian(True)
        fr.apply(pose)

    pdb_base = pose.pdb_info().name().replace(".pdb", "")
    if mutation is not False:
        pose.dump_pdb(f"{pdb_base}_{mutation}_{replica}.pdb")
    else:
        pose.dump_pdb(f"{pdb_base}_{replica}.pdb")

    dG_interface = interface_analyzer(pose, receptor_chains, ligand_chains, not_relax)
    return sfxn.score(pose), dG_interface


def mutate_pose_FR(
    original_pose,
    replica: int,
    mutation,
    receptor_chains: str,
    ligand_chains: str,
    not_relax: bool = False,
    debug: bool = False,
    radius: float = 8.0,
):
    """Clone pose, apply mutation, repack/minimize, return (score, dG_interface)."""
    pose = original_pose.clone()
    if debug:
        print(f"Step mutate_pose_FR: {pose.pdb_info().name()} {mutation}")
    sfxn = score_fxn_checker(not_relax)
    sfxn.score(pose)
    random_seed_pyrosetta(debug)
    return repack_and_minimize(
        pose=pose,
        replica=replica,
        mutation=mutation,
        receptor_chains=receptor_chains,
        ligand_chains=ligand_chains,
        not_relax=not_relax,
        debug=debug,
        radius=radius,
    )


# ---------------------------------------------------------------------------
# Dataframe-level scanning
# ---------------------------------------------------------------------------


def FastRelax_replics(
    pose,
    replics: int,
    mutation: str,
    receptor_chains: str,
    ligand_chains: str,
    not_relax: bool = False,
    debug: bool = False,
    radius: float = 8.0,
) -> list:
    """Run up to *replics* stochastic repack/relax trials for one mutation.

    Stops early (after replica 2) when both complex scores agree within 1 REU
    and full relaxation is enabled.
    """
    if debug:
        print(f"Step FastRelax_replics: {pose.pdb_info().name()}")

    score_mut: list = ["NaN"] * (replics * 2)

    for replica in range(1, replics + 1):
        score_mut[replica - 1], score_mut[replica - 1 + replics] = mutate_pose_FR(
            original_pose=pose,
            replica=replica,
            mutation=mutation,
            receptor_chains=receptor_chains,
            ligand_chains=ligand_chains,
            not_relax=not_relax,
            debug=debug,
            radius=radius,
        )
        if replica == 2:
            if not_relax:
                continue
            if abs(score_mut[0] - score_mut[1]) < 1:  # type: ignore[operator]
                break

    if debug:
        print(mutation, score_mut)
    return score_mut


def df_FastRelax(
    df: pd.DataFrame,
    pose,
    replics: int,
    receptor_chains: str,
    ligand_chains: str,
    not_relax: bool = False,
    debug: bool = False,
    radius: float = 8.0,
) -> pd.DataFrame:
    """Apply FastRelax_replics to every row in *df* and store results."""
    if debug:
        print(f"Step df_FastRelax: {pose.pdb_info().name()}")

    df[list(df.columns[1:])] = df.apply(
        lambda x: FastRelax_replics(
            pose=pose,
            replics=replics,
            mutation=x.Name,
            receptor_chains=receptor_chains,
            ligand_chains=ligand_chains,
            not_relax=not_relax,
            debug=debug,
            radius=radius,
        ),
        axis=1,
    ).to_list()

    return df


def ddG_calculation(x, df: pd.DataFrame):
    """Compute ΔΔG_complex and ΔΔG_interface for one mutant row."""
    wt_name = wildtyper(x.Name)
    wt_match = df.loc[df["Name"] == wt_name, "REU_min"]

    # Handle case where wildtype row might be missing or has NaN values
    if wt_match.empty or pd.isna(wt_match.iat[0]):
        logger.warning(f"Wildtype not found or unevaluated for {x.Name} (expected: {wt_name})")
        return float("nan"), float("nan")

    wt_REU = wt_match.iat[0]
    wt_dG = df.loc[df["Name"] == wt_name, "Min_dG_Interface"].iat[0]
    return x.REU_min - wt_REU, x.Min_dG_Interface - wt_dG


def df_ddG_postprocessing(
    df_concatenated: pd.DataFrame,
    pose,
    replics: int,
    postfix: str = "",
    ambiguity_reports: list[str] | None = None,
) -> pd.DataFrame:
    """Compute REU_min, Min_dG_Interface, ΔΔG values; move PDB files; write CSV.

    If ambiguity_reports are provided, saves them to a separate report file.
    """
    df1 = df_concatenated.copy()
    for col in df1.columns[1:]:
        df1[col] = pd.to_numeric(df1[col], errors="coerce")

    complex_cols = [f"Replica_{i+1}_Complex" for i in range(replics)]

    df1.loc[:, "REU_min"] = df1[complex_cols].min(axis=1)
    df1.loc[:, "Min_structure"] = (
        df1[complex_cols].idxmin(axis=1).str.extract(r"Replica_(\d+)_Complex")[0].astype(int)
    )
    df1.loc[:, "Min_dG_Interface"] = df1.apply(
        lambda x: x[f"Replica_{int(x.Min_structure)}_dG"], axis=1
    )

    df2 = df1.loc[:, ["Name", "REU_min", "Min_dG_Interface"]].copy()
    df2.loc[:, "Is_WT"] = df2.Name.apply(lambda x: x == wildtyper(x))

    pdb_base = pose.pdb_info().name().replace(".pdb", "")

    all_dir = f"all_REU_pdb{postfix}"
    min_dir = f"min_REU_pdb{postfix}"
    if not os.path.exists(all_dir):
        os.mkdir(all_dir)
    os.system(f"mv {pdb_base}_* {all_dir}")

    if not os.path.exists(min_dir):
        os.mkdir(min_dir)
    df1.apply(
        lambda x: os.system(f"cp {all_dir}/{pdb_base}_{x.Name}_{x.Min_structure}.pdb {min_dir}"),
        axis=1,
    )

    ddg_results = df2.apply(
        lambda x: ddG_calculation(x, df2), axis=1
    )
    df2.loc[:, "ddG_complex"] = [x[0] for x in ddg_results]
    df2.loc[:, "ddG_interface"] = [x[1] for x in ddg_results]

    df3 = (
        df2.query("Is_WT != True")
        .loc[:, ["Name", "REU_min", "ddG_complex", "ddG_interface"]]
        .reset_index(drop=True)
        .sort_values("ddG_interface")
    )
    csv_name = f"Rosetta_ddG_mut{postfix}.csv"
    df3.to_csv(csv_name, index=False)
    logger.info("Results written to %s (%d mutant rows)", csv_name, len(df3))

    # Save ambiguity reports if any were detected
    if ambiguity_reports:
        report_name = f"CSV_AMBIGUITIES{postfix}.txt"
        with open(report_name, "w") as f:
            f.write("CSV PARSING AMBIGUITIES DETECTED AND RESOLVED\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"Total ambiguities: {len(ambiguity_reports)}\n\n")
            for i, report in enumerate(ambiguity_reports, 1):
                f.write(f"{i}. {report}\n\n")
            f.write("\nEXPLANATION:\n")
            f.write("-" * 80 + "\n")
            f.write(
                "When a CSV entry has a single trailing character that is a valid amino acid\n"
                "(e.g., 'GH100A'), it creates an ambiguity:\n\n"
                "  - Mutation interpretation: G→A at position H100 (no insertion code)\n"
                "  - Position spec interpretation: All 18 mutations at position H100A (with insertion code A)\n\n"
                "The parser checks if both positions exist in the PDB:\n"
                "  - If BOTH exist: generates mutations for BOTH interpretations\n"
                "  - If only one exists: uses that interpretation\n"
                "  - If neither exist: skips the entry\n\n"
                "This ensures no mutations are missed due to ambiguous CSV format.\n"
            )
        logger.info("Ambiguity report written to %s", report_name)

    return df3


# ---------------------------------------------------------------------------
# Interface residue selection
# ---------------------------------------------------------------------------


def prepare_mut_df(
    pose,
    receptor_chains: str,
    ligand_chains: str,
    replics: int,
    radius: float = 8.0,
    debug: bool = False,
) -> pd.DataFrame:
    """Build a DataFrame of all single-point mutations at the interface.

    Every receptor-chain residue that passes the hybrid interface filter
    (``is_interface_residue``) within *radius* Å is included.  Uses a
    pre-built ``cKDTree`` of ligand atom coordinates for fast queries.
    """
    chains: dict[str, list] = {}
    for i in range(1, pose.num_chains() + 1):
        cs, ce = pose.chain_begin(i), pose.chain_end(i)
        name = pose.pdb_info().chain(cs)
        chains[name] = [i, cs, ce]

    logger.debug("Chain map: %s", chains)

    # Collect all ligand residue indices
    ligand_resids: list[int] = []
    for chain in ligand_chains:
        cs, ce = chains[chain][1], chains[chain][2]
        ligand_resids.extend(range(cs, ce + 1))

    # Build KDTree from all ligand atom coordinates
    ligand_coords = _extract_ligand_coords(pose, ligand_resids)
    ligand_tree = cKDTree(ligand_coords)

    # Collect all receptor residue indices
    receptor_resids: list[int] = []
    for chain in receptor_chains:
        cs, ce = chains[chain][1], chains[chain][2]
        receptor_resids.extend(range(cs, ce + 1))

    # Keep only receptor positions that pass the hybrid interface filter
    receptor_interface = [
        pos
        for pos in receptor_resids
        if is_interface_residue(pose, pos, ligand_tree, cutoff=radius)
    ]

    # Log selected interface positions
    position_names: list[str] = []
    screen_mut: list[list] = []
    for position in receptor_interface:
        res = pose.residue(position)
        from_aa = res.name1()
        pdb_pos = pose.pdb_info().pose2pdb(position).split()
        pdb_number = pdb_pos[0]
        chain = pdb_pos[1]
        insertion_code = pose.pdb_info().icode(position).replace(" ", "")
        position_names.append(f"{from_aa}{chain}{pdb_number}{insertion_code}")
        for to_aa in AA:
            row = [f"{from_aa}{chain}{pdb_number}{insertion_code}{to_aa}"]
            row.extend(["NaN"] * (replics * 2))
            screen_mut.append(row)

    n_positions = len(receptor_interface)
    n_mutations = len(screen_mut)
    logger.info(
        "Interface positions detected: %d (from %d receptor residues, cutoff=%.1f A)",
        n_positions,
        len(receptor_resids),
        radius,
    )
    logger.info("Positions: %s", ", ".join(position_names))
    logger.info(
        "Mutations to evaluate: %d (%d positions x %d amino acids)",
        n_mutations,
        n_positions,
        len(AA),
    )

    col_names = ["Name"]
    col_names += [f"Replica_{i+1}_Complex" for i in range(replics)]
    col_names += [f"Replica_{i+1}_dG" for i in range(replics)]

    return pd.DataFrame(screen_mut, columns=col_names)


def prepare_custom_mut_df(
    pose,
    csv_path,
    receptor_chains: str,
    ligand_chains: str,
    replics: int,
    debug: bool = False,
) -> tuple[pd.DataFrame, list[str]]:
    """Build a DataFrame of mutations from a user-supplied CSV file.

    The CSV should have a single column with entries in one of these formats::

        AH98    — position only; generates all 18 standard mutations
        AH98D   — specific mutation (from=A, chain=H, number=98, to=D)
        AH98bD  — specific mutation with insertion code 'b'
        AH98b   — position with insertion code; generates all 18 mutations

    Trailing-character rules (characters after chain + number digits):
        0 chars  → position, generate 18 mutations
        1 char that is a valid target AA → specific mutation (but see AMBIGUITY below)
        1 char NOT a valid target AA    → insertion code, generate 18 mutations
        2 chars  → first is insertion code, second is to_aa

    AMBIGUITY HANDLING:
        If a single trailing char is a valid AA (e.g., 'A' in 'GH100A'):
        - Checks if BOTH interpretations exist in the PDB:
            1. Mutation: G→A at position H100 (no insertion code)
            2. Position spec: All 18 mutations at position H100A (insertion code)
        - If both exist: generates BOTH sets and marks as ambiguous in report
        - If only one exists: uses that interpretation
        - Reports ambiguity detected and resolved

    Returns:
        tuple of (DataFrame, list of ambiguity messages)

    Validates that from_aa matches the actual residue in the pose.
    Automatically adds wildtype entries for ddG calculation.
    """
    import re
    from pathlib import Path

    csv_path = Path(csv_path)
    lines = csv_path.read_text().strip().splitlines()

    mutation_re = re.compile(r"^([A-Z])([A-Za-z])(\d+)(.*)$")
    all_valid_aa = set("ACDEFGHIKLMNPQRSTVWY")

    entries: list[tuple[str, str, str, str]] = []
    for i, line in enumerate(lines):
        line = line.strip().strip(",")
        if not line:
            continue
        m = mutation_re.match(line)
        if m is None:
            if i == 0:
                continue  # skip header row
            print(f"WARNING: skipping unrecognized entry at line {i + 1}: {line!r}")
            continue
        entries.append(m.groups())  # type: ignore[arg-type]

    screen_mut: list[list] = []
    seen_wt: set[str] = set()
    ambiguity_reports: list[str] = []

    for from_aa, chain, number, trailing in entries:
        number_int = int(number)

        # --- parse trailing chars to get insertion_code + to_aa(s) ---
        if len(trailing) == 0:
            insertion_code = ""
            to_aas = list(AA)
        elif len(trailing) == 1:
            ch = trailing[0]
            if ch.isupper() and ch in all_valid_aa:
                # AMBIGUOUS: could be mutation OR insertion code
                # Check both interpretations
                mutation_interpretation = _check_mutation_interpretation(
                    pose, from_aa, chain, number_int, ch
                )
                insertion_code_interpretation = _check_insertion_code_interpretation(
                    pose, from_aa, chain, number_int, ch
                )

                if mutation_interpretation and insertion_code_interpretation:
                    # Both exist! Process both and report ambiguity
                    ambiguity_msg = (
                        f"AMBIGUITY DETECTED in '{from_aa}{chain}{number}{ch}': "
                        f"Both mutation ({from_aa}→{ch} at {chain}{number}) and "
                        f"position spec (all mutations at {chain}{number}{ch.upper()}) exist in PDB. "
                        f"Processing BOTH interpretations."
                    )
                    print(f"ℹ️  {ambiguity_msg}")
                    ambiguity_reports.append(ambiguity_msg)

                    # Process as BOTH mutation AND position spec
                    _process_mutation_entry(
                        pose, from_aa, chain, number_int, "", [ch], replics, screen_mut, seen_wt
                    )
                    _process_mutation_entry(
                        pose,
                        from_aa,
                        chain,
                        number_int,
                        ch.upper(),
                        list(AA),
                        replics,
                        screen_mut,
                        seen_wt,
                    )
                    continue

                elif mutation_interpretation:
                    # Only mutation interpretation valid
                    insertion_code = ""
                    to_aas = [ch]
                elif insertion_code_interpretation:
                    # Only insertion code interpretation valid
                    insertion_code = ch
                    to_aas = list(AA)
                else:
                    # Neither exists
                    print(
                        f"WARNING: neither interpretation of '{from_aa}{chain}{number}{ch}' "
                        f"found in pose, skipping"
                    )
                    continue
            else:
                insertion_code = ch
                to_aas = list(AA)
        elif len(trailing) == 2:
            insertion_code = trailing[0]
            to_aas = [trailing[1]]
        else:
            print(
                f"WARNING: cannot parse trailing '{trailing}' in "
                f"{from_aa}{chain}{number}{trailing}, skipping"
            )
            continue

        # Process normal (non-ambiguous) cases
        _process_mutation_entry(
            pose, from_aa, chain, number_int, insertion_code, to_aas, replics, screen_mut, seen_wt
        )

    if debug:
        print(f"{len(screen_mut)} mutations generated from custom CSV")
        if ambiguity_reports:
            print(f"\n{len(ambiguity_reports)} ambiguities detected and resolved:")
            for msg in ambiguity_reports:
                print(f"  - {msg}")

    col_names = ["Name"]
    col_names += [f"Replica_{i+1}_Complex" for i in range(replics)]
    col_names += [f"Replica_{i+1}_dG" for i in range(replics)]

    return pd.DataFrame(screen_mut, columns=col_names), ambiguity_reports


def _check_mutation_interpretation(pose, from_aa: str, chain: str, number: int, to_aa: str) -> bool:
    """Check if '{from_aa}{chain}{number}' exists (mutation interpretation)."""
    try:
        rosetta_resid = pose.pdb_info().pdb2pose(chain, number)
        if rosetta_resid == 0:
            return False
        actual_aa = pose.residue(rosetta_resid).name1()
        return actual_aa == from_aa
    except Exception:
        return False


def _check_insertion_code_interpretation(
    pose, from_aa: str, chain: str, number: int, insertion_code: str
) -> bool:
    """Check if '{chain}{number}{insertion_code}' exists (insertion code interpretation)."""
    try:
        rosetta_resid = pose.pdb_info().pdb2pose(chain, number, insertion_code)
        if rosetta_resid == 0:
            return False
        actual_aa = pose.residue(rosetta_resid).name1()
        return actual_aa == from_aa
    except Exception:
        return False


def _process_mutation_entry(
    pose,
    from_aa: str,
    chain: str,
    number: int,
    insertion_code: str,
    to_aas: list[str],
    replics: int,
    screen_mut: list[list],
    seen_wt: set[str],
) -> None:
    """Process a single mutation entry (helper function)."""
    # --- validate position exists in pose ---
    if insertion_code:
        rosetta_resid = pose.pdb_info().pdb2pose(chain, number, insertion_code)
    else:
        rosetta_resid = pose.pdb_info().pdb2pose(chain, number)

    if rosetta_resid == 0:
        print(f"WARNING: residue {chain}{number}{insertion_code} not found " f"in pose, skipping")
        return

    # --- validate from_aa matches actual residue ---
    actual_aa = pose.residue(rosetta_resid).name1()
    if from_aa != actual_aa:
        print(
            f"WARNING: {from_aa}{chain}{number}{insertion_code} — "
            f"expected {from_aa} but pose has {actual_aa}, skipping"
        )
        return

    # --- generate mutation rows ---
    for to_aa in to_aas:
        name = f"{from_aa}{chain}{number}{insertion_code}{to_aa}"
        row = [name]
        row.extend(["NaN"] * (replics * 2))
        screen_mut.append(row)

    # --- ensure wildtype entry exists for this position ---
    wt_name = f"{from_aa}{chain}{number}{insertion_code}{from_aa}"
    if wt_name not in seen_wt:
        seen_wt.add(wt_name)
        # only add if not already generated above
        if from_aa not in to_aas:
            wt_row = [wt_name]
            wt_row.extend(["NaN"] * (replics * 2))
            screen_mut.append(wt_row)


def _parse_dm_component(pose, part: str) -> tuple[str, str, int, str, list[str]] | None:
    """Parse one side of a DM CSV entry (e.g., 'AH99' or 'AH99D').

    Returns:
        (from_aa, chain, number_int, insertion_code, to_aas) or None if parse fails.
    """
    import re

    mutation_re = re.compile(r"^([A-Z])([A-Za-z])(\d+)(.*)$")
    all_valid_aa = set("ACDEFGHIKLMNPQRSTVWY")

    m = mutation_re.match(part)
    if m is None:
        print(f"WARNING: cannot parse DM component '{part}', skipping")
        return None

    from_aa, chain, number, trailing = m.groups()
    number_int = int(number)

    # --- parse trailing chars ---
    if len(trailing) == 0:
        insertion_code = ""
        to_aas = list(AA)
    elif len(trailing) == 1:
        ch = trailing[0]
        if ch.isupper() and ch in all_valid_aa:
            # Ambiguous: check both interpretations
            mutation_interpretation = _check_mutation_interpretation(
                pose, from_aa, chain, number_int, ch
            )
            insertion_code_interpretation = _check_insertion_code_interpretation(
                pose, from_aa, chain, number_int, ch
            )

            if mutation_interpretation and insertion_code_interpretation:
                # Both exist: process mutation (specific) not position spec
                insertion_code = ""
                to_aas = [ch]
            elif mutation_interpretation:
                insertion_code = ""
                to_aas = [ch]
            elif insertion_code_interpretation:
                insertion_code = ch
                to_aas = list(AA)
            else:
                print(f"WARNING: neither interpretation of '{from_aa}{chain}{number}{ch}' found in pose, skipping")
                return None
        else:
            insertion_code = ch
            to_aas = list(AA)
    elif len(trailing) == 2:
        insertion_code = trailing[0]
        to_aas = [trailing[1]]
    else:
        print(f"WARNING: cannot parse trailing '{trailing}' in {from_aa}{chain}{number}{trailing}, skipping")
        return None

    return (from_aa, chain, number_int, insertion_code, to_aas)


def prepare_double_mut_df(
    pose,
    csv_path,
    replics: int,
    debug: bool = False,
) -> pd.DataFrame:
    """Build a DataFrame of double mutations from a user-supplied CSV file.

    The CSV should have a single column with entries in one of these formats::

        AH99_GL101         — position pair; generates all 18×18 combinations
        AH99D_GL101H       — specific double mutation
        AH99_GL101H        — mixed position and specific mutation
        AH99D_GL101        — mixed specific and position

    Returns:
        DataFrame with same column format as `df_FastRelax` input (Name + replicas columns)
    """
    from pathlib import Path

    csv_path = Path(csv_path)
    lines = csv_path.read_text().strip().splitlines()

    double_mut: list[list] = []
    seen_wt: set[str] = set()

    for i, line in enumerate(lines):
        line = line.strip().strip(",")
        if not line:
            continue
        if "_" not in line:
            print(f"WARNING: DM entry at line {i + 1} has no underscore: '{line}', skipping")
            continue

        parts = line.split("_", 1)
        if len(parts) != 2:
            print(f"WARNING: cannot split DM entry '{line}' into two parts, skipping")
            continue

        # Parse both sides
        result1 = _parse_dm_component(pose, parts[0])
        result2 = _parse_dm_component(pose, parts[1])

        if result1 is None or result2 is None:
            continue

        from_aa1, chain1, number1, ic1, to_aas1 = result1
        from_aa2, chain2, number2, ic2, to_aas2 = result2

        # Add wildtype pair FIRST so it's not filtered out by debug mode cutoff
        wt_pair = (
            f"{from_aa1}{chain1}{number1}{ic1}{from_aa1}_"
            f"{from_aa2}{chain2}{number2}{ic2}{from_aa2}"
        )
        if wt_pair not in seen_wt:
            seen_wt.add(wt_pair)
            wt_row = [wt_pair]
            wt_row.extend(["NaN"] * (replics * 2))
            double_mut.append(wt_row)

        # Generate all combinations
        for to_aa1 in to_aas1:
            for to_aa2 in to_aas2:
                name = (
                    f"{from_aa1}{chain1}{number1}{ic1}{to_aa1}_"
                    f"{from_aa2}{chain2}{number2}{ic2}{to_aa2}"
                )
                row = [name]
                row.extend(["NaN"] * (replics * 2))
                double_mut.append(row)

    if debug:
        print(f"{len(double_mut)} double mutations generated from CSV")

    col_names = ["Name"]
    col_names += [f"Replica_{i+1}_Complex" for i in range(replics)]
    col_names += [f"Replica_{i+1}_dG" for i in range(replics)]

    return pd.DataFrame(double_mut, columns=col_names)
