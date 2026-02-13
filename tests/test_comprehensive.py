"""Comprehensive test suite for residue scanning pipeline.

Tests:
1. Interface residue selection on 4ZFO and 5JXE with report generation
2. RS (Residue Scanning) mode on all PDB files with --debug --not_relax
3. CM (Custom Mutations) mode on 7S7I.pdb with input CSV
4. DM (Double Mutation) mode on 7S7I.pdb with --debug --not_relax
5. Verification of output files and data integrity

Folder structure:
    tests/pdbs/                          # PDB files
        4ZFO.pdb, 5JXE.pdb, 7S7I.pdb
    tests/csv/                           # Reference and input CSVs
        manually_selected_4ZFO.csv
        manually_selected_5JXE.csv
        input_csv.csv
    tests/results/                       # Test output directory (created)
        selection/                       # Interface selection reports
        RS/{PDB}/                        # RS mode results
        custom_csv/                      # CM mode results
        double_mut/                      # DM mode results
"""

import json
import os
from pathlib import Path

import pytest

TESTS_DIR = Path(__file__).parent
PDBS_DIR = TESTS_DIR / "pdbs"
CSV_DIR = TESTS_DIR / "csv"
RESULTS_DIR = TESTS_DIR / "results"


# ============================================================================
# TEST 1: Interface Residue Selection
# ============================================================================


def find_pdb_files() -> dict[str, Path]:
    """Find all PDB files in pdbs folder."""
    pdbs = {}
    for pdb_file in sorted(PDBS_DIR.glob("*.pdb")):
        pdb_name = pdb_file.stem
        pdbs[pdb_name] = pdb_file
    return pdbs


def load_reference_positions(csv_file: Path) -> set[str]:
    """Load manually curated reference positions from CSV."""
    positions = set()
    if not csv_file.exists():
        return positions
    for line in csv_file.read_text().strip().splitlines():
        line = line.strip()
        if line:
            positions.add(line)
    return positions


def calculate_metrics(detected: set[str], expected: set[str]) -> dict:
    """Calculate precision, recall, F1 score."""
    tp = len(detected & expected)
    fp = len(detected - expected)
    fn = len(expected - detected)

    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1 = 2 * (precision * recall) / (precision + recall) if (precision + recall) > 0 else 0.0

    return {
        "true_positives": tp,
        "false_positives": fp,
        "false_negatives": fn,
        "precision": round(precision, 3),
        "recall": round(recall, 3),
        "f1_score": round(f1, 3),
        "detected": sorted(list(detected)),
        "missed": sorted(list(expected - detected)),
        "extra": sorted(list(detected - expected)),
    }


@pytest.mark.slow
def test_interface_selection_4zfo_5jxe():
    """Test interface residue selection on 4ZFO and 5JXE.

    Compares detected positions with manually_selected_*.csv files.
    Generates comprehensive report in tests/results/selection/.
    """
    import pyrosetta
    from scipy.spatial import cKDTree

    from residue_scanning.core import (
        PYROSETTA_FLAGS,
        _extract_ligand_coords,
        is_interface_residue,
    )

    pyrosetta.init(PYROSETTA_FLAGS)

    # Create output directory
    selection_dir = RESULTS_DIR / "selection"
    selection_dir.mkdir(parents=True, exist_ok=True)

    # Test both PDBs
    test_cases = {
        "4ZFO": ("4ZFO.pdb", "manually_selected_4ZFO.csv"),
        "5JXE": ("5JXE.pdb", "manually_selected_5JXE.csv"),
    }

    all_metrics = {}

    for pdb_name, (pdb_filename, ref_filename) in test_cases.items():
        pdb_file = PDBS_DIR / pdb_filename
        ref_file = CSV_DIR / ref_filename

        if not pdb_file.exists():
            pytest.skip(f"{pdb_file} not found")

        expected = load_reference_positions(ref_file)
        if not expected:
            pytest.skip(f"No reference data for {pdb_name}")

        # Load PDB
        pose = pyrosetta.pose_from_pdb(str(pdb_file))

        # Build chain map
        chains = {}
        for i in range(1, pose.num_chains() + 1):
            cs, ce = pose.chain_begin(i), pose.chain_end(i)
            name = pose.pdb_info().chain(cs)
            chains[name] = [i, cs, ce]

        # Detect ligand chain (usually 'A')
        ligand_chain = "A" if "A" in chains else min(chains.keys())

        # Build ligand KDTree
        ligand_resids = list(range(chains[ligand_chain][1], chains[ligand_chain][2] + 1))
        ligand_coords = _extract_ligand_coords(pose, ligand_resids)
        ligand_tree = cKDTree(ligand_coords)

        # Get receptor residues (all except ligand)
        receptor_chains = set(chains.keys()) - {ligand_chain}
        receptor_resids = []
        for chain in receptor_chains:
            cs, ce = chains[chain][1], chains[chain][2]
            receptor_resids.extend(range(cs, ce + 1))

        # Run interface detection
        detected = set()
        for resid in receptor_resids:
            if is_interface_residue(pose, resid, ligand_tree, cutoff=8.0):
                res = pose.residue(resid)
                pdb_pos = pose.pdb_info().pose2pdb(resid).split()
                detected.add(f"{res.name1()}{pdb_pos[1]}{pdb_pos[0]}")

        # Calculate metrics
        metrics = calculate_metrics(detected, expected)
        all_metrics[pdb_name] = metrics

        # Save metrics JSON
        pdb_dir = selection_dir / pdb_name
        pdb_dir.mkdir(parents=True, exist_ok=True)
        with open(pdb_dir / "metrics.json", "w") as f:
            json.dump(metrics, f, indent=2)

    # Generate comprehensive report
    report_path = selection_dir / "INTERFACE_SELECTION_REPORT.txt"
    with open(report_path, "w") as f:
        f.write("=" * 90 + "\n")
        f.write("INTERFACE RESIDUE SELECTION REPORT\n")
        f.write("=" * 90 + "\n\n")

        # Summary table
        f.write(
            f"{'PDB':<10} {'TP':>5} {'FP':>5} {'FN':>5} {'Precision':>10} {'Recall':>10} {'F1':>10}\n"
        )
        f.write("-" * 90 + "\n")

        for pdb_name in sorted(all_metrics.keys()):
            m = all_metrics[pdb_name]
            f.write(
                f"{pdb_name:<10} {m['true_positives']:>5} {m['false_positives']:>5} "
                f"{m['false_negatives']:>5} {m['precision']:>10.3f} {m['recall']:>10.3f} {m['f1_score']:>10.3f}\n"
            )

        f.write("\n" + "=" * 90 + "\n")
        f.write("DETAILED RESULTS\n")
        f.write("=" * 90 + "\n\n")

        for pdb_name in sorted(all_metrics.keys()):
            m = all_metrics[pdb_name]
            f.write(f"\n{pdb_name}\n")
            f.write("-" * 90 + "\n")
            f.write(f"Reference positions: {m['true_positives'] + m['false_negatives']}\n")
            f.write(f"Detected positions: {m['true_positives'] + m['false_positives']}\n")
            f.write(f"True Positives: {m['true_positives']}\n")
            f.write(f"False Positives: {m['false_positives']}\n")
            f.write(f"False Negatives: {m['false_negatives']}\n\n")
            f.write(f"Precision: {m['precision']:.3f}\n")
            f.write(f"Recall: {m['recall']:.3f}\n")
            f.write(f"F1 Score: {m['f1_score']:.3f}\n\n")

            if m["detected"]:
                f.write(f"Detected ({len(m['detected'])}): {', '.join(m['detected'][:10])}\n")
                if len(m["detected"]) > 10:
                    f.write(f"  ... and {len(m['detected']) - 10} more\n\n")

            if m["missed"]:
                f.write(f"❌ Missed ({len(m['missed'])}): {', '.join(m['missed'])}\n")

            if m["extra"]:
                f.write(f"⚠️  Extra ({len(m['extra'])}): {', '.join(m['extra'])}\n")

            if not m["missed"] and not m["extra"]:
                f.write("✅ Perfect match with reference\n")

    print(f"✅ Interface selection report: {report_path}")


# ============================================================================
# TEST 2: Residue Scanning (RS) Mode
# ============================================================================


@pytest.mark.slow
def test_residue_scanning_all_pdbs():
    """Run RS (Residue Scanning) mode on all PDB files.

    Uses --debug --not_relax --cpu 4 for quick testing.
    Outputs to tests/results/RS/{PDB}/.
    """
    from residue_scanning.cli import main

    # Create results directory
    rs_results_dir = RESULTS_DIR / "RS"
    rs_results_dir.mkdir(parents=True, exist_ok=True)

    original_cwd = os.getcwd()

    # Test each PDB file
    for pdb_file in sorted(PDBS_DIR.glob("*.pdb")):
        pdb_name = pdb_file.stem
        output_dir = rs_results_dir / pdb_name
        output_dir.mkdir(parents=True, exist_ok=True)

        try:
            print(f"\n{'='*60}")
            print(f"Testing RS mode on {pdb_name}")
            print(f"{'='*60}")

            # All test PDBs have chains: H, L (receptor), A (ligand)
            receptor_chains = "HL"
            ligand_chain = "A"

            main(
                argv=[
                    "-f",
                    str(pdb_file),
                    "-r",
                    receptor_chains,
                    "-l",
                    ligand_chain,
                    "-m",
                    "RS",
                    "--cpu",
                    "4",
                    "--debug",
                    "--not_relax",
                    "-o",
                    str(output_dir),
                ]
            )

            # Verify outputs
            assert (output_dir / "FastRelax.pdb").exists(), f"FastRelax.pdb missing for {pdb_name}"
            assert (
                output_dir / "Rosetta_ddG_mut.csv"
            ).exists(), f"Rosetta_ddG_mut.csv missing for {pdb_name}"

            print(f"✅ RS mode test passed for {pdb_name}")

        finally:
            os.chdir(original_cwd)


# ============================================================================
# TEST 3: Custom Mutations (CM) Mode
# ============================================================================


@pytest.mark.slow
def test_custom_mutations_7s7i():
    """Run CM (Custom Mutations) mode on 7S7I.pdb.

    Uses tests/csv/input_csv.csv as --mutations_csv.
    Outputs to tests/results/custom_csv/.
    """
    from residue_scanning.cli import main

    pdb_file = PDBS_DIR / "7S7I.pdb"
    csv_file = CSV_DIR / "input_csv.csv"

    if not pdb_file.exists():
        pytest.skip(f"{pdb_file} not found")
    if not csv_file.exists():
        pytest.skip(f"{csv_file} not found")

    output_dir = RESULTS_DIR / "custom_csv"
    output_dir.mkdir(parents=True, exist_ok=True)

    original_cwd = os.getcwd()

    try:
        print(f"\n{'='*60}")
        print("Testing CM mode on 7S7I.pdb")
        print(f"{'='*60}")

        main(
            argv=[
                "-f",
                str(pdb_file),
                "-r",
                "HL",
                "-l",
                "A",
                "-m",
                "CM",
                "--mutations_csv",
                str(csv_file),
                "--cpu",
                "4",
                "--not_relax",
                "-o",
                str(output_dir),
            ]
        )

        # Verify outputs
        assert (output_dir / "FastRelax.pdb").exists(), "FastRelax.pdb missing for CM mode"
        assert (
            output_dir / "Rosetta_ddG_mut.csv"
        ).exists(), "Rosetta_ddG_mut.csv missing for CM mode"

        print("✅ CM mode test passed for 7S7I.pdb")

    finally:
        os.chdir(original_cwd)


# ============================================================================
# TEST 4: Double Mutations (DM Mode)
# ============================================================================


@pytest.mark.slow
def test_double_mutations_7s7i_short():
    """Quick test for DM (Double Mutation) mode on 7S7I.pdb.

    Uses --debug to limit to 4 interface positions (~100-200 double mutations).
    Runs both condition-based and CSV-based DM modes.

    Configuration: --debug --not_relax --cpu 8

    Duration: ~3-5 minutes
    Outputs to tests/results/double_mut_short/
    """
    from residue_scanning.cli import main

    pdb_file = PDBS_DIR / "7S7I.pdb"

    if not pdb_file.exists():
        pytest.skip(f"{pdb_file} not found")

    output_dir = RESULTS_DIR / "double_mut_short"
    output_dir.mkdir(parents=True, exist_ok=True)

    original_cwd = os.getcwd()

    try:
        # Sub-test 1: Condition-based DM mode (short)
        print(f"\n{'='*60}")
        print("Testing DM mode (condition-based SHORT) on 7S7I.pdb")
        print(f"{'='*60}")

        output_dir_condition = output_dir / "condition_based"
        output_dir_condition.mkdir(parents=True, exist_ok=True)

        main(
            argv=[
                "-f",
                str(pdb_file),
                "-r",
                "HL",
                "-l",
                "A",
                "-m",
                "DM",
                "--condition",
                "ddG_interface < 0",
                "--cpu",
                "8",
                "--not_relax",
                "--debug",
                "-o",
                str(output_dir_condition),
            ]
        )

        # Verify outputs
        assert (
            output_dir_condition / "FastRelax.pdb"
        ).exists(), "FastRelax.pdb missing for condition-based DM"
        assert (
            output_dir_condition / "Rosetta_ddG_mut_double.csv"
        ).exists(), "Rosetta_ddG_mut_double.csv missing for condition-based DM"

        # Verify double mutation CSV has expected columns
        import pandas as pd

        df_double = pd.read_csv(output_dir_condition / "Rosetta_ddG_mut_double.csv")
        assert not df_double.empty, "Double mutation CSV is empty"
        assert "Name" in df_double.columns, "Missing 'Name' column in double mutations"
        assert (
            "ddG_interface" in df_double.columns
        ), "Missing 'ddG_interface' column in double mutations"

        # Verify double mutation names contain underscore (format: XXX_YYY)
        if len(df_double) > 0:
            for name in df_double["Name"].head():
                assert (
                    "_" in name
                ), f"Double mutation name should contain underscore: {name}"

        print(
            f"✅ Condition-based DM SHORT test passed ({len(df_double)} double mutations)"
        )

        # Sub-test 2: CSV-based DM mode (short)
        print(f"\n{'='*60}")
        print("Testing DM mode (CSV-based SHORT) on 7S7I.pdb")
        print(f"{'='*60}")

        csv_file = CSV_DIR / "input_csv.csv"
        if not csv_file.exists():
            pytest.skip(f"{csv_file} not found")

        output_dir_csv = output_dir / "csv_based"
        output_dir_csv.mkdir(parents=True, exist_ok=True)

        main(
            argv=[
                "-f",
                str(pdb_file),
                "-r",
                "HL",
                "-l",
                "A",
                "-m",
                "DM",
                "--mutations_csv",
                str(csv_file),
                "--cpu",
                "8",
                "--not_relax",
                "--debug",
                "-o",
                str(output_dir_csv),
            ]
        )

        # Verify outputs
        assert (
            output_dir_csv / "FastRelax.pdb"
        ).exists(), "FastRelax.pdb missing for CSV-based DM"
        assert (
            output_dir_csv / "Rosetta_ddG_mut_double.csv"
        ).exists(), "Rosetta_ddG_mut_double.csv missing for CSV-based DM"

        # Verify CSV-based DM results
        df_csv = pd.read_csv(output_dir_csv / "Rosetta_ddG_mut_double.csv")
        assert not df_csv.empty, "CSV-based DM result CSV is empty"
        assert "Name" in df_csv.columns, "Missing 'Name' column in CSV-based DM results"

        # Verify double mutation names contain underscore
        if len(df_csv) > 0:
            for name in df_csv["Name"]:
                assert (
                    "_" in name
                ), f"Double mutation name should contain underscore: {name}"

        print(f"✅ CSV-based DM SHORT test passed ({len(df_csv)} double mutations)")

    finally:
        os.chdir(original_cwd)


@pytest.mark.slow
def test_double_mutations_7s7i_long(cpu: str = "4", debug: bool = False, not_relax: bool = True):
    """Comprehensive test for DM (Double Mutation) mode on 7S7I.pdb.

    Parameters:
    - cpu: Number of CPUs to use (default: "4")
    - debug: Limit to 4 interface positions if True (default: False - evaluate all)
    - not_relax: Skip FastRelax if True (default: True)

    Runs both condition-based and CSV-based DM modes.

    Duration:
    - With --debug: ~5 minutes
    - Without --debug (all mutations): ~30-45 minutes

    Outputs to tests/results/double_mut_long/
    """
    from residue_scanning.cli import main

    pdb_file = PDBS_DIR / "7S7I.pdb"

    if not pdb_file.exists():
        pytest.skip(f"{pdb_file} not found")

    output_dir = RESULTS_DIR / "double_mut_long"
    output_dir.mkdir(parents=True, exist_ok=True)

    original_cwd = os.getcwd()

    try:
        # Sub-test 1: Condition-based DM mode (long)
        print(f"\n{'='*60}")
        print("Testing DM mode (condition-based LONG) on 7S7I.pdb")
        print(f"  CPU cores: {cpu}")
        print(f"  Debug mode: {debug}")
        print(f"  Skip FastRelax: {not_relax}")
        print(f"{'='*60}")

        output_dir_condition = output_dir / "condition_based"
        output_dir_condition.mkdir(parents=True, exist_ok=True)

        argv = [
            "-f",
            str(pdb_file),
            "-r",
            "HL",
            "-l",
            "A",
            "-m",
            "DM",
            "--condition",
            "ddG_interface < 0",
            "--cpu",
            cpu,
        ]
        if not_relax:
            argv.append("--not_relax")
        if debug:
            argv.append("--debug")
        argv.extend(["-o", str(output_dir_condition)])

        main(argv=argv)

        # Verify outputs
        assert (
            output_dir_condition / "FastRelax.pdb"
        ).exists(), "FastRelax.pdb missing for condition-based DM"
        assert (
            output_dir_condition / "Rosetta_ddG_mut_double.csv"
        ).exists(), "Rosetta_ddG_mut_double.csv missing for condition-based DM"

        # Verify double mutation CSV has expected columns
        import pandas as pd

        df_double = pd.read_csv(output_dir_condition / "Rosetta_ddG_mut_double.csv")
        assert not df_double.empty, "Double mutation CSV is empty"
        assert "Name" in df_double.columns, "Missing 'Name' column in double mutations"
        assert (
            "ddG_interface" in df_double.columns
        ), "Missing 'ddG_interface' column in double mutations"

        # Verify double mutation names contain underscore (format: XXX_YYY)
        if len(df_double) > 0:
            for name in df_double["Name"].head():
                assert (
                    "_" in name
                ), f"Double mutation name should contain underscore: {name}"

        print(
            f"✅ Condition-based DM LONG test passed ({len(df_double)} double mutations)"
        )

        # Sub-test 2: CSV-based DM mode (long)
        print(f"\n{'='*60}")
        print("Testing DM mode (CSV-based LONG) on 7S7I.pdb")
        print(f"  CPU cores: {cpu}")
        print(f"  Debug mode: {debug}")
        print(f"  Skip FastRelax: {not_relax}")
        print(f"{'='*60}")

        csv_file = CSV_DIR / "input_csv.csv"
        if not csv_file.exists():
            pytest.skip(f"{csv_file} not found")

        output_dir_csv = output_dir / "csv_based"
        output_dir_csv.mkdir(parents=True, exist_ok=True)

        argv = [
            "-f",
            str(pdb_file),
            "-r",
            "HL",
            "-l",
            "A",
            "-m",
            "DM",
            "--mutations_csv",
            str(csv_file),
            "--cpu",
            cpu,
        ]
        if not_relax:
            argv.append("--not_relax")
        if debug:
            argv.append("--debug")
        argv.extend(["-o", str(output_dir_csv)])

        main(argv=argv)

        # Verify outputs
        assert (
            output_dir_csv / "FastRelax.pdb"
        ).exists(), "FastRelax.pdb missing for CSV-based DM"
        assert (
            output_dir_csv / "Rosetta_ddG_mut_double.csv"
        ).exists(), "Rosetta_ddG_mut_double.csv missing for CSV-based DM"

        # Verify CSV-based DM results
        df_csv = pd.read_csv(output_dir_csv / "Rosetta_ddG_mut_double.csv")
        assert not df_csv.empty, "CSV-based DM result CSV is empty"
        assert "Name" in df_csv.columns, "Missing 'Name' column in CSV-based DM results"

        # Verify double mutation names contain underscore
        if len(df_csv) > 0:
            for name in df_csv["Name"]:
                assert (
                    "_" in name
                ), f"Double mutation name should contain underscore: {name}"

        print(f"✅ CSV-based DM LONG test passed ({len(df_csv)} double mutations)")

    finally:
        os.chdir(original_cwd)


# ============================================================================
# TEST 5: Verify Output Files
# ============================================================================


@pytest.mark.slow
def test_verify_all_outputs():
    """Verify that all output files are created and valid."""
    import pandas as pd

    # Check RS results
    rs_dir = RESULTS_DIR / "RS"
    if rs_dir.exists():
        for pdb_dir in rs_dir.iterdir():
            if pdb_dir.is_dir():
                csv_file = pdb_dir / "Rosetta_ddG_mut.csv"
                if csv_file.exists():
                    df = pd.read_csv(csv_file)
                    assert not df.empty, f"{csv_file} is empty"
                    assert "Name" in df.columns, f"Missing 'Name' column in {csv_file}"
                    assert (
                        "ddG_complex" in df.columns
                    ), f"Missing 'ddG_complex' column in {csv_file}"
                    print(f"✅ Verified {pdb_dir.name} RS results ({len(df)} rows)")

    # Check CM results
    cm_dir = RESULTS_DIR / "custom_csv"
    if cm_dir.exists():
        csv_file = cm_dir / "Rosetta_ddG_mut.csv"
        if csv_file.exists():
            df = pd.read_csv(csv_file)
            assert not df.empty, "CM result CSV is empty"
            assert "Name" in df.columns, "Missing 'Name' column in CM results"
            print(f"✅ Verified CM results ({len(df)} rows)")

    # Check DM (double mutation) results
    dm_dir = RESULTS_DIR / "double_mut"
    if dm_dir.exists():
        csv_file = dm_dir / "Rosetta_ddG_mut_double.csv"
        if csv_file.exists():
            df = pd.read_csv(csv_file)
            assert not df.empty, "DM result CSV is empty"
            assert "Name" in df.columns, "Missing 'Name' column in DM results"
            # Verify double mutation format (contains underscore)
            for name in df["Name"]:
                assert "_" in name, f"Invalid double mutation name format: {name}"
            print(f"✅ Verified DM results ({len(df)} double mutations)")

    # Check interface selection report
    report_file = RESULTS_DIR / "selection" / "INTERFACE_SELECTION_REPORT.txt"
    if report_file.exists():
        content = report_file.read_text()
        assert "4ZFO" in content, "Missing 4ZFO in interface selection report"
        assert "5JXE" in content, "Missing 5JXE in interface selection report"
        print("✅ Verified interface selection report")

    print("\n✅ All output files verified successfully")
