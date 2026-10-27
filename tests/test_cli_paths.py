"""CLI path resolution tests that do not require PyRosetta."""

from pathlib import Path

import pytest

from residue_scanning.cli import _parse_args, _resolve_paths


@pytest.mark.parametrize("output", [None, "results", "./results", "nested/results", "../results"])
@pytest.mark.parametrize("absolute_input", [False, True])
def test_relative_paths_use_launch_directory(tmp_path, monkeypatch, output, absolute_input):
    launch_dir = tmp_path / "launch"
    launch_dir.mkdir()
    monkeypatch.chdir(launch_dir)
    pdb = "data/complex.pdb"
    if absolute_input:
        pdb = str(launch_dir / pdb)
    argv = [
        "-f", pdb, "-m", "CM", "-r", "HL", "-l", "A", "--cpu", "1",
        "--mutations_csv", "inputs/mutations.csv",
    ]
    if output is not None:
        argv.extend(["-o", output])

    input_file, mutations_csv, output_dir = _resolve_paths(_parse_args(argv))

    assert input_file == launch_dir / "data/complex.pdb"
    assert mutations_csv == launch_dir / "inputs/mutations.csv"
    assert output_dir == (launch_dir / (output or "results")).resolve()
    assert output_dir != (input_file.parent / (output or "results")).resolve()

    # Once resolved, all paths remain usable after the pipeline changes directory.
    output_dir.mkdir(parents=True, exist_ok=True)
    monkeypatch.chdir(output_dir)
    assert input_file == launch_dir / "data/complex.pdb"
    assert mutations_csv == launch_dir / "inputs/mutations.csv"
    assert output_dir == Path.cwd()


def test_absolute_paths_are_preserved(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    input_file = tmp_path / "structures/complex.pdb"
    mutations_csv = tmp_path / "inputs/mutations.csv"
    output_dir = tmp_path / "output/results"
    args = _parse_args([
        "-f", str(input_file), "-m", "CM", "-r", "HL", "-l", "A", "--cpu", "1",
        "--mutations_csv", str(mutations_csv), "-o", str(output_dir),
    ])

    assert _resolve_paths(args) == (input_file, mutations_csv, output_dir)


def test_optional_csv_is_none(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    args = _parse_args([
        "-f", "data/complex.pdb", "-m", "FR", "-r", "HL", "-l", "A", "--cpu", "1",
    ])

    assert _resolve_paths(args) == (tmp_path / "data/complex.pdb", None, tmp_path / "results")
