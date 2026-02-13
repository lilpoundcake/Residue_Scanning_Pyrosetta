# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Environment

Conda/Micromamba environment name: `rosetta_RS`. PyRosetta is only available from the RosettaCommons conda channel (credentials embedded in `env.yaml`).

```bash
micromamba env create -f env.yaml
micromamba activate rosetta_RS
pip install -e .          # installs the residue-scan CLI entry point
```

## Running the CLI

```bash
residue-scan -f complex.pdb -o results -r HL -l A -m RS --cpu 8
```

Modes (short/long aliases):
- `FR` / `FastRelax` — minimize input structure only
- `RS` / `Residue_Scanning` — compute ΔΔG for all single-point mutations at the interface
- `DM` / `Double_Mut_Searching` — RS + filter by `--condition`, then evaluate neighboring double mutations
- `CM` / `Custom_Mutations` — compute ΔΔG for user-supplied mutations/positions from a CSV file

Use `--debug` to limit the mutation set to 4 positions for quick testing. Use `--not_relax` to skip FastRelax (faster but less accurate; switches scoring function from `ref2015_cart` to `ref2015`).

### Custom Mutations mode (CM)

```bash
residue-scan -f complex.pdb -o results -r HL -l A -m CM --mutations_csv muts.csv --cpu 8
```

The CSV is a single-column file (no header required) with one entry per row:
- **Position format** `AH98` — generates all 18 standard mutations (excludes Cys/Pro)
- **Mutation format** `AH98D` — specific single mutation (from=A, chain=H, number=98, to=D)
- **With insertion code** `AH98bD` (mutation) or `AH98b` (position)

Validation: `from_aa` is checked against the actual residue in the pose; mismatches are warned and skipped. Wildtype entries are auto-generated for ddG calculation.

## Running tests

```bash
pytest tests/ -v -m slow -s   # all comprehensive tests
```

Test structure:
- `tests/pdbs/` — PDB files (4ZFO, 5JXE, 7S7I)
- `tests/csv/` — Reference and input CSVs
- `tests/results/` — Test output (created by tests)
- `tests/test_comprehensive.py` — Main test suite

See `tests/README.md` for detailed test documentation.

## Linting

```bash
ruff check residue_scanning/ tests/
flake8 residue_scanning/ tests/
mypy residue_scanning/
black --check residue_scanning/ tests/
```

## Building a wheel

```bash
python -m build    # produces dist/residue_scanning-*.whl
```

## Package architecture

```
residue_scanning/
    __init__.py          version string
    cli.py               argparse entry point — main() + _parse_args()
    preprocessing.py     normalize_residue_names(), prepare_pdb()
    core.py              all PyRosetta computation functions
tests/
    5JXE.pdb             antibody-antigen complex (chains H,L=receptor, A=ligand)
    manually_selected.csv 26 curated interface positions (reference for filter calibration)
    test_debug.py        pytest integration tests
pyproject.toml           setuptools build config + entry point
env.yaml                 micromamba environment
.flake8                  flake8 configuration
```

**Module responsibilities:**

- `preprocessing.py` — Python-only PDB text processing. `normalize_residue_names()` handles non-standard protonated residue names (HID/HIE/HIP/HSP/HSD → HIS, ASPH/ASPP/ASH → ASP, GLUH/GLUP/GLH → GLU, LYSH/LYP → LYS, CYX/CYM → CYS, etc.) by rewriting the correct PDB columns for ATOM/HETATM lines. `prepare_pdb()` calls `pyrosetta.toolbox.cleanATOM` after normalization.

- `core.py` — All PyRosetta functions. Imports `from pyrosetta import *` at module level (safe before `pyrosetta.init()`). Key functions: `prepare_mut_df()`, `prepare_custom_mut_df()`, `df_FastRelax()`, `FastRelax_replics()`, `mutate_pose_FR()`, `df_ddG_postprocessing()`, `is_interface_residue()`, `_extract_ligand_coords()`.

- `cli.py` — Calls `pyrosetta.init()` inside `main()` (not at module level), then imports from `core` and `preprocessing`. Orchestrates the FR → RS → DM / CM pipeline sequentially.

**Execution flow:**
1. **PDB preprocessing** — `normalize_residue_names()` rewrites protonated residue names in ATOM/HETATM columns, then `cleanATOM` strips non-ATOM records → `clear.pdb`
2. **FastRelax block** — `REPLICS` parallel workers via `multiprocessing.Pool`, picks lowest-REU pose → `FastRelax.pdb`
3. **Residue Scanning block** — `prepare_mut_df()` builds a `cKDTree` from ligand atom coordinates, then runs `is_interface_residue()` on every receptor residue to identify the interface. Generates all mutations (18 AA: excludes Cys and Pro). Parallel `pool.starmap(df_FastRelax, ...)` → `Rosetta_ddG_mut.csv`
4. **Custom Mutations block** — `prepare_custom_mut_df()` reads a user CSV, parses mutations/positions, validates from_aa, auto-generates wildtype entries → same parallel pipeline as RS
5. **Double Mutation block** — filters single mutations by `--condition`, finds mutual spatial neighbours, generates combinations → `Rosetta_ddG_mut_double.csv`

**Interface filter** (`is_interface_residue`): two-tier hybrid filter using `scipy.spatial.cKDTree` for fast nearest-neighbour queries against pre-computed ligand atom coordinates. Tier 1: SC atom within 4.0 Å = unconditional accept. Tier 2: SC atom within 5.0 Å = accept if CB oriented toward ligand (angle < 90°) AND backbone CA within 9.0 Å. Glycine: CA within 4.0 Å. Calibrated against `tests/manually_selected.csv` (26 positions, zero false positives).

**Mutation name format:** `{from_aa}{chain}{pdb_number}{insertion_code?}{to_aa}` — e.g., `GH122AL`. Double mutations use `_` separator. Wildtype entries have the last character equal to the first.

**ddG calculation:** subtracts the matched wildtype row's `REU_min` and `Min_dG_Interface`. Wildtype row identified by `wildtyper()`.

**Early convergence:** `FastRelax_replics` stops after replica 2 if `|score[0] − score[1]| < 1` REU (skipped when `--not_relax`).

**Multiprocessing:** All parameters are passed explicitly to worker functions (no module-level globals) so they are picklable by `pool.starmap`.
