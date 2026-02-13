# Residue Scanning PyRosetta

Protein-protein interface residue scanning and ΔΔG prediction using PyRosetta.

---

## Installation

> **PyRosetta requirement**: PyRosetta is only available from the RosettaCommons conda
> channel and **cannot** be installed with pip or poetry alone. All methods below
> require `micromamba` or `conda` to already be installed on your system.

### Method 1: micromamba / conda (recommended)

This is the only method that installs PyRosetta automatically.
Works on macOS (Intel and Apple Silicon) and Linux. Requires Python 3.10+.

```bash
git clone https://github.com/lilpoundcake/Residue_Scanning_Pyrosetta.git
cd Residue_Scanning_Pyrosetta

# Create the environment (Python 3.10, PyRosetta, all dependencies)
micromamba env create -f env.yaml   # or: conda env create -f env.yaml
micromamba activate rosetta_RS      # or: conda activate rosetta_RS

# Install the residue-scanning package
pip install -e .                    # editable / development install
```

After this, `residue-scan` is available everywhere inside the activated environment.

**What `env.yaml` provides:**

| Package | Source |
|---------|--------|
| Python 3.10 | conda-forge |
| PyRosetta | RosettaCommons channel (credentials in `env.yaml`) |
| numpy, pandas, scipy, biopython, python-blosc | conda-forge |

---

### Method 2: pip (editable install for development)

Use this after the conda environment already exists (e.g. from Method 1).

```bash
micromamba activate rosetta_RS   # activate the environment with PyRosetta

# Editable install — code changes take effect immediately
pip install -e .

# Verify
residue-scan --help
```

---

### Method 3: install from a pre-built wheel

Use this to install a specific released version without a source checkout.

```bash
micromamba activate rosetta_RS   # PyRosetta must already be present

pip install residue_scanning-0.2.1-py3-none-any.whl
```

**Build the wheel yourself from source:**

```bash
micromamba activate rosetta_RS
pip install build
python -m build              # produces dist/residue_scanning-*.whl and .tar.gz
pip install dist/*.whl
```

---

### Method 4: poetry (within the conda environment)

Poetry cannot install PyRosetta (conda-only package), so the conda environment
must be created first. Poetry is then used to manage the remaining Python
dependencies from `pyproject.toml`.

```bash
# 1. Create the conda environment (includes PyRosetta)
micromamba env create -f env.yaml
micromamba activate rosetta_RS

# 2. Install poetry into the activated environment
pip install poetry

# 3. Install the package and dev dependencies via poetry
poetry install

# 4. Run commands through poetry
poetry run residue-scan --help
```

> **Note:** `poetry install` reads `pyproject.toml`. Because PyRosetta is not on PyPI,
> it must remain in the conda environment — poetry will not manage it.

---

## Usage

```bash
residue-scan -f complex.pdb -o results -r HL -l A -m RS --cpu 8
```

### Required arguments

| Flag | Description |
|------|-------------|
| `-f` / `--file` | Path to the input PDB file |
| `-m` / `--mode` | Analysis mode (see table below) |
| `-r` / `--receptor` | Receptor chain IDs, concatenated (e.g. `HL` for chains H and L) |
| `-l` / `--ligand` | Ligand chain IDs, concatenated (e.g. `A`) |
| `--cpu` | Number of CPU cores |

### Modes

| Short | Long | Description |
|-------|------|-------------|
| `FR` | `FastRelax` | Minimize input structure only |
| `RS` | `Residue_Scanning` | Compute ΔΔG for all single-point mutations at the interface |
| `DM` | `Double_Mut_Searching` | RS + evaluate double-mutation combinations filtered by `--condition` |
| `CM` | `Custom_Mutations` | Compute ΔΔG for user-supplied mutations/positions from a CSV file |

### Optional arguments

| Flag | Default | Description |
|------|---------|-------------|
| `-o` / `--output` | `results` | Output directory |
| `--not_relax` | off | Skip FastRelax (uses PackRotamers; faster but less accurate) |
| `--replics` | `5` | Max stochastic replicas per mutation |
| `--radius` | `8.0` | Neighbor shell radius in Å for repacking |
| `--condition` | `ddG_complex < 0.5 and ddG_interface < 0.5` | Pandas query to select single mutations in DM mode |
| `--debug` | off | Limit to 4 interface positions (quick smoke test) |
| `--mutations_csv` | — | Path to CSV with custom mutations/positions (required for CM mode) |

---

## Examples

Antibody-antigen residue scanning (Fv chains H+L, antigen chain A):

```bash
residue-scan -f 5JXE.pdb -o RS_results -r HL -l A -m RS --cpu 110
```

Double-mutation search with custom condition:

```bash
residue-scan -f complex.pdb -o results -r HL -l A \
    -m DM --cpu 110 \
    --radius 10 --condition "ddG_complex < 0 and ddG_interface < 0"
```

Quick run without minimization:

```bash
residue-scan -f complex.pdb -o results -r HL -l A -m RS --cpu 10 --not_relax
```

Custom mutations from CSV:

```bash
residue-scan -f 5JXE.pdb -o CM_results -r HL -l A -m CM --mutations_csv muts.csv --cpu 8
```

The CSV is a single-column file (no header required) with one entry per row:
- **Position format** `AH98` — generates all 18 standard mutations (excludes Cys/Pro)
- **Mutation format** `AH98D` — specific single mutation (from=A, chain=H, number=98, to=D)
- **With insertion code** `AH98AD` (mutation at H98A) or `AH98A` (position H98A)

**⚠️ Important:** Insertion codes use UPPERCASE letters (A, B, C) to match PDB format.

See `CSV_FORMAT_GUIDE.md` for complete CSV format documentation with examples.

Debug / smoke test:

```bash
residue-scan -f tests/5JXE.pdb -r HL -l A -m RS --cpu 2 --debug --not_relax -o tests/results
```

---

## Output files

All files are written to the output directory:

| File / folder | Description |
|---------------|-------------|
| `FastRelax.pdb` | Lowest-REU relaxed structure |
| `Rosetta_ddG_mut.csv` | ΔΔG table for all single mutations (`ddG_complex`, `ddG_interface`) |
| `Rosetta_ddG_mut_double.csv` | ΔΔG table for double mutations (DM mode only) |
| `min_REU_pdb/` | Lowest-energy structures for each single mutation |
| `min_REU_pdb_double/` | Lowest-energy structures for each double mutation |

---

## Running tests

```bash
micromamba activate rosetta_RS
pip install -e .

pytest tests/ -v -m slow -s
```

The comprehensive test suite includes:
1. **Interface Selection** — Validates KDTree filter on 4ZFO and 5JXE
2. **Residue Scanning** — RS mode on all available PDBs
3. **Custom Mutations** — CM mode on 7S7I with user CSV
4. **Output Verification** — Ensures all results are valid

See `tests/README.md` for details and expected results.

---

## Linting

```bash
ruff check residue_scanning/ tests/
flake8 residue_scanning/ tests/
mypy residue_scanning/
black --check residue_scanning/ tests/
```

---

## Input PDB preprocessing

The pipeline automatically handles non-standard protonated residue names before
passing the file to PyRosetta:

| Non-standard | Canonical |
|---|---|
| HID, HIE, HIP, HSD, HSE, HSP, HISB | HIS |
| ASH, ASPH, ASPP | ASP |
| GLH, GLUH, GLUP | GLU |
| LYSH, LYP, LYN | LYS |
| ARGN | ARG |
| CYSH, CYX, CYM | CYS |

---

## Residue selection logic

Interface residues are identified using a two-tier hybrid filter (`is_interface_residue`) with a `scipy.spatial.cKDTree` built from all ligand atom coordinates:

1. **Tier 1 — direct contact** (SC ≤ 4.0 Å): any sidechain atom within 4.0 Å of the nearest ligand atom. Accepted unconditionally — proven interaction.
2. **Tier 2 — oriented contact** (4.0 < SC ≤ 5.0 Å): sidechain atom within 5.0 Å, accepted only if CB is oriented toward the ligand (CA→CB vs CA→nearest-ligand angle < 90°) AND the backbone CA is within 9.0 Å of the ligand. Eliminates residues pointing away or reaching across from distant backbone.
3. **Glycine proxy** — CA within 4.0 Å of any ligand atom (because mutating FROM glycine adds a sidechain).

The default thresholds capture all 26 manually curated interface positions on the 5JXE test case with zero false positives.

---

## Citation

Protocol inspired by:
**"Prediction of Protein Mutational Free Energy: Benchmark and Sampling Improvements Increase Classification Accuracy"**
DOI: [10.3389/fbioe.2020.558247](https://www.frontiersin.org/journals/bioengineering-and-biotechnology/articles/10.3389/fbioe.2020.558247/full)
