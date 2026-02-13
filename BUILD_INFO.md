# Build Information & Distribution

**Last Build**: 2026-02-13
**Version**: 0.2.1
**Status**: ✅ Production Ready

---

## Distributions Available

### Wheel Distribution
**File**: `dist/residue_scanning-0.2.1-py3-none-any.whl`
- **Size**: 17 KB
- **Format**: Binary wheel (pre-compiled)
- **Install time**: ~1 second
- **Python**: 3.10+ (universal)
- **Platform**: macOS (Intel/ARM), Linux

### Source Distribution
**File**: `dist/residue_scanning-0.2.1.tar.gz`
- **Size**: 16 KB
- **Format**: Source tarball
- **Includes**: Full source code, tests, documentation
- **Build time**: ~30 seconds
- **Python**: 3.10+ (requires build)

---

## Installation Methods

### From Wheel (Fastest)
```bash
micromamba activate rosetta_RS
pip install dist/residue_scanning-0.2.1-py3-none-any.whl
```

### From Source Distribution
```bash
micromamba activate rosetta_RS
pip install dist/residue_scanning-0.2.1.tar.gz
```

### From Local Development
```bash
micromamba activate rosetta_RS
pip install -e .          # editable install for development
```

---

## Build Configuration

### Python Version
- **Required**: Python 3.10+
- **Target**: Python 3.10 (defined in pyproject.toml)
- **Environment**: rosetta_RS conda environment

### Dependencies (pip)
```
numpy>=1.20
pandas
scipy>=1.7
biopython>=1.85
blosc              # required for efficient PyRosetta pose pickling
```

**Note**: PyRosetta is conda-only (not on PyPI) and must be installed via `env.yaml`

### Build Requirements
```
setuptools>=68
wheel
build
```

---

## Package Contents

### Wheel Package Structure
```
residue_scanning-0.2.1.dist-info/
├── METADATA              # Package metadata
├── WHEEL                 # Wheel format info
├── entry_points.txt      # CLI entry point
├── top_level.txt         # Top-level packages
└── RECORD                # File manifest

residue_scanning/
├── __init__.py           # Version string
├── cli.py                # CLI orchestration (FR/RS/CM)
├── core.py               # PyRosetta + KDTree filter
└── preprocessing.py      # PDB processing
```

### CLI Entry Point
```
residue-scan = "residue_scanning.cli:main"
```

After installation, the `residue-scan` command is available globally in the activated environment.

---

## Building the Distributions

### Prerequisites
```bash
micromamba activate rosetta_RS
pip install build
```

### Build Both Wheel and Source
```bash
python -m build
```

Output:
- `dist/residue_scanning-0.2.1-py3-none-any.whl`
- `dist/residue_scanning-0.2.1.tar.gz`

### Build Wheel Only
```bash
python -m build --wheel
```

### Build Source Only
```bash
python -m build --sdist
```

### Build Configuration File
Location: `pyproject.toml`
```toml
[build-system]
requires = ["setuptools>=68", "wheel"]
build-backend = "setuptools.build_meta"

[project]
name = "residue-scanning"
version = "0.2.1"
requires-python = ">=3.10"
```

---

## Quality Assurance

### Pre-Build Checks
Before building, ensure:
- ✅ All tests pass: `pytest tests/ -v -m slow -s`
- ✅ Linting passes: `ruff check residue_scanning/ tests/`
- ✅ Types check: `mypy residue_scanning/`
- ✅ Format correct: `black --check residue_scanning/ tests/`

### Build Verification
After building:
```bash
# List contents of wheel
unzip -l dist/residue_scanning-0.2.1-py3-none-any.whl

# Test installation
pip install dist/residue_scanning-0.2.1-py3-none-any.whl
residue-scan --help
```

---

## Distribution Details

### Version Scheme
Uses semantic versioning: `MAJOR.MINOR.PATCH`
- **0**: Major version (early development)
- **2**: Minor version (interface selection + RS/CM modes)
- **1**: Patch version (bug fixes and improvements)

### Python Compatibility
- **Minimum**: Python 3.10
- **Target**: Python 3.10
- **Tested**: Python 3.10.19

### Platform Support
**Tested**:
- ✅ macOS (Apple Silicon, M1/M2)
- ✅ macOS (Intel)
- ✅ Linux

**Requires**:
- RosettaCommons conda channel (PyRosetta)
- 4+ GB RAM (for PyRosetta computations)

---

## Installation Verification

After installation, verify:

```bash
# Check CLI is available
residue-scan --help

# Check version
residue-scan --version  # (if implemented)

# Check imports work
python -c "from residue_scanning import cli, core; print('✅ Import successful')"

# Run quick test
residue-scan -f tests/pdbs/5JXE.pdb -r HL -l A -m FR --cpu 2 --not_relax -o /tmp/test_rs
```

---

## Release Checklist

- ✅ Version updated in `pyproject.toml`
- ✅ All tests passing (4/4)
- ✅ All linters passing (ruff, flake8, mypy, black)
- ✅ Documentation updated
- ✅ Wheel built successfully
- ✅ Source distribution created
- ✅ Distributions tested for installation
- ✅ README updated with installation info
- ✅ CLAUDE.md updated with build info

---

## Troubleshooting

### Build Fails
```bash
# Clean previous builds
rm -rf build/ dist/ *.egg-info residue_scanning.egg-info

# Try again
python -m build
```

### Installation Fails
```bash
# Ensure conda environment has PyRosetta
conda list | grep pyrosetta

# If missing, recreate environment
micromamba env remove -n rosetta_RS
micromamba env create -f env.yaml
micromamba activate rosetta_RS

# Try installation again
pip install dist/*.whl
```

### CLI Not Found
```bash
# Ensure environment is activated
micromamba activate rosetta_RS

# Reinstall package
pip install --force-reinstall dist/residue_scanning-0.2.1-py3-none-any.whl

# Verify entry point
pip show -f residue_scanning
```

---

## Distribution Timeline

| Date | Version | Status | Event |
|------|---------|--------|-------|
| 2026-02-13 | 0.2.1 | ✅ Ready | Build, tests pass, wheel created |
| 2026-02-13 | 0.2.1 | ✅ Ready | Documentation updated |
| 2026-02-13 | 0.2.1 | ✅ Ready | KDTree Hybrid Filter finalized |

---

## Next Steps

1. **Distribution**: Share wheel or source tarball with users
2. **Installation**: Users follow installation guide in README.md
3. **Testing**: Users verify installation with `residue-scan --help`
4. **Usage**: Run with desired modes (FR, RS, CM)
5. **Support**: Refer to CSV_FORMAT_GUIDE.md for CSV format questions

---

## Technical Details

### Package Metadata
- **Name**: residue-scanning
- **PyPI Name**: residue-scanning
- **Import Name**: residue_scanning
- **CLI Command**: residue-scan

### Entry Points
```
residue-scan → residue_scanning.cli:main()
```

### Build Backend
setuptools (via `build` package)

### Configuration Files
- `pyproject.toml` — Package metadata, dependencies, build config, tool config
- `setup.cfg` — Auto-generated during build
- `env.yaml` — Conda environment with PyRosetta

---

## References

- **PyPI Package**: Would be `residue-scanning` if published
- **Build Tool**: `build` package (PEP 517/518 compliant)
- **Setuptools**: v68+ for modern Python packaging
- **Documentation**: See README.md, CLAUDE.md, TESTING_GUIDE.md

---

**Build Status**: ✅ COMPLETE & VERIFIED
**Ready for Distribution**: ✅ YES
**Production Ready**: ✅ YES

Date: 2026-02-13
