# Project Memory

## Environment
- PyRosetta is in the `pyrosetta` conda env (Python 3.9), NOT `rosetta_RS`
- Use `/Users/lilpoundcake/micromamba/envs/pyrosetta/bin/python` and `pip`
- `blosc` must be installed for efficient PyRosetta pickle (pose serialization)
  - Without it, multiprocessing hangs on macOS (spawn method, poses too large to pickle)
  - Fix: `pip install blosc` in the pyrosetta env

## Multiprocessing on macOS
- macOS uses `spawn` start method by default since Python 3.8
- Workers start fresh — `pyrosetta.init()` is NOT called automatically
- Solution: `mp.Pool(n, initializer=_worker_init)` where `_worker_init` calls `pyrosetta.init(PYROSETTA_FLAGS)`
- `_worker_init` is defined in `core.py` and imported by `cli.py`

## Package structure
- Package: `residue_scanning/` (cli.py, core.py, preprocessing.py, __init__.py)
- Entry point: `residue-scan` → `residue_scanning.cli:main`
- pyproject.toml: `requires-python = ">=3.9"`, `numpy>=1.20`
- Tests: `tests/test_debug.py` uses `@pytest.mark.slow`, calls `main()` directly

## Running tests
```bash
/Users/lilpoundcake/micromamba/envs/pyrosetta/bin/pytest tests/ -v -m slow -s
```
Test takes ~12 min (1BRS, debug mode, --not_relax, 2 CPUs).
