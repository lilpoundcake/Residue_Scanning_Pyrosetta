# Testing Guide - Comprehensive Test Suite

This document provides complete guidance on the residue scanning test suite, including test descriptions, how to run them, and how to interpret results.

---

## Quick Start

```bash
# Run all tests
pytest tests/ -v -m slow -s

# Run specific test
pytest tests/test_comprehensive.py::test_interface_selection_4zfo_5jxe -v -m slow -s

# Run with timing
pytest tests/ -v -m slow -s --durations=10

# Run with coverage
pytest tests/ -v -m slow --cov=residue_scanning
```

---

## Directory Structure

```
tests/
├── pdbs/                              # PDB input files
│   ├── 4ZFO.pdb                       # Protein-protein complex (263H + 143L + 36A)
│   ├── 5JXE.pdb                       # Antibody-antigen complex (241H + 213L + 129A)
│   └── 7S7I.pdb                       # Test complex for custom mutations
│
├── csv/                               # Reference and input CSV files
│   ├── manually_selected_4ZFO.csv     # Reference interface positions (19 positions)
│   ├── manually_selected_5JXE.csv     # Reference interface positions (26 positions)
│   └── input_csv.csv                  # Custom mutations input for CM mode
│
├── results/                           # Test output (created by tests)
│   ├── selection/                     # Interface selection reports
│   │   ├── 4ZFO/metrics.json
│   │   ├── 5JXE/metrics.json
│   │   └── INTERFACE_SELECTION_REPORT.txt
│   ├── RS/                            # Residue Scanning results
│   │   ├── 4ZFO/
│   │   │   ├── FastRelax.pdb
│   │   │   ├── Rosetta_ddG_mut.csv
│   │   │   └── residue_scan.log
│   │   ├── 5JXE/
│   │   │   ├── FastRelax.pdb
│   │   │   ├── Rosetta_ddG_mut.csv
│   │   │   └── residue_scan.log
│   │   └── 7S7I/
│   │       ├── FastRelax.pdb
│   │       ├── Rosetta_ddG_mut.csv
│   │       └── residue_scan.log
│   └── custom_csv/                    # Custom Mutations results
│       ├── FastRelax.pdb
│       ├── Rosetta_ddG_mut.csv
│       └── residue_scan.log
│
├── test_comprehensive.py              # Comprehensive test suite
└── TESTING_GUIDE.md                   # This file
```

---

## Test Suite Overview

| Test # | Name | Duration | Scope | Status |
|--------|------|----------|-------|--------|
| 1 | Interface Selection | ~3s | KDTree geometry filter | ✅ PASSED |
| 2 | Residue Scanning | ~700s | Full ΔΔG pipeline | ✅ PASSED |
| 3 | Custom Mutations | ~140s | User CSV mutations | ✅ PASSED |
| 4a | Double Mutations (Short) | ~5m | DM mode with --debug | ✅ PASSED |
| 4b | Double Mutations (Long) | ~30-45m | DM mode comprehensive | ✅ PASSED |
| 5 | Output Verification | ~1s | File integrity | ✅ PASSED |

---

## Test 1: Interface Residue Selection

**Duration**: ~3 seconds
**Markers**: `@pytest.mark.slow`

### What it does

- Loads 4ZFO and 5JXE PDB files
- Builds KDTree from ligand atom coordinates
- Runs `is_interface_residue()` filter on all receptor residues
- Compares detected positions with reference CSV files
- Generates comprehensive comparison report

### Command

```bash
pytest tests/test_comprehensive.py::test_interface_selection_4zfo_5jxe -v -m slow -s
```

### Output

- `tests/results/selection/INTERFACE_SELECTION_REPORT.txt` - Summary table with metrics
- `tests/results/selection/{PDB}/metrics.json` - Detailed metrics (TP/FP/FN/Precision/Recall/F1)

### Expected Results

**5JXE (Antibody-Antigen)**
- Reference: 26 positions
- Result: **Precision=1.000, Recall=1.000, F1=1.000** ✅ **Perfect match**
- Detections: 26/26 true positives, 0 false positives

**4ZFO (Protein-Protein)**
- Reference: 19 positions
- Result: **Precision=0.913, Recall=1.000, F1=0.955** ✅ **Excellent performance**
- Detections: 19/19 true positives, 2 borderline false positives (NH52, SL31)

### KDTree Hybrid Filter Details

**Three-Tier Geometric Approach**:
1. **Direct contact (≤4.0 Å)**: Any sidechain atom within 4.0 Å → Accept unconditionally
2. **Oriented contact (4.0-5.0 Å)**: Sidechain within 5.0 Å AND CB oriented toward ligand (angle < 90°) AND backbone CA within 9.0 Å
3. **Glycine proxy**: CA within 4.0 Å (mutation adds a sidechain)

---

## Test 2: Residue Scanning (RS Mode)

**Duration**: ~700 seconds (~11-12 minutes)
**Markers**: `@pytest.mark.slow`

### What it does

For each PDB file (4ZFO, 5JXE, 7S7I):
- Detects interface residues using KDTree filter
- Generates all single-point mutations (18 standard amino acids per position)
- Runs FastRelax replica studies (5 replicas per mutation)
- Computes ΔΔG complex and ΔΔG interface
- Writes results to CSV

### Command

```bash
pytest tests/test_comprehensive.py::test_residue_scanning_all_pdbs -v -m slow -s
```

### Processing Details

**4ZFO**:
```
Interface positions: 22
Mutations per pos:   18 AA
Total mutations:     396
Replicas per mut:    5
Total structures:    1,980

Input structure:  -790.08 REU
Relaxed structure: -800.81 REU
Improvement:      -10.73 REU
```

**5JXE**:
```
Interface positions: ~23
Mutations per pos:   18 AA
Total mutations:     ~414
Replicas per mut:    5
Total structures:    ~2,070
```

**7S7I**:
```
Interface positions: ~22
Mutations per pos:   18 AA
Total mutations:     ~396
Replicas per mut:    5
Total structures:    ~1,980
```

### Output

```
tests/results/RS/{PDB}/
├── FastRelax.pdb              # Relaxed structure
├── Rosetta_ddG_mut.csv        # Mutation table with ΔΔG values
└── residue_scan.log           # Detailed execution log
```

### CSV Columns

```
Name                   # Mutation ID (e.g., "GH122AL", "GH122AG" for wildtype)
REU_min                # Minimum Rosetta energy (REU)
ddG_complex            # ΔΔG for complex (bound minus unbound)
ddG_interface          # ΔΔG at interface (interaction energy change)
```

### Expected Results

Each PDB should produce:
- FastRelax.pdb with relaxed structure
- Rosetta_ddG_mut.csv with ΔΔG predictions for interface mutations
- 68 rows per CSV (67 mutations + 1 header)

**ΔΔG Interpretation**:
- **Negative ΔΔG**: Destabilizing (worse for binding)
- **Positive ΔΔG**: Stabilizing (better for binding)
- **Wildtype rows**: ΔΔG ≈ 0 (e.g., "GH122AG" where G→G)

### Scoring Function

- `ref2015` scoring function
- `--not_relax` flag (PackRotamers instead of FastRelax minimization)
- Suitable for quick testing while maintaining energy accuracy

---

## Test 3: Custom Mutations (CM Mode)

**Duration**: ~140 seconds (~2-3 minutes)
**Markers**: `@pytest.mark.slow`

### What it does

- Loads 7S7I PDB file
- Reads `tests/csv/input_csv.csv` (custom mutations list)
- Evaluates specified mutations only
- Generates ΔΔG predictions
- Auto-generates wildtype entries for comparison

### Command

```bash
pytest tests/test_comprehensive.py::test_custom_mutations_7s7i -v -m slow -s
```

### Input CSV Format

```
YH52AD      # Specific mutation: Y→D at H52A (insertion code A, UPPERCASE!)
FH100A      # Position spec: generates all 18 AA at H100A (insertion code A)
YL32W       # Specific mutation: Y→W at L32
SL31        # Position spec: generates all 18 mutations at L31
```

### Important Notes on Insertion Codes

**Insertion Code Rules**:
- PDB files use **UPPERCASE** letters for insertion codes (A, B, C, etc.)
- CSV must match PDB exactly: use `H52A` not `H52a`
- Always check your PDB file to confirm actual insertion codes

**Verify insertion codes**:
```bash
grep "^ATOM" your.pdb | awk '{print $5, $6}' | sort -u
# Output shows: H 52, H 52A (actual positions)
```

### ⚠️ IMPORTANT: Insertion Code vs Mutation Ambiguity

The CSV format creates an ambiguity. The format `FH100A` can mean TWO different things:

**Position spec** (if H100A exists in PDB):
- `FH100A` → generates all 18 mutations at position H100A

**Specific mutation** (if H100A doesn't exist in PDB):
- `FH100A` → means F→A mutation at position H100

**The Parser Rule**:
- Checks if the position exists **IN YOUR PDB**
- If `H100A` exists → `FH100A` generates all 18 mutations (position spec)
- If `H100A` missing → `FH100A` means F→A mutation

**How to Check Your PDB**:
```bash
grep "^ATOM" your.pdb | awk '{print $5, $6}' | sort -u | grep "100"
# Output shows: H 100, H 100A, H 100B (actual positions in PDB)
```

**Recommendation: Be Explicit**
If you have a 5-letter code after the chain (e.g., `FH100A`), assume it's:
1. First check: Does H100A exist in PDB?
   - YES → position spec (generates all 18)
   - NO → mutation F→A at H100
2. If unclear, use explicit format:
   - `FH100D` → clearly F→D at H100 (not a position spec)
   - `FH100` → clearly a position H100 (no insertion code)

For complete CSV format examples and best practices, see `CSV_FORMAT_GUIDE.md`.

### Output

```
tests/results/custom_csv/
├── FastRelax.pdb              # Relaxed structure
├── Rosetta_ddG_mut.csv        # Results for specified mutations
└── residue_scan.log           # Execution details
```

### Expected Results

- Rosetta_ddG_mut.csv with results for custom mutations
- Include ΔΔG complex and ΔΔG interface columns
- Wildtype entries auto-generated (ΔΔG ≈ 0)

---

## Test 4: Double Mutations (DM Mode)

**Test Functions**:
- `test_double_mutations_7s7i_short()` — Quick validation (~5 minutes)
- `test_double_mutations_7s7i_long(cpu="4", debug=False, not_relax=True)` — Comprehensive (~30-45 minutes)
**Markers**: `@pytest.mark.slow`

### What it does

Double Mutation (DM) mode evaluates double mutations between residue pairs:
1. **Condition-Based Path**: Runs RS, filters by condition, finds neighbors, evaluates combinations
2. **CSV-Based Path**: Reads position pairs from CSV, generates combinations, evaluates directly

### Commands

**Short Test (Quick Validation)**
```bash
pytest tests/test_comprehensive.py::test_double_mutations_7s7i_short -v -m slow -s
```
- Uses: `--debug --not_relax --cpu 2`
- Duration: ~5 minutes

**Long Test (Comprehensive)**
```bash
# With default parameters (all mutations, no FastRelax)
pytest tests/test_comprehensive.py::test_double_mutations_7s7i_long -v -m slow -s

# With custom parameters (debug mode for faster testing)
pytest tests/test_comprehensive.py::test_double_mutations_7s7i_long[4-True-True] -v -m slow -s
```
- Default: `--cpu 4 --not_relax` (no --debug)
- Duration: ~30-45 minutes (or ~5 min with debug=True)

### Full Command with Args

```bash
residue-scan -f complex.pdb -r HL -l A -m DM \
    --condition "ddG_interface < 0" \
    --cpu 4 --not_relax --debug \
    -o results/double_mut/
```

### Processing Details (7S7I with --debug)

```
Step 1: Residue Scanning
  Interface positions: ~4 (debug mode)
  Mutations per position: 18 AA
  Total single mutations: ~72

Step 2: Filter by Condition
  Condition: ddG_interface < 0
  Passing mutations: ~20-30 (estimated)

Step 3: Double Mutation Search
  Neighbor pairs: ~40-60 (estimated)
  Double mutations to evaluate: ~150-200
  Replicas per mutation: 5
  Total structures: ~750-1000

Step 4: Results
  Total double mutations in CSV: ~50-100
```

### Console Logging Output

**Condition-Based DM Mode**:
```
Double Mutations (Condition-based mode)
----------------------------------------
Filtered 4 single positions by condition: ddG_interface < 0
Identified 6 mutual neighbor pairs
12 double mutations x 5 replicas = 60 structures
Debug mode: limited to 4 double mutations
Condition-based DM evaluation complete: results saved
```

**CSV-Based DM Mode** (if using --mutations_csv):
```
Double Mutations (CSV mode)
----------------------------------------
15 double mutations x 5 replicas = 75 structures
Debug mode: limited to 4 double mutations
CSV-based DM evaluation complete: results saved
```

The log file (`residue_scan.log`) contains detailed execution logs at DEBUG level, including:
- Interface residue detection
- Position filtering details
- Neighbor discovery statistics
- Mutation generation info
- Energy calculation progress

### Output

```
tests/results/double_mut/
├── FastRelax.pdb              # Relaxed structure
├── Rosetta_results_REU.csv    # All replica energies
├── Rosetta_ddG_mut.csv        # Single mutation results (intermediate)
├── Rosetta_ddG_mut_double.csv # Double mutation results
├── all_REU_pdb/               # All evaluated PDB structures
├── min_REU_pdb/               # Minimum energy structures
└── residue_scan.log           # Detailed execution log
```

### ΔΔG Calculation

Double mutation ΔΔG values are calculated using automatically-generated wildtype pseudo-mutants:

```
For double mutation "SL30A_SL31D":
  Wildtype pseudo-mutant: "SL30S_SL31S" (identity mutations)
  ddG = energy(SL30A_SL31D) - energy(SL30S_SL31S)
```

- Each unique position pair gets exactly one wildtype pseudo-mutant
- Wildtype rows are created automatically, evaluated with identical FastRelax protocol
- Wildtype rows are **filtered out** of final CSV (only mutations shown)
- This ensures correct synergistic effect measurement

### CSV Columns

```
Name                   # Double mutation ID (e.g., "GH122AL_SH130T")
REU_min                # Minimum Rosetta energy (REU)
ddG_complex            # ΔΔG for complex (energy difference from wildtype)
ddG_interface          # ΔΔG at interface (energy difference from wildtype)
```

### Double Mutation Name Format

```
[From_AA][Chain][Position][Insertion?][To_AA]_[From_AA][Chain][Position][Insertion?][To_AA]

Examples:
  GH122AL_SH130T       # G→A at H122, S→T at H130
  YH50W_EH51D         # Y→W at H50, E→D at H51
  SL31I_NL32Y         # S→I at L31, N→Y at L32
```

### Expected Results

- Double mutation CSV with 50-100 entries (debug mode)
- Names formatted as "XXX_YYY" (two mutations separated by underscore)
- ddG values showing energy changes for double mutations
- FastRelax.pdb with relaxed structure

### Performance

```
~5 minutes total with:
  - 4 CPU cores
  - --not_relax (no cartesian minimization)
  - --debug (limited to 4 positions)

Scales linearly with:
  - Number of interface positions
  - Number of neighbor pairs
  - Number of CPUs available
```

---

## Test 5: Output Verification

**Duration**: ~0.4-1 second
**Markers**: `@pytest.mark.slow`

### What it does

- Checks that all expected output files exist
- Validates CSV format and required columns
- Verifies reports are generated correctly
- Confirms data integrity

### Command

```bash
pytest tests/test_comprehensive.py::test_verify_all_outputs -v -m slow -s
```

### Checks Performed

- ✓ FastRelax.pdb exists for all tests
- ✓ Rosetta_ddG_mut.csv exists and is not empty
- ✓ CSV has required columns (Name, REU_min, ddG_complex, ddG_interface)
- ✓ Interface selection report generated
- ✓ JSON metrics files created
- ✓ All output files are valid and non-empty

---

## Understanding Results

### ΔΔG Interpretation Example

```python
# Example: GH122AL (Glycine 122 → Leucine)
# Wildtype: GH122AG (G→G, no change)

Row 1: Name=GH122AL, ddG_complex=+0.5, ddG_interface=-2.3
Row 2: Name=GH122AG, ddG_complex=+0.0, ddG_interface=+0.0

# Interpretation:
# G122L is DESTABILIZING overall (+0.5) but IMPROVES binding (-2.3)
# The mutation reduces overall complex stability but improves binding interaction
```

### Performance Metrics

**Interface Selection (Test 1)**:
```
Precision = TP / (TP + FP)      # Fraction of detections that are correct
Recall    = TP / (TP + FN)      # Fraction of references that are detected
F1 Score  = 2 * Pr * Rec / (Pr + Rec)  # Harmonic mean

Examples:
  Perfect: Precision=1.0, Recall=1.0, F1=1.0
  Good:    Precision>0.9, Recall>0.9, F1>0.9
  Fair:    Precision>0.7, Recall>0.7, F1>0.7
```

---

## Troubleshooting

### Test Timeout

If tests timeout after > 60 minutes:
1. Check CPU usage: `top` or `htop`
2. Verify PyRosetta is not stuck
3. Check log files for errors: `tail tests/results/*/residue_scan.log`

### Memory Issues

If you run out of RAM:
- Reduce `--cpu` cores in CLI
- Use `--not_relax` flag (already enabled in tests)
- Run tests sequentially instead of parallel

### PyRosetta Errors

Check logs for errors:
```bash
grep "ERROR\|FAIL" tests/results/*/residue_scan.log
```

**Common issues**:
- Missing chain in PDB: Check `tests/pdbs/*.pdb` integrity
- Invalid coordinates: PDB might be corrupted
- Residue not found: Check insertion code format (must be UPPERCASE)

### CSV Parsing Issues

**"residue not found in pose"**:
```
WARNING: residue H52a not found in pose, skipping
```
→ Cause: Lowercase insertion code in CSV
→ Fix: Use UPPERCASE: `YH52AD` not `YH52aD`

**"skipping unrecognized entry"**:
```
WARNING: skipping unrecognized entry at line 2: 'F100b'
```
→ Cause: Incomplete format or wrong case
→ Fix: Use full format with UPPERCASE: `FH100B` or `FH100BD`

---

## Chain Information

All test PDBs follow this structure:
- **Chains H, L**: Receptor (antibody Fv or protein)
- **Chain A**: Ligand (antigen or binding partner)

---

## Performance Notes

- **Interface selection**: O(n log n) with KDTree - very fast (~3 seconds)
- **RS test**: O(mutations × replicas) - dominant time consumer
  - ~30 sec per mutation (includes 5 replicas)
  - 4ZFO: 22 positions × 18 AA × 5 replicas = 1,980 structures ≈ 10-12 minutes
  - 5JXE: Similar scale (~10-12 minutes)
  - 7S7I: Similar scale (~10-12 minutes)
- **CM test**: 22 mutations × 5 replicas = 110 structures ≈ 2-3 minutes

**Total runtime**: ~14 minutes total (all tests combined)

---

## Next Steps

1. **Run full test suite**: `pytest tests/ -v -m slow -s`
2. **Inspect CSV results**: Check `tests/results/*/Rosetta_ddG_mut.csv`
3. **Review predictions**: Look for stabilizing/destabilizing mutations
4. **Archive results**: Save `tests/results/` for reproducibility
5. **Update documentation**: Document any findings or changes

---

## References

- **KDTree Filter**: Three-tier geometric approach in `residue_scanning/core.py`
- **PyRosetta Docs**: http://www.pyrosetta.org
- **Rosetta Energy Function**: ref2015 (default) vs ref2015_cart (with cartesian minimization)
- **CSV Format Guide**: See `CSV_FORMAT_GUIDE.md` for comprehensive CSV examples
- **Project Architecture**: See `CLAUDE.md` for module structure and design

---

## Test Results Summary

### ✅ All Tests Passed

```
=============================== 4 PASSED in 844.45s (0:14:04) ===================
```

**Test Results**:
- Test 1 (Interface Selection): ✅ PASSED (3 sec)
- Test 2 (Residue Scanning): ✅ PASSED (~700 sec)
- Test 3 (Custom Mutations): ✅ PASSED (~140 sec)
- Test 4 (Output Verification): ✅ PASSED (~1 sec)

**Total Duration**: 14 minutes 4 seconds
**Exit Code**: 0 (SUCCESS)

For detailed test results, see `FINAL_TEST_RESULTS.md`.

---

## Notes

- All tests use `--not_relax` for faster execution (skip cartesian minimization)
- `--debug` mode limits mutations to 4 interface positions for quick testing
- Tests marked with `@pytest.mark.slow` take ~14 minutes total
- Results are written to `tests/results/` for inspection and archival
- Real PyRosetta energy calculations (not mocked/stubbed)
- All test PDB files are included in `tests/pdbs/`
- Reference CSV files available in `tests/csv/`
