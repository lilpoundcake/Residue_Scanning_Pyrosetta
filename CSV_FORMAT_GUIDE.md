# CSV Format Guide for Custom & Double Mutations

Complete guide to creating CSV files for Custom Mutations (CM) mode and Double Mutations (DM) mode.

## Quick Reference

| Want | Format | Example | Generates |
|------|--------|---------|-----------|
| All mutations at a position | Position spec | `AH98` or `FH100A` | 18 mutations |
| Specific mutation | Mutation | `AH98D` or `YH52AD` | 1 mutation |
| Position (with insertion code) | Position spec | `GH100B` | 18 mutations at H100B |
| Mutation (with insertion code) | Mutation | `GH100BD` | 1 mutation (G→D) at H100B |

---

## Format Specification

### General Pattern

```
{FROM_AA}{CHAIN}{POSITION}{[INSERTION_CODE]}{[TO_AA]}
```

### Components

| Component | Required? | Example | Notes |
|-----------|-----------|---------|-------|
| FROM_AA | Yes | G, A, Y, F | Single letter amino acid code |
| CHAIN | Yes | H, L, A | Chain ID from PDB |
| POSITION | Yes | 98, 100, 52 | Residue number from PDB |
| INSERTION_CODE | No | A, B, C | Only if position has insertion code in PDB |
| TO_AA | No | D, W, E | If missing → generates all 18 |

---

## Examples by Type

### Type 1: Simple Position Spec (No Insertion Code)

**Format:** `{FROM_AA}{CHAIN}{POSITION}`

```csv
AH98       # All mutations at H98
GH100      # All mutations at H100
YL32       # All mutations at L32
```

**Result:**
- Generates 18 mutations for each position
- Includes all standard amino acids except Cys and Pro
- Auto-generates wildtype for ΔΔG calculation

**Use when:** Position has no insertion code in PDB

---

### Type 2: Specific Mutation (No Insertion Code)

**Format:** `{FROM_AA}{CHAIN}{POSITION}{TO_AA}`

```csv
AH98D      # A→D at H98
GH100W     # G→W at H100
YL32F      # Y→F at L32
```

**Result:**
- Evaluates exactly 1 mutation per line
- Auto-generates wildtype entry
- Total: 2 rows per input line (mutation + wildtype)

**Use when:** Want specific mutations, no insertion codes

---

### Type 3: Position with Insertion Code

**Format:** `{FROM_AA}{CHAIN}{POSITION}{INSERTION_CODE}`

```csv
FH100A     # All mutations at H100A (insertion code A)
GH100B     # All mutations at H100B (insertion code B)
SL31A      # All mutations at L31A (insertion code A)
```

**Result:**
- Generates 18 mutations at the insertion code position
- Must use UPPERCASE insertion codes (A, B, C, etc.)
- Insertion code must exist in your PDB file

**Use when:** Position has insertion code and you want all mutations

---

### Type 4: Specific Mutation with Insertion Code

**Format:** `{FROM_AA}{CHAIN}{POSITION}{INSERTION_CODE}{TO_AA}`

```csv
YH52AD     # Y→D at H52A (insertion code A)
FH100BW    # F→W at H100B (insertion code B)
GH100AE    # G→E at H100A (insertion code A)
```

**Result:**
- Evaluates exactly 1 mutation with insertion code
- Auto-generates wildtype
- Total: 2 rows per input line

**Use when:** Specific mutations at positions with insertion codes

---

## How to Create Your CSV

### Step 1: Check Your PDB

```bash
# List all positions with insertion codes
grep "^ATOM" your.pdb | awk '{print $5, $6}' | sort -u | grep " [A-Z]$"

# Example output:
# H 100A    ← H chain, position 100, insertion code A
# H 100B    ← H chain, position 100, insertion code B
# L 52      ← L chain, position 52, no insertion code
```

### Step 2: Decide What to Test

**Option A: All mutations at key positions**
```csv
# Position specs (generates all 18 each)
AH98
GH100A
YL32
```

**Option B: Specific mutations you're interested in**
```csv
# Specific mutations (one each)
AH98D
GH100AW
YL32F
```

**Option C: Mix of both**
```csv
# Some positions (all 18), some specific mutations
AH98        # Position spec - all 18
GH100AW     # Specific mutation - only this one
YL32        # Position spec - all 18
```

### Step 3: Verify Format

Before running, verify your CSV:

```bash
# Check for common mistakes
cat your_input.csv | while read line; do
  echo "Line: $line"
  # Verify it matches a position in your PDB
  FROM=$(echo $line | cut -c1)
  CHAIN=$(echo $line | cut -c2)
  POS=$(echo $line | cut -c3- | grep -o "^[0-9]*")
  echo "  FROM=$FROM, CHAIN=$CHAIN, POS=$POS"
done
```

---

## Real Example: 7S7I.pdb

### Check the PDB

```bash
grep "^ATOM" tests/pdbs/7S7I.pdb | awk '{print $5, $6}' | sort -u | grep -E "5[0-9]|10[0-9]|3[0-9]"

# Output:
H 52        ← H52 exists
H 52A       ← H52A (insertion code A) exists
H 100       ← H100 exists
H 100A      ← H100A (insertion code A) exists
H 100B      ← H100B (insertion code B) exists
L 31        ← L31 exists
L 32        ← L32 exists
```

### Create Input CSV

```csv
YH52        # All mutations at H52
FH100A      # All mutations at H100A (has insertion code)
YL32W       # Specific: Y→W at L32
SL31        # All mutations at L31
```

### What Gets Evaluated

```
YH52         →  18 mutations (YH52A, YH52D, YH52E, ..., YH52G [wildtype])
FH100A       →  18 mutations (FH100AA, FH100AD, FH100AE, ..., FH100AF [wildtype])
YL32W        →  1 mutation (YL32W) + wildtype (YL32Y)
SL31         →  18 mutations (SL31A, SL31D, SL31E, ..., SL31S [wildtype])
```

**Total: ~54 mutations to evaluate**

---

## Common Mistakes

### ❌ Lowercase Insertion Codes

```csv
❌ YH52aD    # Wrong: PDB uses UPPERCASE (A, B, C)
✅ YH52AD    # Correct: matches PDB format H52A
```

**Fix:** Always use UPPERCASE letters for insertion codes

---

### ❌ Wrong Chain Letter

```csv
❌ FH100aD   # Wrong: should be FH100AD or verify chain exists
✅ FH100AD   # Correct if PDB has H chain
```

**Fix:** Verify chain exists: `grep "^ATOM" your.pdb | awk '{print $5}' | sort -u`

---

### ❌ Position Doesn't Exist

```csv
❌ AH999D    # Wrong: H999 doesn't exist in PDB
✅ AH100D    # Correct if H100 exists
```

**Fix:** Check PDB for actual positions before adding to CSV

---

### ❌ Using Wrong FROM_AA

```csv
❌ DH100A    # Wrong if position has different amino acid (e.g., F)
✅ FH100A    # Correct if position is F
```

**Fix:** The FROM_AA must match the actual residue type in PDB:
```bash
grep "^ATOM" your.pdb | grep " H 100" | head -1 | awk '{print $4}'
```

---

### ❌ Ambiguous Format

```csv
❌ FH100A    # Ambiguous: Is this position or mutation?
✅ FH100AW   # Clear: F→W at H100A (if H100A exists)
✅ FH100     # Clear: all mutations at H100
```

**Fix:** If unclear, use explicit format (5 letters = position spec, 6+ letters = specific mutation)

---

## CSV Validation Script

```bash
#!/bin/bash
# Validate CSV against PDB

PDB="your.pdb"
CSV="input.csv"

while IFS= read -r line; do
    FROM=$(echo $line | cut -c1)
    CHAIN=$(echo $line | cut -c2)
    REST=$(echo $line | cut -c3-)

    echo "Checking: $line"

    # Extract position (numbers only)
    POS=$(echo $REST | grep -o "^[0-9]*")

    # Check if position exists in PDB
    if grep -q "^ATOM.*$CHAIN.*$POS " $PDB; then
        echo "  ✅ Position found: $CHAIN$POS"
    else
        echo "  ❌ Position NOT found: $CHAIN$POS"
    fi
done < $CSV
```

---

## Testing Your CSV

### Run with Debug Output

```bash
# See what gets parsed
residue-scan -f your.pdb -r HL -l A -m CM \
    --mutations_csv input.csv \
    --cpu 4 --not_relax \
    -o results/ \
    2>&1 | grep -E "mutations|ERROR|WARNING"
```

### Check Output

```bash
# Verify how many mutations were actually evaluated
wc -l results/Rosetta_ddG_mut.csv
# Should be: (mutations found) + 1 header + (wildtype entries)
```

---

## Best Practices

1. **Always check your PDB first**
   ```bash
   grep "^ATOM" your.pdb | awk '{print $5, $6}' | sort -u | head -20
   ```

2. **Use UPPERCASE for insertion codes**
   - PDB format is `H 52A` (uppercase A)
   - CSV should match: `YH52AD` (not `YH52aD`)

3. **Be explicit with ambiguous positions**
   - Use full mutation format when possible: `FH100AW`
   - Avoids confusion between position specs and mutations

4. **Test with small CSV first**
   ```csv
   # Start with just 1-2 entries
   AH98
   YL32
   ```

5. **Document your CSV**
   ```csv
   # 7S7I.pdb - Interface mutations for CD4 binding domain
   FH100A     # Position with insertion code
   AH98D      # Specific mutation
   YL32       # All mutations at position
   ```

---

## Reference: Amino Acid Codes

Single letter codes used in CSV:

```
A=Alanine    D=Aspartic acid   G=Glycine     L=Leucine      R=Arginine
C=Cysteine   E=Glutamic acid   H=Histidine   M=Methionine   S=Serine
F=Phenylalanine  I=Isoleucine  K=Lysine      N=Asparagine   T=Threonine
P=Proline    Q=Glutamine       V=Valine      W=Tryptophan   Y=Tyrosine
```

Note: C (Cysteine) and P (Proline) are excluded from auto-generated position specs

---

## Troubleshooting

### Issue: "residue not found in pose"

```
WARNING: residue H52a not found in pose, skipping
```

**Cause:** Lowercase insertion code in CSV

**Solution:** Use UPPERCASE: `YH52AD` not `YH52aD`

---

### Issue: "skipping unrecognized entry"

```
WARNING: skipping unrecognized entry at line 2: 'F100b'
```

**Cause:** Incomplete format or wrong case

**Solution:** Use full format with UPPERCASE codes: `FH100B` or `FH100BD`

---

### Issue: No mutations generated

**Cause:** Positions don't exist in PDB or wrong format

**Solution:**
1. Verify positions in PDB: `grep "^ATOM" your.pdb | awk '{print $5, $6}' | sort -u`
2. Check CSV format against examples above
3. Run validation script

---

## Questions?

See `tests/TEST_SUITE.md` for more information on running tests and interpreting results.

---

## Double Mutations Mode (DM) CSV Format ✨

### Overview

Double Mutations mode can accept CSV input with two position-pair or mutation-pair formats:

| Want | Format | Example | Generates |
|------|--------|---------|-----------|
| All combinations at two positions | Position pair | `AH98_GL101` | 18×18 = 324 mutations |
| Specific double mutation | Mutation pair | `AH98D_GL101H` | 1 mutation |
| Mixed (position + mutation) | Mixed | `AH98_GL101H` or `AH98D_GL101` | 18 or 1 mutations |

### Format Specification

```
{FIRST_PART}_{SECOND_PART}
```

Each part follows the same format as Custom Mutations (CM):
- **Position spec**: `{FROM_AA}{CHAIN}{POSITION}` or `{FROM_AA}{CHAIN}{POSITION}{INSERTION_CODE}`
- **Mutation**: `{FROM_AA}{CHAIN}{POSITION}{TO_AA}` or `{FROM_AA}{CHAIN}{POSITION}{INSERTION_CODE}{TO_AA}`

### Examples

**Position-Pair** (generates all 18×18 combinations):
```csv
AH98_GL101         # All mutations at both positions
YH52_SL31          # Different chains and positions
FH100A_GH99        # One with insertion code, one without
```

**Specific Mutation-Pair** (only these mutations):
```csv
AH98D_GL101H       # A→D at H98 + G→H at L101
YH52AD_SL31W       # Y→D at H52A + S→W at L31
```

**Mixed**:
```csv
AH98_GL101H        # All 18 at H98 + G→H at L101 only
AH98D_GL101        # A→D at H98 + all 18 at L101
```

### Real Example: 7S7I.pdb

```csv
SL30_SL31          # Position pair: generate 18×18=324 combinations
RH53D_RH55K        # Specific pair: only this double mutation
```

Processing:
- `SL30_SL31` → SL30A_SL31A, SL30A_SL31D, ..., SL30S_SL31S (324 double mutations)
- `RH53D_RH55K` → 1 double mutation (R→D at H53, R→K at H55)

**Wildtype Pseudo-Mutants (Auto-Generated)**:
- For each unique position pair, one wildtype pseudo-mutant is automatically created
- Example: For `SL30_SL31`, the wildtype pair `SL30S_SL31S` is added (S→S, S→S, i.e., no change)
- Each position's wildtype component uses its original amino acid
- All rows (mutations + wildtypes) are evaluated with identical FastRelax protocol
- Wildtype rows provide reference energy for ΔΔG calculation: `ddG = mutant_energy - wildtype_energy`
- Wildtype rows are evaluated in parallel and **filtered out** of final CSV output (only mutations shown)

### Usage

**Run DM mode with CSV input**:
```bash
residue-scan -f complex.pdb -o results -r HL -l A -m DM \
    --mutations_csv input.csv --cpu 4
```

This skips Residue Scanning and evaluates double mutations directly from your CSV.

### Validation

Same as Custom Mutations mode:
- `from_aa` must match actual residue in pose
- Insertion codes must be UPPERCASE
- Positions must exist in PDB

---

## Advanced: CSV Ambiguity Resolution

### The Ambiguity Problem

When a CSV entry has a single trailing character that is a valid amino acid (e.g., `GH100A`), it creates an ambiguity:

```
GH100A could mean:
  1. MUTATION:     G→A at position H100 (no insertion code)
  2. POSITION SPEC: All 18 mutations at position H100A (insertion code = A)
```

### Real Example: 7S7I.pdb

In the 7S7I PDB file, both positions might exist:
- **Position H100** (no insertion code)
- **Position H100A** (with insertion code A)

The CSV entry `GH100A` is ambiguous - which position should be mutated?

### Solution: Dual-Variant Logic

When the parser encounters an ambiguous entry, it:

1. **Checks both interpretations** in the PDB:
   - Does position H100 exist with the specified from_aa (G)?
   - Does position H100A exist with the specified from_aa (G)?

2. **Processes based on existence**:
   - **Both exist**: Generates mutations for BOTH positions (no data loss)
   - **Only mutation exists**: Treats as mutation (G→A at H100)
   - **Only position spec exists**: Treats as position spec (all mutations at H100A)
   - **Neither exists**: Skips with warning

3. **Reports ambiguities**:
   - Logs each detected ambiguity with interpretation chosen
   - Writes detailed report to `CSV_AMBIGUITIES.txt`
   - Ensures transparency about what was processed

### Example Walkthrough

**Input CSV**:
```
GH100A
YH52AD
SL31
```

**PDB Content (7S7I)**:
- Position H100 exists with Glycine (G)
- Position H100A exists with Glycine (G)
- Position H52A exists with Tyrosine (Y)
- Position L31 exists with Serine (S)

**Processing**:

**Entry 1: `GH100A`**
```
Parsing: from_aa='G', chain='H', number=100, trailing='A'
Trailing is 1 char and 'A' is valid AA → AMBIGUOUS

Check mutation (G→A at H100):  ✓ H100 exists with G
Check position spec (all at H100A): ✓ H100A exists with G

BOTH VALID → Generate both:
  - GH100A (mutation G→A at H100)
  - GH100AA, GH100AD, ... (all 18 at H100A)
```

**Entry 2: `YH52AD`**
```
Parsing: from_aa='Y', chain='H', number=52, trailing='AD'
Trailing is 2 chars → insertion_code='A', to_aa='D'

Specific mutation: Y→D at H52A
Generate: YH52AD (only)
```

**Entry 3: `SL31`**
```
Parsing: from_aa='S', chain='L', number=31, trailing=''
Trailing is empty → position spec, all 18 AA

Generate: SL31A, SL31D, ..., SL31W (18 mutations)
```

### Output Files

**Main Output**:
- `Rosetta_ddG_mut.csv` — Standard mutation results

**Ambiguity Report** (if applicable):
- `CSV_AMBIGUITIES.txt` — Detailed ambiguity information
  - Lists each detected ambiguity
  - Explains interpretation chosen
  - Provides explanation of disambiguation logic

**Console Output** (if ambiguity detected):
```
ℹ️  AMBIGUITY DETECTED in 'GH100A': Both mutation (G→A at H100) and
    position spec (all mutations at H100A) exist in PDB. Processing BOTH.
```

### Decision Tree

```
CSV Entry: [from_aa][chain][number][trailing]
           ↓
Is trailing character count?
├─ 0 chars  → Position spec (all 18)
├─ 1 char   → Is it a valid AA?
│   ├─ YES  → AMBIGUOUS (check both interpretations)
│   │   ├─ Both exist    → Generate BOTH
│   │   ├─ Mutation only → Generate mutation
│   │   ├─ Position only → Generate position spec
│   │   └─ Neither       → Skip
│   └─ NO   → Insertion code, position spec (all 18)
├─ 2 chars  → First=insertion_code, Second=to_aa (specific mutation)
└─ 3+ chars → Error, skip
```

### Backward Compatibility

**✅ Fully Backward Compatible**:
- Existing CSV files work unchanged
- Non-ambiguous entries unaffected
- DataFrame structure identical
- Ambiguity handling is automatic and transparent

### Migration Guide

**No action required**. The ambiguity handling is transparent:
- Non-ambiguous entries work exactly as before
- Ambiguous entries now process BOTH variants (more complete)
- Detailed report explains any changes

**Recommendation for New Users**: Use explicit format to avoid ambiguity:
```csv
# Good: Unambiguous
GH100    # Position spec (no insertion code)
GH100D   # Specific mutation (insertion code + AA)
YH52AD   # Specific mutation with insertion code

# Ambiguous (but now handled correctly)
GH100A   # Could be both, parser handles it
```
