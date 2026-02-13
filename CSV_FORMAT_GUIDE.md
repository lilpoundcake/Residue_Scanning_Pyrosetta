# CSV Format Guide for Custom Mutations

Complete guide to creating `input_csv.csv` for Custom Mutations (CM) mode.

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

See `TESTING_GUIDE.md` for more information on running tests and interpreting results.
