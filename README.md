# Residue Scanning Pyrosetta

Script for checking mutations in a protein-protein interface

## Installation

Clone the repo and create Conda/Micromamba environment:

```bash
git clone https://github.com/lilpoundcake/Residue_Scanning_Pyrosetta.git

micromamba env create -f env.yaml
```

## Usage

Activate the environment:

```bash
micromamba activate rosetta_RS
```
And use the script:

```bash
python rosetta_RS.py \
-f {PATH_TO_PDB_FILE} \
-o {PATH_TO_OUTPUT_FOLDER} \
-r {NAMES_OF_RECEPTOR_CHAINS} \
-l {NAMES_OF_LIGAND_CHAINS} \
-mode {MODE} \
--cpu {CPU_COUNT}
```

### Required keys

 - `-f`, `--file` - Path to PDB file containing the protein complex
 - `-m`, `--mode` - Analysis protocol to execute (short or long versions accepted)
    - `FastRelax`,`FR` - protocol for relaxation of the input structure only (without Residue Scanning)
    - `Residue_Scanning`, `RS` - protocol for ddG calculation for every amino acid in the interface
    - `Double_Mut_Searching`, `DM` - for searching double mutations in the interface. Performs `Residue_Scanning`, filters single mutations using the [pandas-like](https://pandas.pydata.org/docs/reference/api/pandas.DataFrame.query.html) `condition` (Default is `ddG_complex < 0.5 and ddG_interface < 0.5`), then combines nearby mutations within the specified `radius`
 - `-r`, `--receptor` and `-l`, `--ligand` - chain names of receptor and ligand
 - `--cpu` - number of CPUs for parallelisation

### Optional keys

 - `--output` - directory for output files (Default: folder "results" will be created)
 - `--not_relax` - skip the FastRelax minimization step (~10 minutes per structure). **Note**: This significantly reduces prediction accuracy
 - `--condition` - filtering criteria for single mutations in `Double_Mut_Searching` mode (Default: `ddG_complex < 0.5 and ddG_interface < 0.5`)
 - `--replics` - maximum number of conformational replicas for side-chain repacking
 - `--radius` - distance cutoff (Å) for considering amino acid pairs as neighbors in `Double_Mut_Searching` mode

## Examples

Common Residue Scanning task:

```bash
python rosetta_RS.py \
-f complex.pdb -o RS_results \
-r HL -l A \
-mode RS --cpu 110
```

Double mutation search with custom conditions:

```bash
python rosetta_RS.py \
-f complex.pdb -o results \
-r HL -l A \
-mode DM --cpu 110 \
--radius 10 --condition "ddG_complex < 0 and ddG_interface < 0"
```

Quick analysis without minimization:

```bash
python rosetta_RS.py \
-f complex.pdb -o results \
-r HL -l A \
-mode DM --cpu 10 \
--not_relax 
```

## Output files and folders

The following files and folders will be generated in your specified output directory:

 - `FastRelax.pdb` - the energetically minimized structure with the lowest Rosetta Energy Units (REU) from the FastRelax protocol
 - `Rosetta_ddG_mut.csv` - CSV table containing ΔΔG values for all single-point mutations analyzed:
    - `ddG_complex`: free energy change for the entire complex. Indicates stability
    - `ddG_interface`: binding free energy change
 - `Rosetta_ddG_mut_double.csv` (Only in `Double_Mut_Searching` mode) - CSV table containing ΔΔG values for all double mutations analyzed, with the same columns as the single mutation table

 - `min_REU_pdb` - lowest-energy structures for each single mutation
 - `min_REU_pdb_double` - lowest-energy structures for each double mutation

## Citation

 Protocol inspired by an article **"Prediction of Protein Mutational Free Energy: Benchmark and Sampling Improvements Increase Classification Accuracy"** <br></br>
 DOI: [10.3389/fbioe.2020.558247](https://www.frontiersin.org/journals/bioengineering-and-biotechnology/articles/10.3389/fbioe.2020.558247/full)
