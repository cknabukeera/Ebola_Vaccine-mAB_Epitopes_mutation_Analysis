# Antibody-GP Complex Binding Energy Analysis
## Inmazeb individual component analysis

Since this is a complex of three antibodies, I need to extract the individual GP proteins and its corresponding antibody, to have the GP interacting individually with each antibody

## Background & Rationale

### The Problem

**Bundibugyo ebolavirus (BDBV)** is one of the six known species of the genus *Ebolavirus* and has been responsible for multiple outbreaks, including the 2007-2008 outbreak in Uganda and the 2012 outbreak in the Democratic Republic of Congo. Despite its pathogenicity, there are currently **no licensed therapeutics specifically designed against BDBV**.

### The Antibody Cocktail: Inmazeb

**Inmazeb** (REGN-EB3) is a triple monoclonal antibody cocktail developed by Regeneron Pharmaceuticals and approved by the FDA in October 2020 for the treatment of **Zaire ebolavirus (ZEBOV)** infection. Inmazeb consists of three antibodies:
- **REGN3470**
- **REGN3471** 
- **REGN3479**

These antibodies were developed against the **ZEBOV GP** (glycoprotein) and have demonstrated potent neutralization activity against Zaire ebolavirus.

### The Knowledge Gap

While Inmazeb is highly effective against ZEBOV, the **efficacy of these antibodies against BDBV remains poorly characterized**. Key questions include:

1. **Do REGN3470, REGN3471, and REGN3479 bind effectively to BDBV GP?**
2. **What is the impact of BDBV-specific mutations on antibody-GP binding affinity?**
3. **Can Inmazeb or its individual components serve as therapeutics against BDBV?**

### Our Approach

To address these questions, we have taken a structure-based approach:

1. **Crystal Structure**: We utilized the cryo-EM structure of the Inmazeb-GP complex (PDB: 7TN9), which contains ZEBOV GP bound to the three antibodies.

2. **Extraction of Individual Complexes**: We extracted three separate complexes from 7TN9:
   - **GP-REGN3470**: ZEBOV GP with REGN3470 (heavy: G,I; light: H,J)
   - **GP-REGN3471**: ZEBOV GP with REGN3471 (heavy: A,C,E; light: B,D,F)
   - **GP-REGN3479**: ZEBOV GP with REGN3479 (heavy: M,O,Q; light: N,P,R)

3. **Introduction of BDBV Mutations**: We introduced mutations that differ between ZEBOV and BDBV GP into the ZEBOV GP structure. These mutations are located at key positions in the GP that may affect antibody binding.

4. **Binding Energy Calculations**: Using PyRosetta's Flex ddG protocol, we calculated the change in binding free energy (ΔΔG) upon mutation to predict how these substitutions affect antibody-GP interactions.

### BDBV Mutations Analyzed

#### GP-REGN3470 Complex (Chain T)
- **S263N**: Serine to Asparagine at position 263
- **K265R**: Lysine to Arginine at position 265
- **E280T**: Glutamic Acid to Threonine at position 280

#### GP-REGN3471 Complex (Chain T)
- **E112D**: Glutamic Acid to Aspartic Acid at position 112
- **P116A**: Proline to Alanine at position 116

#### GP-REGN3479 Complex (Chain V)
- **I504T**: Isoleucine to Threonine at position 504
- **V505L**: Valine to Leucine at position 505
- **N507S**: Asparagine to Serine at position 507
- **N507R**: Asparagine to Arginine at position 507

### Expected Outcomes

- **ΔΔG > 0**: Mutation destabilizes antibody-GP binding → BDBV may escape neutralization
- **ΔΔG < 0**: Mutation stabilizes antibody-GP binding → BDBV may be more susceptible
- **ΔΔG ≈ 0**: Mutation has minimal effect on binding

These results will inform whether Inmazeb or its individual components could be repurposed as therapeutics against BDBV, or if antibody engineering is required to enhance cross-reactivity.

---

## Requirements

- PyRosetta 2026+
- Python 3.11+
- SLURM HPC cluster (optional, for high-throughput calculations)

---

### Cleaning the 7TN9 PDB file 
 - First, we need to clean the PDB file to remove the water molecules, ligands and non protein atoms
```
## PDB Cleaning Script

### `clean_pdb.py`

The following script removes water, ligands, and non-protein atoms from any PDB file. It uses a pure text-filtering approach that works on any PDB without rebuilding bugs.

```python
#!/usr/bin/env python3

"""
clean_pdb.py - Remove water, ligands, and non-protein atoms from any PDB.
Optionally keep only specific chains.

Pure text-filtering approach -> works on ANY PDB, no rebuilding bugs.
"""
import argparse
import sys

# Standard 20 amino acids (3-letter codes)
STANDARD_AAS = {
    "ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"
}

def main():
    parser = argparse.ArgumentParser(description="Clean any PDB for Rosetta")
    parser.add_argument("-i", "--input", required=True, help="Input PDB file")
    parser.add_argument("-o", "--output", required=True, help="Output cleaned PDB")
    parser.add_argument("-c", "--keep_chains", default=None,
                        help="Comma-separated chains to keep (e.g., A,B). If not set, keeps all.")
    parser.add_argument("--remove_water", action="store_true", default=True,
                        help="Remove water (HOH)")
    parser.add_argument("--remove_ligands", action="store_true", default=True,
                        help="Remove non-canonical ligands/ions")
    parser.add_argument("--keep_hetatm", action="store_true", default=False,
                        help="Keep HETATM lines (by default, only ATOM lines are kept)")
    
    args = parser.parse_args()

    # Parse chain filter
    keep_chains = set()
    if args.keep_chains:
        keep_chains = set(args.keep_chains.upper().split(','))
        print(f"Keeping only chains: {', '.join(keep_chains)}")

    filtered_lines = []
    total_lines = 0
    kept_lines = 0

    with open(args.input, 'r') as f:
        for line in f:
            total_lines += 1
            
            # If not keeping HETATM, only accept ATOM lines
            if not args.keep_hetatm and not line.startswith("ATOM"):
                continue
            
            # Keep TER lines (they mark chain ends)
            if line.startswith("TER"):
                filtered_lines.append(line)
                kept_lines += 1
                continue

            # Only process ATOM or HETATM lines
            if not line.startswith(("ATOM", "HETATM")):
                continue

            # Extract chain (PDB column 22, 0-based index 21)
            if len(line) < 22:
                continue  # malformed line, skip
            chain = line[21]
            
            # Extract residue name (columns 18-20, 0-based index 17-20)
            resname = line[17:20].strip().upper()

            # --- Chain filter ---
            if keep_chains and chain not in keep_chains:
                continue

            # --- Water filter ---
            if args.remove_water and resname == "HOH":
                continue

            # --- Ligand filter ---
            if args.remove_ligands and resname not in STANDARD_AAS:
                # Check if it's a modified residue we might want to keep (e.g., MSE -> MET)
                # By default, we skip it to be safe.
                continue

            # If we made it here, keep this atom
            filtered_lines.append(line)
            kept_lines += 1

    # Write the cleaned PDB
    with open(args.output, 'w') as f:
        f.writelines(filtered_lines)

    print(f"Input: {args.input} ({total_lines} lines)")
    print(f"Kept: {kept_lines} lines")
    print(f"Cleaned PDB saved to: {args.output}")

if __name__ == "__main__":
    main()
```

### Usage

```bash
# Clean full PDB (remove water and ligands)
python3 clean_pdb.py -i 7TN9.pdb -o 7TN9_clean.pdb

# Keep only specific chains
python3 clean_pdb.py -i 7TN9.pdb -o 7TN9_GP_REGN3470.pdb -c G,H,I,J,S,T,U,V,W,X
```

### Options

| Option | Description |
|--------|-------------|
| `-i, --input` | Input PDB file (required) |
| `-o, --output` | Output cleaned PDB file (required) |
| `-c, --keep_chains` | Comma-separated chains to keep (e.g., A,B,C) |
| `--remove_water` | Remove water molecules (HOH) (default: True) |
| `--remove_ligands` | Remove non-canonical ligands (default: True) |
| `--keep_hetatm` | Keep HETATM lines (default: False) |
```

```
### 1. Using PyMOL to extract the Individual antibody chains and combine them with the Zaire GP
#### GP_REGN3479 chains 

```
# Load the cleaned PDB
load 7TN9_clean.pdb

# Select GP1 (S,T,U) and GP2 (V,W,X)
select GP1, chain S+T+U
select GP2, chain V+W+X
select GP_full, GP1 or GP2

# Select REGN3479 antibody (heavy: M,O,Q; light: N,P,R)
select REGN3479, chain M+O+Q+N+P+R

# Combine them
select GP_REGN3479_full, GP_full or REGN3479

# Visualise
hide everything
show cartoon, GP_REGN3479_full
color blue, chain S+T+U+V+W+X
color red, chain M+O+Q
color yellow, chain N+P+R

# Extract and save
extract complex_obj, GP_REGN3479_full
save 7TN9_GP_REGN3479_with_GP2.pdb, complex_obj
```

### GP_REGN3470 and GP_REGN3471
```# Load your cleaned PDB
load 7TN9_clean.pdb

# --- GP full (GP1 + GP2) ---
select GP1, chain S+T+U
select GP2, chain V+W+X
select GP, GP1 or GP2

# --- REGN3471 (heavy: A,C,E; light: B,D,F) ---
select REGN3471, chain A+B+C+D+E+F
select GP_REGN3471, GP or REGN3471 # combine both
extract obj_REGN3471, GP_REGN3471 
save 7TN9_GP_REGN3471.pdb, obj_REGN3471

# --- REGN3470 (heavy: G,I; light: H,J) ---
select REGN3470, chain G+H+I+J
select GP_REGN3470, GP or REGN3470
extract obj_REGN3470, GP_REGN3470
save 7TN9_GP_REGN3470.pdb, obj_REGN3470
```

### 2. In silico mutagenesis(introducing BDBV mutations in Zaire) and ddG energy claculations with pyrosetta
Here, we introduce BDBV epitope mutations into Zaire GP antibody complexes to access the impact on the binding energy (The cost on the energy, if the antibody is to bind the BDBV which comprises these mutations in its epitope.

```
##The ddg calculation script is located at this folder;
./ddg_calcV3.py 
````
### Performing individual antibogy_GP complex mutation introduction and binding energy calculations in the mutant Vs wildType
#### GP_REGN3479
 - Description: Ebola GP with REGN3479 antibody (GP1: chains S,T,U; GP2: chains V,W,X; REGN3479 heavy: M,O,Q; REGN3479 light: N,P,R)
 - PDB: 7TN9_GP_REGN3479_with_GP2.pdb
 - Mutations on GP (chain V):I504T, V505L, N507S, N507R. These mutations are located on Chain V
```
PDB="/path_to_file/7TN9_GP_REGN3479_with_GP2.pdb"
# Run mutations on the GP (chain T) in the GP_REGN3479 complex
# New mutations from REGN3470 (assuming on chain T)
python3 ../ddg_calcV3.py -p "$PDB" -r 504 -c V -m T --relax_rounds 5 -o result_I504T_subset.csv
python3 ../ddg_calcV3.py -p "$PDB" -r 505 -c V -m L --relax_rounds 5 -o result_V505L_subset.csv
python3 ../ddg_calcV3.py -p "$PDB" -r 507 -c V -m S --relax_rounds 5 -o result_N507S_subset.csv
python3 ../ddg_calcV3.py -p "$PDB" -r 507 -c V -m R --relax_rounds 5 -o result_N507R_subset.csv
```
##### Running with trajectories and control
```
PDB="/path_to_file/7TN9_GP_REGN3479_with_GP2.pdb"

# Run mutations on GP (chain V) in the GP_REGN3479 complex
python3 ddg_calcV3_with_errors.py -p "$PDB" -r 504 -c V -m T --n_trajectories 35 --run_controls -o result_I504T_35traj.csv --output_controls controls_I504T.csv
python3 ddg_calcV3_with_errors.py -p "$PDB" -r 505 -c V -m L --n_trajectories 35 --run_controls -o result_V505L_35traj.csv --output_controls controls_V505L.csv
python3 ddg_calcV3_with_errors.py -p "$PDB" -r 507 -c V -m S --n_trajectories 35 --run_controls -o result_N507S_35traj.csv --output_controls controls_N507S.csv
python3 ddg_calcV3_with_errors.py -p "$PDB" -r 507 -c V -m R --n_trajectories 35 --run_controls -o result_N507R_35traj.csv --output_controls controls_N507R.csv

OR
sbatch GP_REGN3479_seed.sh
```
#### GP_REGN3470
 - Description: Ebola GP with REGN3470 antibody (GP1: chains S,T,U; GP2: chains V,W,X; REGN3470 heavy: G,I; REGN3470 light: H,J)
 - PDB: 7TN9_GP_REGN3470.pdb
 - Mutations on GP (chain T): S263N, K265R, E280T
```
PDB="/path/to/7TN9_GP_REGN3470.pdb"

python3 ddg_calcV3.py -p "$PDB" -r 263 -c T -m N --relax_rounds 5 -o result_S263N_subset.csv
python3 ddg_calcV3.py -p "$PDB" -r 265 -c T -m R --relax_rounds 5 -o result_K265R_subset.csv
python3 ddg_calcV3.py -p "$PDB" -r 280 -c T -m T --relax_rounds 5 -o result_E280T_subset.csv
```
##### With Trajectories and Control
```
PDB="/path/to/7TN9_GP_REGN3470.pdb"

python3 ddg_calcV3_with_errors.py -p "$PDB" -r 263 -c T -m N --n_trajectories 35 --run_controls -o result_S263N_35traj.csv --output_controls controls_S263N.csv
python3 ddg_calcV3_with_errors.py -p "$PDB" -r 265 -c T -m R --n_trajectories 35 --run_controls -o result_K265R_35traj.csv --output_controls controls_K265R.csv
python3 ddg_calcV3_with_errors.py -p "$PDB" -r 280 -c T -m T --n_trajectories 35 --run_controls -o result_E280T_35traj.csv --output_controls controls_E280T.csv
OR
sbatch GP_REGN3470_seed.sh
```
#### GP_REGN3471
 - Description: Ebola GP with REGN3471 antibody (GP1: chains S,T,U; GP2: chains V,W,X; REGN3471 heavy: A,C,E; REGN3471 light: B,D,F)
 - PDB: 7TN9_GP_REGN3471.pdb
 - Mutations on GP (chain T): E112D, P116A
```
PDB="/path/to/7TN9_GP_REGN3471.pdb"

python3 ddg_calcV3.py -p "$PDB" -r 112 -c T -m D --relax_rounds 5 -o result_E112D_subset.csv
python3 ddg_calcV3.py -p "$PDB" -r 116 -c T -m A --relax_rounds 5 -o result_P116A_subset.csv
```
##### With trajetories and controls (35 trajectories with error bars)
```
PDB="/path/to/7TN9_GP_REGN3471.pdb"

python3 ddg_calcV3_with_errors.py -p "$PDB" -r 112 -c T -m D --n_trajectories 35 --run_controls -o result_E112D_35traj.csv --output_controls controls_E112D.csv
python3 ddg_calcV3_with_errors.py -p "$PDB" -r 116 -c T -m A --n_trajectories 35 --run_controls -o result_P116A_35traj.csv --output_controls controls_P116A.csv
OR
sbatch GP_REGN3471_seed.sh
```
#### ddG Binding energy Results

## Full 7TN9 Complex Analysis (All mAbs Together)

### Background

While analyzing individual antibody-GP complexes provides valuable per-mAb binding information, the **full 7TN9 complex** (containing all three mAbs - REGN3470, REGN3471, and REGN3479 - bound to GP simultaneously) offers insights into:

1. **Crosstalk between antibodies**: How mutations affect binding when all mAbs are present
2. **Cooperativity**: Whether antibodies bind independently or cooperatively
3. **Overall complex stability**: The cumulative effect of BDBV mutations on the entire Inmazeb-GP complex

### Structure Overview

The full 7TN9 complex contains:
- **GP**: chains S, T, U, V, W, X (GP1 + GP2 trimer)
- **REGN3470**: chains G, H, I, J (heavy: G,I; light: H,J)
- **REGN3471**: chains A, B, C, D, E, F (heavy: A,C,E; light: B,D,F)
- **REGN3479**: chains M, N, O, P, Q, R (heavy: M,O,Q; light: N,P,R)

**Total**: ~3000+ residues

### Mutations Analyzed

#### GP1 Mutations (chain T)
| Residue | Wild Type | Mutant | Mutation Type |
|---------|-----------|--------|---------------|
| 112 | E (Glutamic Acid) | D (Aspartic Acid) | Conservative |
| 116 | P (Proline) | A (Alanine) | Conservative |
| 263 | S (Serine) | N (Asparagine) | Conservative |
| 265 | K (Lysine) | R (Arginine) | Conservative |

#### GP2 Mutations (chain V)
| Residue | Wild Type | Mutant | Mutation Type |
|---------|-----------|--------|---------------|
| 504 | I (Isoleucine) | T (Threonine) | Conservative |
| 505 | V (Valine) | L (Leucine) | Conservative |
| 507 | N (Asparagine) | T (Threonine) | Conservative |

### Single-Trajectory Calculations (Quick Testing)

```
PDB="/path/to/7TN9_clean.pdb"

#### GP1 mutations (chain T)
python3 ddg_calcV3_full.py -p "$PDB" -r 112 -c T -m D --relax_rounds 5 -o result_E112D_full.csv
python3 ddg_calcV3_full.py -p "$PDB" -r 116 -c T -m A --relax_rounds 5 -o result_P116A_full.csv
python3 ddg_calcV3_full.py -p "$PDB" -r 263 -c T -m N --relax_rounds 5 -o result_S263N_full.csv
python3 ddg_calcV3_full.py -p "$PDB" -r 265 -c T -m R --relax_rounds 5 -o result_K265R_full.csv
```

#### GP2 mutations (chain V)

```
python3 ddg_calcV3_full.py -p "$PDB" -r 504 -c V -m T --relax_rounds 5 -o result_I504T_full.csv
python3 ddg_calcV3_full.py -p "$PDB" -r 505 -c V -m L --relax_rounds 5 -o result_V505L_full.csv
python3 ddg_calcV3_full.py -p "$PDB" -r 507 -c V -m T --relax_rounds 5 -o result_N507T_full.csv
```
##### GP1 mutations (chain T) - 35 trajectories each
```
python3 ddg_calcV3_full_with_errors.py -p "$PDB" -r 112 -c T -m D --n_trajectories 35 --run_controls -o result_E112D_full_35traj.csv --output_controls controls_E112D_full.csv
python3 ddg_calcV3_full_with_errors.py -p "$PDB" -r 116 -c T -m A --n_trajectories 35 --run_controls -o result_P116A_full_35traj.csv --output_controls controls_P116A_full.csv
python3 ddg_calcV3_full_with_errors.py -p "$PDB" -r 263 -c T -m N --n_trajectories 35 --run_controls -o result_S263N_full_35traj.csv --output_controls controls_S263N_full.csv
python3 ddg_calcV3_full_with_errors.py -p "$PDB" -r 265 -c T -m R --n_trajectories 35 --run_controls -o result_K265R_full_35traj.csv --output_controls controls_K265R_full.csv
```
#### GP2 mutations (chain V) - 35 trajectories each

```
python3 ddg_calcV3_full_with_errors.py -p "$PDB" -r 504 -c V -m T --n_trajectories 35 --run_controls -o result_I504T_full_35traj.csv --output_controls controls_I504T_full.csv
python3 ddg_calcV3_full_with_errors.py -p "$PDB" -r 505 -c V -m L --n_trajectories 35 --run_controls -o result_V505L_full_35traj.csv --output_controls controls_V505L_full.csv
python3 ddg_calcV3_full_with_errors.py -p "$PDB" -r 507 -c V -m T --n_trajectories 35 --run_controls -o result_N507T_full_35traj.csv --output_controls controls_N507T_full.csv
```
## Citation
 - Rosetta/PyRosetta: Chaudhury, S., et al. (2010). PyRosetta: a script-based interface for implementing molecular modeling algorithms using Rosetta. Bioinformatics, 26(5), 689-691.

 - Flex ddG Protocol: Barlow, K. A., et al. (2018). Flex ddG: Rosetta ensemble-based estimation of changes in protein-protein binding affinity upon mutation. The Journal of Physical Chemistry B, 122(21), 5389-5399.

 - Inmazeb Structure: PDB 7TN9 - Cryo-EM structure of the Ebola virus glycoprotein in complex with the Inmazeb antibodies.
