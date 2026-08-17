## Inmazeb individual component analysis

Since this is a complex of three antibodies, I need to extract the individual GP proteins and its corresponding antibody, to have the GP interacting individually with each antibody
# Antibody-GP Complex Binding Energy Analysis

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
# GP2 mutations (chain V) - 35 trajectories each

```
python3 ddg_calcV3_full_with_errors.py -p "$PDB" -r 504 -c V -m T --n_trajectories 35 --run_controls -o result_I504T_full_35traj.csv --output_controls controls_I504T_full.csv
python3 ddg_calcV3_full_with_errors.py -p "$PDB" -r 505 -c V -m L --n_trajectories 35 --run_controls -o result_V505L_full_35traj.csv --output_controls controls_V505L_full.csv
python3 ddg_calcV3_full_with_errors.py -p "$PDB" -r 507 -c V -m T --n_trajectories 35 --run_controls -o result_N507T_full_35traj.csv --output_controls controls_N507T_full.csv
```
