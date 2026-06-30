## Inmazeb individual component analysis

Since this is a complex of three antibodies, I need to extract the individual GP proteins and its corresponding antibody, to have the GP interacting individually with each antibody

### 1. Using PyMOL to extract the GP_REGN3479 chains
#### Step 1: Loading the cleaned PDB structure:7TN9

```
#Loading the PDB structure
load 7TN9_clean.pdb
```
##### Creating the GP_REGN3479
```
# select the GP
select GP, chain S+T+U
# The GP chains are S, T, U

# Select the REGN3479
select REGN3479, chain M+O+Q+N+P+R
#Heavy chains: M, O, Q
#Light chains: N, P, R

#Select the combined GP_REGN3479
select GP_REGN3479, chain S+T+U+M+O+Q+N+P+R
```

```
#View the Selection Individually
hide everything
show cartoon, GP_REGN3479
color blue, chain S+T+U
color red, chain M+O+Q
color yellow, chain N+P+R
```
```
# Extract the Selection into a New Object
extract complex_obj, GP_REGN3479
```

```
# Save the Extracted PDB
save 7TN9_GP_REGN3479.pdb, complex_obj
```
```
#save the GP and REGN3479 differently
extract GP_obj, GP
save 7TN9_GP_only.pdb, GP_obj
extract REGN3479_obj, REGN3479
save 7TN9_REGN3479.pdb, REGN3479_obj
```

### 2. In silico mutagenesis(introducing BDBV mutations in Zaire) and ddG energy claculations with pyrosetta
 #### ddG script
```
#ddg calculation script

```
``
python3 ../ddg_calcV3.py -p "$PDB" -r 504 -c T -m T --relax_rounds 5 -o result_I504T_subset.csv
python3 ../ddg_calcV3.py -p "$PDB" -r 505 -c T -m L --relax_rounds 5 -o result_V505L_subset.csv
python3 ../ddg_calcV3.py -p "$PDB" -r 507 -c T -m S --relax_rounds 5 -o result_N507S_subset.csv
```

python3 ../ddg_calcV3.py -p "$PDB" -r 507 -c T -m R --relax_rounds 5 -o result_N507R_subset.csv

```
