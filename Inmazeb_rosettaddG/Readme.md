## Inmazeb individual component analysis

Since this is a complex of three antibodies, I need to extract the individual GP proteins and its corresponding antibody, to have the GP interacting individually with each antibody

### 1. Using PyMOL to extract the GP_REGN3479 chains

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


### 2. In silico mutagenesis(introducing BDBV mutations in Zaire) and ddG energy claculations with pyrosetta
 #### ddG script
```
#ddg calculation script
```

```
PDB="/path_to_file/7TN9_GP_REGN3479_with_GP2.pdb"
# Run mutations on the GP (chain T) in the GP_REGN3479 complex
# New mutations from REGN3470 (assuming on chain T)
python3 ../ddg_calcV3.py -p "$PDB" -r 504 -c V -m T --relax_rounds 5 -o result_I504T_subset.csv
python3 ../ddg_calcV3.py -p "$PDB" -r 505 -c V -m L --relax_rounds 5 -o result_V505L_subset.csv
python3 ../ddg_calcV3.py -p "$PDB" -r 507 -c V -m S --relax_rounds 5 -o result_N507S_subset.csv
python3 ../ddg_calcV3.py -p "$PDB" -r 507 -c V -m R --relax_rounds 5 -o result_N507R_subset.csv
```


