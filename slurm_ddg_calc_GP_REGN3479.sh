#!/bin/bash
#SBATCH --job-name=ddg_GP_REGN3479
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=8:00:00
#SBATCH --partition=compute
#SBATCH --output=subset_%j.out
#SBATCH --error=subset_%j.err

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

mkdir -p GP_REGN3479
cd GP_REGN3479 || exit 1
PDB="/etc/ace-data/home/cknabukeera/pyrosetta/7TN9/7TN9_GP_REGN3479_with_GP2.pdb"
# Run mutations on the GP (chain T) in the GP_REGN3479 complex
# New mutations from REGN3470 (assuming on chain T)
python3 ../ddg_calcV3.py -p "$PDB" -r 504 -c V -m T --relax_rounds 5 -o result_I504T_subset.csv
python3 ../ddg_calcV3.py -p "$PDB" -r 505 -c V -m L --relax_rounds 5 -o result_V505L_subset.csv
python3 ../ddg_calcV3.py -p "$PDB" -r 507 -c V -m S --relax_rounds 5 -o result_N507S_subset.csv
python3 ../ddg_calcV3.py -p "$PDB" -r 507 -c V -m R --relax_rounds 5 -o result_N507R_subset.csv


# Combine subset results
echo "resnum,chain,wt,mut,ddG,native_score,mutant_score" > all_subset_results.csv
tail -n +2 result_*_subset.csv >> all_subset_results.csv
