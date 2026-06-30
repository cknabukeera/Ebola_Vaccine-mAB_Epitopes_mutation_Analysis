#!/bin/bash
#SBATCH --job-name=ddg_7TN9
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --time=12:00:00
#SBATCH --partition=compute
#SBATCH --output=ddg_%j.out
#SBATCH --error=ddg_%j.err

export OMP_NUM_THREADS=${SLURM_CPUS_PER_TASK}
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

# Run mutations on the GP (chain T)
python3 ddg_calcV3.py -p 7TN9_clean.pdb -r 112 -c T -m D --relax_rounds 5 -o result_E112D.csv
python3 ddg_calcV3.py -p 7TN9_clean.pdb -r 116 -c T -m A --relax_rounds 5 -o result_P116A.csv
python3 ddg_calcV3.py -p 7TN9_clean.pdb -r 263 -c T -m N --relax_rounds 5 -o result_S263N.csv
python3 ddg_calcV3.py -p 7TN9_clean.pdb -r 265 -c T -m R --relax_rounds 5 -o result_K265R.csv

# Combine results
echo "resnum,chain,wt,mut,ddG,native_score,mutant_score" > all_results.csv
tail -n +2 result_*.csv >> all_results.csv
