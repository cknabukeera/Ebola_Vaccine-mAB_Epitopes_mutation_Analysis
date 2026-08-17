#!/bin/bash
#SBATCH --job-name=REGN3471_ddg
#SBATCH --output=REGN3471_%j.out
#SBATCH --error=REGN3471_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=shared

# ============================================
# GP_REGN3471 - Ebola GP with REGN3471 antibody
# Mutations: E112D, P116A on chain T
# ============================================
PDB="/etc/ace-data/home/cknabukeera/pyrosetta/7TN9/7TN9_GP_REGN3471.pdb"
SCRIPT_DIR="/etc/ace-data/home/cknabukeera/pyrosetta/7TN9"
OUTPUT_DIR="/etc/ace-data/home/cknabukeera/pyrosetta/7TN9/GP_REGN3471/REGN3471_seeds"

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

echo "=========================================="
echo "GP_REGN3471 - Rosetta ddG Calculations"
echo "PDB: $PDB"
echo "Output: $OUTPUT_DIR"
echo "Date: $(date)"
echo "=========================================="

# ============================================
# RUN MUTATIONS WITH 35 TRAJECTORIES
# Mutations: E112D, P116A on chain T
# ============================================

echo ""
echo "=== Running E112D (35 trajectories) ==="
python ${SCRIPT_DIR}/ddg_calcV3_with_errors.py \
    -p "$PDB" -r 112 -c T -m D \
    --n_trajectories 35 --run_controls \
    -o result_E112D_35traj.csv \
    --output_controls controls_E112D.csv

if [ $? -ne 0 ]; then
    echo "❌ ERROR: E112D failed"
    exit 1
fi
echo "✅ E112D completed"

echo ""
echo "=== Running P116A (35 trajectories) ==="
python ${SCRIPT_DIR}/ddg_calcV3_with_errors.py \
    -p "$PDB" -r 116 -c T -m A \
    --n_trajectories 35 --run_controls \
    -o result_P116A_35traj.csv \
    --output_controls controls_P116A.csv

if [ $? -ne 0 ]; then
    echo "❌ ERROR: P116A failed"
    exit 1
fi
echo "✅ P116A completed"

# ============================================
# COMBINE RESULTS
# ============================================

echo ""
echo "=== Combining results ==="
echo "mutation,resnum,chain,wt,mut,trajectories,mean_ddG,std_ddG,sem_ddG" > combined_results.csv
for file in result_*_35traj.csv; do
    if [ -f "$file" ]; then
        tail -n +2 "$file" >> combined_results.csv
    fi
done

echo "Combined results saved to: combined_results.csv"

# ============================================
# ALL COMPLETED SUCCESSFULLY
# ============================================

echo ""
echo "=========================================="
echo "ALL GP_REGN3471 MUTATIONS COMPLETED!"
echo "Date: $(date)"
echo "Results saved in: $OUTPUT_DIR"
echo "=========================================="
echo ""
echo "Results files:"
ls -lh *.csv
