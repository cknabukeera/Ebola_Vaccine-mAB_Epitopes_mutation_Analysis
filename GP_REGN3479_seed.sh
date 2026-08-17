#!/bin/bash
#SBATCH --job-name=REGN3479_ddg
#SBATCH --output=REGN3479_%j.out
#SBATCH --error=REGN3479_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=shared

# ============================================
# GP_REGN3479 - Ebola GP with REGN3479 antibody
# ============================================
PDB="/etc/ace-data/home/cknabukeera/pyrosetta/7TN9/7TN9_GP_REGN3479_with_GP2.pdb"
SCRIPT_DIR="/etc/ace-data/home/cknabukeera/pyrosetta/7TN9"
OUTPUT_DIR="/etc/ace-data/home/cknabukeera/pyrosetta/7TN9/GP_REGN3479/REGN3479_seeds"

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

echo "=========================================="
echo "GP_REGN3479 - Rosetta ddG Calculations"
echo "PDB: $PDB"
echo "Output: $OUTPUT_DIR"
echo "Date: $(date)"
echo "=========================================="

# ============================================
# RUN MUTATIONS WITH 35 TRAJECTORIES
# Mutations: I504T, V505L, N507S, N507R on chain V
# ============================================

echo ""
echo "=== Running I504T (35 trajectories) ==="
python ${SCRIPT_DIR}/ddg_calcV3_with_errors.py \
    -p "$PDB" -r 504 -c V -m T \
    --n_trajectories 35 --run_controls \
    -o result_I504T_35traj.csv \
    --output_controls controls_I504T.csv

if [ $? -ne 0 ]; then
    echo "❌ ERROR: I504T failed"
    exit 1
fi
echo "✅ I504T completed"

echo ""
echo "=== Running V505L (35 trajectories) ==="
python ${SCRIPT_DIR}/ddg_calcV3_with_errors.py \
    -p "$PDB" -r 505 -c V -m L \
    --n_trajectories 35 --run_controls \
    -o result_V505L_35traj.csv \
    --output_controls controls_V505L.csv

if [ $? -ne 0 ]; then
    echo "❌ ERROR: V505L failed"
    exit 1
fi
echo "✅ V505L completed"

echo ""
echo "=== Running N507S (35 trajectories) ==="
python ${SCRIPT_DIR}/ddg_calcV3_with_errors.py \
    -p "$PDB" -r 507 -c V -m S \
    --n_trajectories 35 --run_controls \
    -o result_N507S_35traj.csv \
    --output_controls controls_N507S.csv

if [ $? -ne 0 ]; then
    echo "❌ ERROR: N507S failed"
    exit 1
fi
echo "✅ N507S completed"

echo ""
echo "=== Running N507R (35 trajectories) ==="
python ${SCRIPT_DIR}/ddg_calcV3_with_errors.py \
    -p "$PDB" -r 507 -c V -m R \
    --n_trajectories 35 --run_controls \
    -o result_N507R_35traj.csv \
    --output_controls controls_N507R.csv

if [ $? -ne 0 ]; then
    echo "❌ ERROR: N507R failed"
    exit 1
fi
echo "✅ N507R completed"

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
echo "ALL GP_REGN3479 MUTATIONS COMPLETED!"
echo "Date: $(date)"
echo "Results saved in: $OUTPUT_DIR"
echo "=========================================="
echo ""
echo "Results files:"
ls -lh *.csv
