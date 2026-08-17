#!/bin/bash
#SBATCH --job-name=REGN3470_ddg
#SBATCH --output=REGN3470_%j.out
#SBATCH --error=REGN3470_%j.err
#SBATCH --time=48:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=shared

# ============================================
# GP_REGN3470 - Ebola GP with REGN3470 antibody
# ============================================
PDB="/etc/ace-data/home/cknabukeera/pyrosetta/7TN9/7TN9_GP_REGN3470.pdb"
SCRIPT_DIR="/etc/ace-data/home/cknabukeera/pyrosetta/7TN9"
OUTPUT_DIR="/etc/ace-data/home/cknabukeera/pyrosetta/7TN9/GP_REGN3470/REGN3470_seeds"

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

echo "=========================================="
echo "GP_REGN3470 - Rosetta ddG Calculations"
echo "PDB: $PDB"
echo "Output: $OUTPUT_DIR"
echo "Date: $(date)"
echo "=========================================="

# ============================================
# RUN MUTATIONS WITH 35 TRAJECTORIES
# Mutations: S263N, K265R, E280T on chain T
# ============================================

echo ""
echo "=== Running S263N (35 trajectories) ==="
python ${SCRIPT_DIR}/ddg_calcV3_with_errors.py \
    -p "$PDB" -r 263 -c T -m N \
    --n_trajectories 35 --run_controls \
    -o result_S263N_35traj.csv \
    --output_controls controls_S263N.csv

if [ $? -ne 0 ]; then
    echo "❌ ERROR: S263N failed"
    exit 1
fi
echo "✅ S263N completed"

echo ""
echo "=== Running K265R (35 trajectories) ==="
python ${SCRIPT_DIR}/ddg_calcV3_with_errors.py \
    -p "$PDB" -r 265 -c T -m R \
    --n_trajectories 35 --run_controls \
    -o result_K265R_35traj.csv \
    --output_controls controls_K265R.csv

if [ $? -ne 0 ]; then
    echo "❌ ERROR: K265R failed"
    exit 1
fi
echo "✅ K265R completed"

echo ""
echo "=== Running E280T (35 trajectories) ==="
python ${SCRIPT_DIR}/ddg_calcV3_with_errors.py \
    -p "$PDB" -r 280 -c T -m T \
    --n_trajectories 35 --run_controls \
    -o result_E280T_35traj.csv \
    --output_controls controls_E280T.csv

if [ $? -ne 0 ]; then
    echo "❌ ERROR: E280T failed"
    exit 1
fi
echo "✅ E280T completed"

# ============================================
# ALL COMPLETED SUCCESSFULLY
# ============================================

echo ""
echo "=========================================="
echo "ALL GP_REGN3470 MUTATIONS COMPLETED!"
echo "Date: $(date)"
echo "Results saved in: $OUTPUT_DIR"
echo "=========================================="
echo ""
echo "Results files:"
ls -lh *.csv
