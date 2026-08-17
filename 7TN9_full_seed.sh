#!/bin/bash
#SBATCH --job-name=7TN9_full_ddg
#SBATCH --output=7TN9_full_%j.out
#SBATCH --error=7TN9_full_%j.err
#SBATCH --time=120:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --partition=shared

# ============================================
# 7TN9 Full Complex - All mAbs Together
# Using 7TN9_clean.pdb
# ============================================

PDB="/etc/ace-data/home/cknabukeera/pyrosetta/7TN9/7TN9_clean.pdb"
SCRIPT_DIR="/etc/ace-data/home/cknabukeera/pyrosetta/7TN9"
OUTPUT_DIR="/etc/ace-data/home/cknabukeera/pyrosetta/7TN9/7TN9_FULL"

mkdir -p $OUTPUT_DIR
cd $OUTPUT_DIR

echo "=========================================="
echo "7TN9 Full Complex - Rosetta ddG Calculations"
echo "PDB: $PDB"
echo "Output: $OUTPUT_DIR"
echo "Date: $(date)"
echo "=========================================="

echo ""
echo "⚠️ WARNING: Full complex ~3000+ residues"
echo "35 trajectories = 15-30 hours per mutation"
echo "=========================================="

# ============================================
# GP1 mutations (chain T)
# ============================================

echo ""
echo "=== GP1 mutations (chain T) - 35 trajectories ==="

python ${SCRIPT_DIR}/7TN9_ddg_full_with_errors.py -p "$PDB" -r 112 -c T -m D --n_trajectories 35 --run_controls -o result_E112D_7TN9_35traj.csv --output_controls controls_E112D_7TN9.csv && \
python ${SCRIPT_DIR}/7TN9_ddg_full_with_errors.py -p "$PDB" -r 116 -c T -m A --n_trajectories 35 --run_controls -o result_P116A_7TN9_35traj.csv --output_controls controls_P116A_7TN9.csv && \
python ${SCRIPT_DIR}/7TN9_ddg_full_with_errors.py -p "$PDB" -r 263 -c T -m N --n_trajectories 35 --run_controls -o result_S263N_7TN9_35traj.csv --output_controls controls_S263N_7TN9.csv && \
python ${SCRIPT_DIR}/7TN9_ddg_full_with_errors.py -p "$PDB" -r 265 -c T -m R --n_trajectories 35 --run_controls -o result_K265R_7TN9_35traj.csv --output_controls controls_K265R_7TN9.csv && \
python ${SCRIPT_DIR}/7TN9_ddg_full_with_errors.py -p "$PDB" -r 280 -c T -m T --n_trajectories 35 --run_controls -o result_E280T_7TN9_35traj.csv --output_controls controls_E280T_7TN9.csv

# ============================================
# GP2 mutations (chain V)
# ============================================

echo ""
echo "=== GP2 mutations (chain V) - 35 trajectories ==="

python ${SCRIPT_DIR}/7TN9_ddg_full_with_errors.py -p "$PDB" -r 504 -c V -m T --n_trajectories 35 --run_controls -o result_I504T_7TN9_35traj.csv --output_controls controls_I504T_7TN9.csv && \
python ${SCRIPT_DIR}/7TN9_ddg_full_with_errors.py -p "$PDB" -r 505 -c V -m L --n_trajectories 35 --run_controls -o result_V505L_7TN9_35traj.csv --output_controls controls_V505L_7TN9.csv && \
python ${SCRIPT_DIR}/7TN9_ddg_full_with_errors.py -p "$PDB" -r 507 -c V -m S --n_trajectories 35 --run_controls -o result_N507S_7TN9_35traj.csv --output_controls controls_N507S_7TN9.csv && \
python ${SCRIPT_DIR}/7TN9_ddg_full_with_errors.py -p "$PDB" -r 507 -c V -m R --n_trajectories 35 --run_controls -o result_N507R_7TN9_35traj.csv --output_controls controls_N507R_7TN9.csv

# ============================================
# COMBINE RESULTS
# ============================================

echo ""
echo "=== Combining results ==="

echo "mutation,resnum,chain,wt,mut,trajectories,mean_ddG,std_ddG,sem_ddG" > combined_7TN9_results.csv

for file in result_*_7TN9_35traj.csv; do
    if [ -f "$file" ]; then
        tail -n +2 "$file" >> combined_7TN9_results.csv
    fi
done

echo "=========================================="
echo "ALL 7TN9 MUTATIONS COMPLETED!"
echo "Date: $(date)"
echo "Results: $OUTPUT_DIR"
echo "=========================================="
ls -lh *.csv
EOF

chmod +x 7TN9_full_seed.sh
