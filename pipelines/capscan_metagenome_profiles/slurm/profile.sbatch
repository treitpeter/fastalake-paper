#!/bin/bash -l
#SBATCH --job-name=capscan_dna
#SBATCH --partition=p.hpcl92
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=04:00:00
#SBATCH --nice=10000
#SBATCH -o slurm-%x-%j.out
# Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_WORK_DIR (working directory).
# Python code addresses inputs as LAB_ROOT/<relative path>, so link it into the working directory.
cd "${FASTALAKE_WORK_DIR:?set FASTALAKE_WORK_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
set -euo pipefail
source analyses/capscan_metag_analysis_2026-09-10/environment.sh
${FASTALAKE_PYTHON:-python} -u analyses/capscan_metag_analysis_2026-09-10/run_sample.py --index "$SLURM_ARRAY_TASK_ID" --stage "${1:-all}" --threads "$SLURM_CPUS_PER_TASK"
