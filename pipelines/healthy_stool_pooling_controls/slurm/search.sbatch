#!/bin/bash
#SBATCH --job-name=fig3_pool_controls
#SBATCH --partition=p.hpcl93
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --array=0-23%6
#SBATCH --output=${FASTALAKE_DATA_DIR}/%x_%A_%a.out
#SBATCH --error=${FASTALAKE_DATA_DIR}/%x_%A_%a.err
set -euo pipefail
export RAYON_NUM_THREADS="$SLURM_CPUS_PER_TASK"
export RUST_MIN_STACK=16777216
ulimit -s unlimited
${FASTALAKE_DATA_DIR}/python ${FASTALAKE_DATA_DIR}/run_task.py "$SLURM_ARRAY_TASK_ID"
