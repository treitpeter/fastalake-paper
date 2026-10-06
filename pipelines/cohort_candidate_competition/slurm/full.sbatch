#!/bin/bash
# Original: job 6025474, 1 CPU, 6G, 03:00:00.
#SBATCH --job-name=FL_review_candidates_full
#SBATCH --partition=p.hpcl8,p.hpcl92,p.hpcl93
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1
#SBATCH --mem=6G
#SBATCH --time=03:00:00
set -euo pipefail
bash "${PIPELINE_DIR}/slurm/stage.sh"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
${FASTALAKE_PYTHON:-python3} "${FASTALAKE_OUT_DIR}/candidates/full_gate.py" 
