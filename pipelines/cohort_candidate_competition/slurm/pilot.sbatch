#!/bin/bash
# Original: job 6025473, 1 CPU, 4G, 00:15:00.
#SBATCH --job-name=FL_review_candidates_pilot
#SBATCH --partition=p.hpcl8,p.hpcl92,p.hpcl93
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:15:00
set -euo pipefail
bash "${PIPELINE_DIR}/slurm/stage.sh"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
${FASTALAKE_PYTHON:-python3} "${FASTALAKE_OUT_DIR}/candidates/worker.py" pilot
