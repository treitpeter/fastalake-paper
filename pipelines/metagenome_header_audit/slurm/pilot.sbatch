#!/bin/bash
# Original: job 6025421 (pilot, three specimens per cohort), 1 CPU, 12G, 20 min.
#SBATCH --job-name=FL_review_accum_pilot
#SBATCH --partition=p.hpcl8,p.hpcl92,p.hpcl93
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1
#SBATCH --mem=12G
#SBATCH --time=00:20:00
set -euo pipefail
bash "${PIPELINE_DIR}/slurm/stage.sh"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
${FASTALAKE_PYTHON:-python3} "${FASTALAKE_OUT_DIR}/worker.py" PREDICT pilot
${FASTALAKE_PYTHON:-python3} "${FASTALAKE_OUT_DIR}/worker.py" CAPSCAN pilot
