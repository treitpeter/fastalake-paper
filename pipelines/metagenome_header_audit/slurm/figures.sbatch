#!/bin/bash
# Original: job 6025428 (--dependency=afterok on both full jobs), 1 CPU, 20G, 1 h.
#SBATCH --job-name=FL_review_accum_figures
#SBATCH --partition=p.hpcl8,p.hpcl92,p.hpcl93
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1
#SBATCH --mem=20G
#SBATCH --time=01:00:00
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
${FASTALAKE_PYTHON:-python3} "${FASTALAKE_OUT_DIR}/plot.py"
