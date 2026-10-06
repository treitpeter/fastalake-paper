#!/bin/bash
# Original: job 6025427, partitions p.hpcl8,p.hpcl92,p.hpcl93, nice 10000, 1 CPU, 48G request, time 06:00:00.
#SBATCH --job-name=FL_review_accum_full_CAPSCAN
#SBATCH --partition=p.hpcl8,p.hpcl92,p.hpcl93
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --time=06:00:00
set -euo pipefail
bash "${PIPELINE_DIR}/slurm/stage.sh"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
${FASTALAKE_PYTHON:-python3} "${FASTALAKE_OUT_DIR}/worker.py" CAPSCAN full
