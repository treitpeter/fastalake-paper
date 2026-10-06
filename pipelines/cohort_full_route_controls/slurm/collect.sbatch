#!/bin/bash
#SBATCH --job-name=FLfull_collect
#SBATCH --partition=p.hpcl93
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=${FASTALAKE_DATA_DIR}/collect_%j.log
set -euo pipefail
export OPENBLAS_NUM_THREADS=1
${FASTALAKE_DATA_DIR}/python ${FASTALAKE_DATA_DIR}/collect_cohort.py "${1:?cohort}"
