#!/bin/bash
#SBATCH --job-name=FL_full_route
#SBATCH --time=02:00:00
#SBATCH --output=${FASTALAKE_DATA_DIR}/%x_%A_%a.log
set -euo pipefail
R=${FASTALAKE_DATA_DIR}/full_route_entrapment_2026-09-06
STAGE=${1:?stage}; COHORT_NAME=${2:?cohort}
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
ARGS=(--config "$R/experiment.json" "$STAGE" --cohort "$COHORT_NAME")
if [ "$STAGE" = sample ]; then ARGS+=(--index "${SLURM_ARRAY_TASK_ID:?array task}"); fi
${FASTALAKE_DATA_DIR}/python "$R/src/full_route_entrapment.py" "${ARGS[@]}"
