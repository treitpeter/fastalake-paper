#!/bin/bash
set -euo pipefail
export PYTHONNOUSERSITE=1
export RAYON_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
TASK_MODE=$1
TASK_WORK=$2
exec /mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python "$TASK_WORK/pipeline.py" "$TASK_MODE" "${SLURM_ARRAY_TASK_ID:-0}"
