#!/bin/bash
set -euo pipefail
bash "${FASTALAKE_OUT_DIR}/frozen/slurm/stage_frozen.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
${FASTALAKE_PYTHON:-python3} -B ${FASTALAKE_OUT_DIR}/frozen/aggregate.py ${FASTALAKE_OUT_DIR} "${SLURM_ARRAY_TASK_ID}"
