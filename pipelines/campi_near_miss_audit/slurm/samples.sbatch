#!/bin/bash
set -euo pipefail
bash "${PIPELINE_DIR}/slurm/stage_frozen.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR=${FASTALAKE_OUT_DIR}/results/.matplotlib
${FASTALAKE_PYTHON:-python3} -B ${FASTALAKE_OUT_DIR}/frozen/analysis.py ${FASTALAKE_OUT_DIR} sample "${SLURM_ARRAY_TASK_ID}"
