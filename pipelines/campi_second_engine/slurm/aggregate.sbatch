#!/bin/bash
set -euo pipefail
bash "${PIPELINE_DIR}/slurm/stage_frozen.sh"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
${FASTALAKE_PYTHON:-python3} ${FASTALAKE_OUT_DIR}/frozen/collect.py ${FASTALAKE_OUT_DIR} engine
