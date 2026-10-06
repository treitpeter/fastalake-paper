#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
${FASTALAKE_PYTHON:-python3} -B ${PIPELINE_DIR}/collect_all.py ${FASTALAKE_OUT_DIR}
