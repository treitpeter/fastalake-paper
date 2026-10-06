#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
cd "${FASTALAKE_CODE_DIR}"
git rev-parse HEAD
date -Is
exec ${FASTALAKE_PYTHON:-python3} ${PIPELINE_DIR}/collect.py --manifest "$1"
