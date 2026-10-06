#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=2
${FASTALAKE_PYTHON:-python3} ${PIPELINE_DIR}/search.py pilot
${FASTALAKE_PYTHON:-python3} ${PIPELINE_DIR}/collect.py pilot
