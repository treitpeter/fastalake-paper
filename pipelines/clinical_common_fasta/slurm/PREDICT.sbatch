#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python -B ${FASTALAKE_DATA_DIR}/MM_review_analysis/end_of_day_2026-09-15/clinical_quant_review/code/summarize.py PREDICT
