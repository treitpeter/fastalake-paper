#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR=${FASTALAKE_DATA_DIR}/MM_review_analysis/A02_4_2_lookup_specificity/results/.matplotlib
/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python -B ${FASTALAKE_DATA_DIR}/MM_review_analysis/A02_4_2_lookup_specificity/frozen/worker.py ${FASTALAKE_DATA_DIR}/MM_review_analysis/A02_4_2_lookup_specificity lookup
