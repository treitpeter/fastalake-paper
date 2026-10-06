#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR=${FASTALAKE_DATA_DIR}/MM_review_analysis/A01_4_1_database_sizes/results/.matplotlib
/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python -B ${FASTALAKE_DATA_DIR}/MM_review_analysis/A01_4_1_database_sizes/frozen/audit_v2.py ${FASTALAKE_DATA_DIR}/MM_review_analysis/A01_4_1_database_sizes
