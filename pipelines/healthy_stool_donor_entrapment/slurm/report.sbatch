#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python ${FASTALAKE_DATA_DIR}/MM_review_analysis/A11_4_2_control_reconciliation/frozen/analysis.py ${FASTALAKE_DATA_DIR}/MM_review_analysis/A11_4_2_control_reconciliation collect
