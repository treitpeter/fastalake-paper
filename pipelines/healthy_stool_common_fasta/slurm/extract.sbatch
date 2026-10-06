#!/bin/bash
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python -B ${FASTALAKE_DATA_DIR}/MM_review_analysis/A21_7_hstool_quant_comparison/frozen/search.py ${FASTALAKE_DATA_DIR}/MM_review_analysis/A21_7_hstool_quant_comparison "${SLURM_ARRAY_TASK_ID}"
