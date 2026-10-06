#!/bin/bash
# Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_WORK_DIR (working directory).
# Python code addresses inputs as LAB_ROOT/<relative path>, so link it into the working directory.
cd "${FASTALAKE_WORK_DIR:?set FASTALAKE_WORK_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python -B ${FASTALAKE_LAB_ROOT}/MM_review_analysis/A16_4_5_public_comparison/frozen/search.py ${FASTALAKE_LAB_ROOT}/MM_review_analysis/A16_4_5_public_comparison "${SLURM_ARRAY_TASK_ID}"
