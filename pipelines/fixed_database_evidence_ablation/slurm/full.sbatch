#!/bin/bash
# Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_WORK_DIR (working directory).
# Python code addresses inputs as LAB_ROOT/<relative path>, so link it into the working directory.
cd "${FASTALAKE_WORK_DIR:?set FASTALAKE_WORK_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python -B ${FASTALAKE_LAB_ROOT}/MM_review_analysis/end_of_day_2026-09-15/evidence_ablation/code/partition.py ${FASTALAKE_LAB_ROOT}/MM_review_analysis/end_of_day_2026-09-15/evidence_ablation full 0
