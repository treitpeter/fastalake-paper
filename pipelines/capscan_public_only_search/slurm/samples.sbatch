#!/bin/bash
# Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_WORK_DIR (working directory).
# Python code addresses inputs as LAB_ROOT/<relative path>, so link it into the working directory.
cd "${FASTALAKE_WORK_DIR:?set FASTALAKE_WORK_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=8
${FASTALAKE_PYTHON:-python3} -B ${FASTALAKE_LAB_ROOT}/MM_review_analysis/A08_4_5_capscan_public_only/frozen/analysis.py ${FASTALAKE_LAB_ROOT}/MM_review_analysis/A08_4_5_capscan_public_only sample "${SLURM_ARRAY_TASK_ID}"
