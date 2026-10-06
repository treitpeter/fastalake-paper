#!/bin/bash
# A25 case audit: usage  sbatch slurm/run_task.sbatch pilot|full
# Required: FASTALAKE_PIPELINES_DIR, FASTALAKE_DATA_DIR (lab project root), FASTALAKE_PYTHON.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
: "${FASTALAKE_PIPELINES_DIR:?}" "${FASTALAKE_DATA_DIR:?}"
export MM_REVIEW_ROOT="$FASTALAKE_DATA_DIR/MM_review_analysis"
exec "${FASTALAKE_PYTHON:-python3}" -B "$FASTALAKE_PIPELINES_DIR/predict_histidine_case_audit/audit_predict_spectrum_cases.py" --phase "${1:?pilot or full}"
