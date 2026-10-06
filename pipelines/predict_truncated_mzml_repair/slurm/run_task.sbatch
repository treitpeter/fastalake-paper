#!/bin/bash
# A25 repair: search one repaired PREDICT acquisition (array index into raw_recovery/manifest.json tasks).
# Required: FASTALAKE_PIPELINES_DIR, FASTALAKE_DATA_DIR (lab project root holding the recovery directories), FASTALAKE_PYTHON.
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=8
: "${FASTALAKE_PIPELINES_DIR:?}" "${FASTALAKE_DATA_DIR:?}"
exec "${FASTALAKE_PYTHON:-python3}" -B -u "$FASTALAKE_PIPELINES_DIR/predict_truncated_mzml_repair/repaired_predict_searches_20260915.py" "${SLURM_ARRAY_TASK_ID}"
