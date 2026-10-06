#!/bin/bash
# A32 array task: two rank-5 Sage searches and the audit for one acquisition.
# Required environment: FASTALAKE_PIPELINES_DIR, FASTALAKE_OUT_DIR (holds inputs/manifest.json and inputs/samples_296.txt),
# FASTALAKE_SAGE (Sage 0.14.6 binary), FASTALAKE_PYTHON.
set -euo pipefail
: "${FASTALAKE_PIPELINES_DIR:?}" "${FASTALAKE_OUT_DIR:?}" "${FASTALAKE_SAGE:?}"
S=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$FASTALAKE_OUT_DIR/inputs/samples_296.txt")
[ -n "$S" ] || { echo "no sample at index $SLURM_ARRAY_TASK_ID"; exit 1; }
exec bash "$FASTALAKE_PIPELINES_DIR/predict_fl_only_spectrum_audit/task.sh" "$S"
