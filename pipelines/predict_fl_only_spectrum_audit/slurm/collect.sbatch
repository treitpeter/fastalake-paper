#!/bin/bash
# Summary job of A32 (8 GB, 1 s in the original).
set -euo pipefail
: "${FASTALAKE_PIPELINES_DIR:?}" "${FASTALAKE_OUT_DIR:?}"
exec "${FASTALAKE_PYTHON:-python3}" "$FASTALAKE_PIPELINES_DIR/predict_fl_only_spectrum_audit/collect.py"
