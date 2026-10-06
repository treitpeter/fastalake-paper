#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
${FASTALAKE_DATA_DIR}/python -u ${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/hstool_precision.py "$1"
