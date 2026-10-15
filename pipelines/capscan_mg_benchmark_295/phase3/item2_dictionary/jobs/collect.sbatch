#!/bin/bash
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
exec /mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item2_dictionary/code/dispatch.py collect "$1/manifest.json"
