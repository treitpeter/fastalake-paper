#!/bin/bash
# args: <run dir> <n> <arm>; export_assignments.py copy (diff: EXPORT_CHANGES.diff = output/input locations only); production interpreter
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export P3_EXPORT_OUT="$1/export" P3_GROUPING_RUN="$1"
exec ${FASTALAKE_PYTHON:-python3} LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item2_dictionary/code/export_assignments.py "CAPSCAN_$2_$3"
