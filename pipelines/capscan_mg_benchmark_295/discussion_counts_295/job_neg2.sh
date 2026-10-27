#!/bin/bash
set -uo pipefail
cd LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/discussion_counts_295
/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python build_neg2.py || exit 2
bash job_worker.sh neg2 accum; rc=$?
echo "worker exit code $rc (expected nonzero)"
if [ $rc -ne 0 ]; then echo NEGATIVE_CONTROL_PASS; exit 0; else echo NEGATIVE_CONTROL_FAIL; exit 1; fi
