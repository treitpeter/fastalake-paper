#!/bin/bash
#SBATCH --cpus-per-task=1
set -euo pipefail
A=${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/reruns_347/A2_minlen7
/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python $A/scripts/collect.py
cat $A/RESULTS.md
