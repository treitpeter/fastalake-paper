#!/bin/bash
#SBATCH -J i3_build -p p.hpcl92 --nice=10000 -c 1 --mem=2G -t 0:30:00
#SBATCH -o LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item3_collision/logs/build_%j.log
set -euo pipefail
/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item3_collision/build_manifest.py
