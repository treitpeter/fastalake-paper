#!/bin/bash
#SBATCH --job-name=mg295_db13
#SBATCH --partition=p.hpcl92
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1 --mem=4G --time=01:00:00
#SBATCH --output=LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/logs/db13_%j.out
#SBATCH --error=LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/logs/db13_%j.err
set -euo pipefail
cd LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/db
exec /mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python -B build_mg_fasta.py build13
