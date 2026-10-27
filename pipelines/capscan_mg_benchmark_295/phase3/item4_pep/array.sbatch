#!/bin/bash
#SBATCH --job-name=P3i4_pep13
#SBATCH --partition=p.hpcl92
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1
#SBATCH --mem=800M
#SBATCH --time=00:30:00
#SBATCH --array=1-12
#SBATCH --output=LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item4_pep/logs/array_%A_%a.out
#SBATCH --error=LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item4_pep/logs/array_%A_%a.err
# --mem: pilot 6085666 MaxRSS 696,280 kB x1.05 = 731 MB -> 800M
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
/usr/bin/time -v ${FASTALAKE_PYTHON:-python3} LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item4_pep/run_pep.py new $SLURM_ARRAY_TASK_ID
