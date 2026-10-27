#!/bin/bash
#SBATCH --job-name=P3i4_collect
#SBATCH --partition=p.hpcl92
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=00:30:00
#SBATCH --output=LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item4_pep/logs/collect_%j.out
#SBATCH --error=LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item4_pep/logs/collect_%j.err
# --mem: controls 6085691 MaxRSS 146 MB (analyze) x1.05 -> 1G floor
set -uo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
PY=${FASTALAKE_PYTHON:-python3}
cd LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item4_pep
$PY conf_summary.py negctrl
$PY analyze.py n295
