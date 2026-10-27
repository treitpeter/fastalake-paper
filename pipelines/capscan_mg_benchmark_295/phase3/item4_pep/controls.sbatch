#!/bin/bash
#SBATCH --job-name=P3i4_ctrl
#SBATCH --partition=p.hpcl92
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item4_pep/logs/controls_%j.out
#SBATCH --error=LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item4_pep/logs/controls_%j.err
set -uo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
PY=${FASTALAKE_PYTHON:-python3}
cd LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item4_pep
/usr/bin/time -v $PY conf_summary.py control
/usr/bin/time -v $PY conf_summary.py n295
/usr/bin/time -v $PY conf_summary.py negctrl
/usr/bin/time -v $PY analyze.py control
for f in CF_PEP_PER_ACQUISITION.tsv CF_PEP_SUMMARY.tsv CF_PEP_SUMMARY_BY_REGION.tsv; do
  if cmp -s control/$f LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C1_pep_primary/crossfit/$f; then echo "ANALYZE_CONTROL $f IDENTICAL"; else echo "ANALYZE_CONTROL $f DIFFERENT"; fi
done
