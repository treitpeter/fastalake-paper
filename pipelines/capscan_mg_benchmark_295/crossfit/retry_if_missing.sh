#!/bin/bash
# OOM safety net: rerun a cross-fit task only if its COMPLETE.json is missing (same code, more memory)
I=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/crossfit/design/ids_p15_rest.txt)
S=$(awk -F'\t' -v i=$I '$1==i{print $3}' LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/disjoint/roster_mg295.tsv)
[ -f LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/crossfit/outputs/CAPSCAN/$S/COMPLETE.json ] && { echo "SKIP $S complete"; exit 0; }
exec bash LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/crossfit/design/run.sbatch crossfit $I
