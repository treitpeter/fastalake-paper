#!/bin/bash
# runs the production second_engine_all/code/task.sh (unchanged) on the mg295 RUN dir; index = SLURM_ARRAY_TASK_ID or $1
exec LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/code/task.sh LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/second_engine "${1:-${SLURM_ARRAY_TASK_ID}}"
