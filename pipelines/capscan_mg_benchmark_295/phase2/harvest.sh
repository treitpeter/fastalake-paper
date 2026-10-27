#!/bin/bash
cd LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase2
for u in figures/Fig4 figures/SFig06 figures/SFig09 text_statistics/supplementary_methods_6_7 text_statistics/methods text_statistics/abstract_discussion; do
  MPLBACKEND=Agg ${FASTALAKE_PYTHON:-python3} -B extract.py $1 $u 2>&1 | tail -n 2
done
