#!/bin/bash
# usage: job_producers.sh <c282|r295>  -- runs a_accumulation.py and e_candidates.py e1 (byte copies) with common.py ACC -> that tree
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cd LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/discussion_counts_295
cp common_$1.py common.py; rm -rf __pycache__
sha256sum common.py a_accumulation.py e_candidates.py
/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python a_accumulation.py
/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python e_candidates.py e1
mkdir -p out_$1
mv A_accumulation_readout.tsv A_per_specimen_ids.tsv A_per_specimen_medians.tsv A_RESULT.json E1_per_acquisition.tsv E1_summary.tsv E1_RESULT.json out_$1/
cp common_original.py common.py; rm -rf __pycache__
