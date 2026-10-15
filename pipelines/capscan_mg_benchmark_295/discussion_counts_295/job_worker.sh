#!/bin/bash
# usage: job_worker.sh <tree> <accum|cand>   -- runs the byte-identical archived worker inside a mirror tree
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
S=LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/discussion_counts_295/trees/$1/review/v15_integration/analyses/sequence_accumulation_2026-09-29
sha256sum $S/worker.py $S/candidates/worker.py $S/manifest.json
if [ "$2" = accum ]; then /mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python $S/worker.py CAPSCAN full; else /mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python $S/candidates/worker.py full; fi
