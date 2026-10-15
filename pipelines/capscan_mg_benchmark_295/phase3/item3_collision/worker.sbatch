#!/bin/bash
# usage: sbatch -J i3_<run> --mem=<M> -o logs/<run>_%j.log worker.sbatch <run282|run295|neg295>
#SBATCH -p p.hpcl92 --nice=10000 -c 1 -t 3:00:00
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
A=LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item3_collision/lab/review/v15_integration/analyses
/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python "$A/$1/worker.py" CAPSCAN full
