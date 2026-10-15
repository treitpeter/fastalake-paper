#!/bin/bash
# usage: sbatch [opts] run.sbatch <mode> <index or id-list file in design/>
set -euo pipefail
export PYTHONNOUSERSITE=1
export RAYON_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
cd LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/crossfit/design
exec /mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python worker.py "$1" "$2"
