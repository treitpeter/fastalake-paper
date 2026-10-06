#!/bin/bash
# usage: sbatch [opts] run.sbatch <mode> <index or id-list file in design/>
# Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_WORK_DIR (working directory).
# Python code addresses inputs as LAB_ROOT/<relative path>, so link it into the working directory.
cd "${FASTALAKE_WORK_DIR:?set FASTALAKE_WORK_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
set -euo pipefail
export PYTHONNOUSERSITE=1
export RAYON_NUM_THREADS=${SLURM_CPUS_PER_TASK:-8}
cd ${FASTALAKE_LAB_ROOT}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration4/crossfit/design
exec ${FASTALAKE_PYTHON:-python3} worker.py "$1" "$2"
