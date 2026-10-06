#!/bin/bash
#SBATCH --partition=p.hpcl8,p.hpcl92,p.hpcl93
#SBATCH --nice=10000
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --job-name=FL_cf_final_verify_analyze
# Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_WORK_DIR (working directory).
# Python code addresses inputs as LAB_ROOT/<relative path>, so link it into the working directory.
cd "${FASTALAKE_WORK_DIR:?set FASTALAKE_WORK_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
set -euo pipefail
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
cd ${FASTALAKE_LAB_ROOT}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration4/crossfit/run
${FASTALAKE_PYTHON:-python} verify.py
${FASTALAKE_PYTHON:-python} analyze_final.py > ANALYZE_STDOUT.txt
