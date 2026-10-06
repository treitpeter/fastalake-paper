#!/bin/bash
#SBATCH -J C3_region_test
#SBATCH -p p.hpcl92
#SBATCH -c 1
#SBATCH --mem=4G
#SBATCH -t 00:30:00
#SBATCH -o slurm-%x-%j.out
# Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_WORK_DIR (working directory).
# Python code addresses inputs as LAB_ROOT/<relative path>, so link it into the working directory.
cd "${FASTALAKE_WORK_DIR:?set FASTALAKE_WORK_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
set -euo pipefail
${FASTALAKE_PYTHON:-python} ${FASTALAKE_LAB_ROOT}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C3_region_test/region_test.py ${FASTALAKE_LAB_ROOT}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C3_region_test/../../../iteration4/crossfit/run/PER_ACQUISITION.tsv ${FASTALAKE_LAB_ROOT}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C3_region_test > ${FASTALAKE_LAB_ROOT}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C3_region_test/STDOUT.txt 2>&1
