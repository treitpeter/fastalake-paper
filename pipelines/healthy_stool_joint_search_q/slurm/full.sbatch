#!/bin/bash
#SBATCH --job-name=FL_D3_full
#SBATCH --partition=p.hpcl8,p.hpcl93,p.hpcl92
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --array=0-2
#SBATCH --output=${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/slurm/%x_%A_%a.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ARGS=("original fileq ${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/out/original_fileq" "joint_corrected fileq ${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/out/joint_fileq" "joint_corrected sageq ${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/out/joint_sageq_repro")
/usr/bin/time -v ${FASTALAKE_DATA_DIR}/python -B ${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/code/extract_fileq.py ${ARGS[$SLURM_ARRAY_TASK_ID]}
