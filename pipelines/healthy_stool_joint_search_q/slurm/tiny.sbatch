#!/bin/bash
#SBATCH --job-name=FL_D3_tiny
#SBATCH --partition=p.hpcl8,p.hpcl93,p.hpcl92
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/slurm/%x_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
P="${FASTALAKE_DATA_DIR}/python -B ${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/code/extract_fileq.py"
/usr/bin/time -v $P original fileq ${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/out/tiny_original_fileq --pilot
/usr/bin/time -v $P joint_corrected fileq ${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/out/tiny_joint_fileq --pilot --nrows 2000000
/usr/bin/time -v $P joint_corrected sageq ${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/out/tiny_joint_sageq --pilot --nrows 2000000
