#!/bin/bash
#SBATCH --job-name=FL_D3_compare
#SBATCH --partition=p.hpcl8,p.hpcl93,p.hpcl92
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=01:00:00
#SBATCH --output=${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/slurm/%x_%j.log
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
/usr/bin/time -v ${FASTALAKE_DATA_DIR}/python -B ${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/code/compare_fileq.py ${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/out/compare ${FASTALAKE_DATA_DIR}/MM_review_analysis/A21_7_hstool_quant_comparison/results/original ${FASTALAKE_DATA_DIR}/MM_review_analysis/A21_7_hstool_quant_comparison/results/union_individual ${FASTALAKE_DATA_DIR}/MM_review_analysis/A21_7_hstool_quant_comparison/results/joint_corrected ${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/out/joint_sageq_repro ${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/out/joint_fileq ${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/out/original_fileq
