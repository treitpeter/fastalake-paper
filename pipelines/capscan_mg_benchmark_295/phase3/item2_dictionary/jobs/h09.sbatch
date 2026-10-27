#!/bin/bash
# args: <run dir> <n>; frozen h09_evidence_ablation.py (sha256 161aefe8) + production apt_env interpreter, then SFig08 b-d table
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=LAB_ROOT/reviews/quantification_mechanisms_2026-09-13/apt_env/bin/python
$PY -B LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item2_dictionary/code/make_h09_manifest.py "$1" "$2"
$PY -B LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item2_dictionary/code/h09_evidence_ablation.py "$1/h09" full 0
mkdir -p "$1/values"
$PY -B LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item2_dictionary/code/ablation_values.py "$1/h09/full_CAPSCAN/acquisition_summary.tsv" "$2" LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item2_dictionary/../../phase2/repo_base/source_data/SFig08/clinical_fixed_dictionary_ablation.tsv "$1/values"
