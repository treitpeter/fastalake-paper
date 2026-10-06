#!/bin/bash
# Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_WORK_DIR (working directory).
# Python code addresses inputs as LAB_ROOT/<relative path>, so link it into the working directory.
cd "${FASTALAKE_WORK_DIR:?set FASTALAKE_WORK_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
A=${FASTALAKE_LAB_ROOT}/MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19
B=${FASTALAKE_LAB_ROOT}
RC=$B/MM_review_analysis/recovery_2026-09-15/biological_grouping_recovery/full/cohorts
PY=${FASTALAKE_PYTHON:-python}
for f in "$RC/PREDICT_noclust/direct/feature_assignments.tsv.gz" \
         "$RC/PREDICT_noclust/direct/peptide_to_study_group.tsv" \
         "$RC/PREDICT_noclust/direct/input_manifest.json" \
         "$RC/PREDICT_noclust/reaggregated/representative_annotations.tsv.gz" \
         "$RC/PREDICT_noclust/reaggregated/changed_feature_assignments.tsv.gz" \
         "$A/revised/build_variants.py"; do
  [ -e "$f" ] || { echo "PREFLIGHT FAIL missing: $f"; exit 2; }
done
echo "PREFLIGHT OK repaired PREDICT"
[ "${PREFLIGHT_ONLY:-0}" = 1 ] && { echo "PREFLIGHT ONLY - stopping"; exit 0; }
echo "HOST $(hostname) JOB ${SLURM_JOB_ID:-none}"
export FL_GROUPING_COHORTS="$RC"
"$PY" -B "$A/revised/build_variants.py" PREDICT "$A/inputs/matrices_repaired"
