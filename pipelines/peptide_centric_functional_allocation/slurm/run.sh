#!/bin/bash
# Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_WORK_DIR (working directory).
# Python code addresses inputs as LAB_ROOT/<relative path>, so link it into the working directory.
cd "${FASTALAKE_WORK_DIR:?set FASTALAKE_WORK_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
set -uo pipefail
SH="$1"; COH="$2"; PY=$HOME/miniconda3/bin/python
B=${FASTALAKE_LAB_ROOT}
G=$B/simplified_grouping_results/manuscript_cohorts_2026-09-10_2200/cohorts/${COH}_noclust
fail(){ echo "PREFLIGHT FAIL $1"; exit 2; }
for f in "$G/direct/feature_assignments.tsv.gz" "$G/direct/peptide_to_study_group.tsv" \
         "$G/direct/input_manifest.json" "$G/reaggregated/representative_annotations.tsv.gz" \
         "$B/canonical/publication_data/arms/noclust_2026-08-07/cohorts/$COH/11_eggnog/${COH}_eggnog_accession_keyed.tsv" \
         "$PY" "$SH/build2/build_variants.py"; do
  [ -e "$f" ] || fail "missing: $f"
done
"$PY" -c "import ahocorasick,numpy" 2>/dev/null || fail "python deps"
echo "PREFLIGHT OK $COH"
[ "${PREFLIGHT_ONLY:-0}" = 1 ] && { echo "PREFLIGHT ONLY - stopping"; exit 0; }
echo "HOST $(hostname) JOB ${SLURM_JOB_ID:-none} COHORT $COH"
"$PY" -B "$SH/build2/build_variants.py" "$COH" "$SH/build2/out"
