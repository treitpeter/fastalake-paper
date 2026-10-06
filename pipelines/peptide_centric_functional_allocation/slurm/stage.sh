#!/bin/bash
# Copy of A33 revised/run.sh (producer/run.sh). Changes, roster only:
#  - grouping run chosen by MODE (old = Sept-10 manuscript run, 331 roster; 347 = C2 Sept-10-driver run on fixed_roster_347)
#    and passed to the unmodified producer via its own FL_GROUPING_COHORTS variable
#  - producer = byte-identical copy producer/build_variants.py (sha256 asserted); output folder per MODE
#  - added gate: acquisition set of the grouping run must equal the roster; output must not already be complete
# Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_WORK_DIR (working directory).
# Python code addresses inputs as LAB_ROOT/<relative path>, so link it into the working directory.
cd "${FASTALAKE_WORK_DIR:?set FASTALAKE_WORK_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
set -uo pipefail
MODE="$1"; COH=CAPSCAN; PY=$HOME/miniconda3/bin/python
S=${FASTALAKE_LAB_ROOT}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/suppfigs_n15/SF12/alloc347
B=${FASTALAKE_LAB_ROOT}
C2=$B/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15
case "$MODE" in
  old) export FL_GROUPING_COHORTS=${FL_GROUPING_COHORTS:-$B/simplified_grouping_results/manuscript_cohorts_2026-09-10_2200/cohorts}; ROSTER=${ROSTER:-$S/roster/roster_331.txt}; OUT=${OUT:-$S/repro331};;
  347) export FL_GROUPING_COHORTS=${FL_GROUPING_COHORTS:-$C2/fig5/grouping/cohorts}; ROSTER=${ROSTER:-$S/roster/roster_347.txt}; OUT=${OUT:-$S/run347};;
  *) echo "usage: stage.sh old|347"; exit 2;;
esac
G=$FL_GROUPING_COHORTS/${COH}_noclust
fail(){ echo "PREFLIGHT FAIL $1"; exit 2; }
for f in "$G/direct/feature_assignments.tsv.gz" "$G/direct/peptide_to_study_group.tsv" \
         "$G/direct/input_manifest.json" "$G/reaggregated/representative_annotations.tsv.gz" \
         "$G/reaggregated/changed_feature_assignments.tsv.gz" "$G/COMPLETE.json" \
         "$B/canonical/publication_data/arms/noclust_2026-08-07/cohorts/$COH/11_eggnog/${COH}_eggnog_accession_keyed.tsv" \
         "$PY" "$S/producer/build_variants.py" "$ROSTER"; do
  [ -e "$f" ] || fail "missing: $f"
done
[ "$(sha256sum < $S/producer/build_variants.py | cut -d' ' -f1)" = 57ea24e6d2376f029bee295054a809fccb139d696daab51d1dda5ca2f25e92a9 ] || fail "producer drifted"
"$PY" -c "import ahocorasick,numpy" 2>/dev/null || fail "python deps"
cmp -s <(cut -f1 "$G/direct/acquisition_summary.tsv" | tail -n +2 | sort) <(sort "$ROSTER") || fail "grouping acquisitions != roster $ROSTER"
[ -e "$OUT/${COH}_variant_stats.json" ] && fail "output already complete: $OUT"
echo "PREFLIGHT OK $COH $MODE roster $(wc -l < $ROSTER) grouping $G"
[ "${PREFLIGHT_ONLY:-0}" = 1 ] && { echo "PREFLIGHT ONLY - stopping"; exit 0; }
echo "HOST $(hostname) JOB ${SLURM_JOB_ID:-none} COHORT $COH MODE $MODE"
"$PY" -B "$S/producer/build_variants.py" "$COH" "$OUT"
