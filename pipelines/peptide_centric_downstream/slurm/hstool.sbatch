#!/bin/bash
# Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_WORK_DIR (working directory).
# Python code addresses inputs as LAB_ROOT/<relative path>, so link it into the working directory.
cd "${FASTALAKE_WORK_DIR:?set FASTALAKE_WORK_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
A=${FASTALAKE_LAB_ROOT}/MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19
B=${FASTALAKE_LAB_ROOT}
PY=${FASTALAKE_PYTHON:-python}
MODE="$1"
for f in "$A/revised/hstool_pc.py" "$A/revised/derive_hstool_pc.py" "$A/revised/common_pepcentric.py" \
         "$B/release_prep/integration_2026-09-14/repo/manuscript/current/code/derive_hstool.py" \
         "$B/release_prep/integration_2026-09-14/repo/manuscript/current/code/common.py" \
         "$B/release_prep/integration_2026-09-14/repo/manuscript/current/source_data/Figure3/pair_summary_local_candidate.tsv" \
         "$A/inputs/matrices/HSTOOL_KO_local_candidate_all.tsv.gz" \
         "$A/inputs/matrices/HSTOOL_KO_species.npz" \
         "$A/inputs/matrices/HSTOOL_KO_species_original.npz"; do
  [ -e "$f" ] || { echo "PREFLIGHT FAIL missing: $f"; exit 2; }
done
case "$MODE" in baseline|simplified|original) ;; *) echo "PREFLIGHT FAIL bad mode: $MODE"; exit 2;; esac
echo "PREFLIGHT OK hstool $MODE"
[ "${PREFLIGHT_ONLY:-0}" = 1 ] && { echo "PREFLIGHT ONLY - stopping"; exit 0; }
echo "HOST $(hostname) JOB ${SLURM_JOB_ID:-none} MODE $MODE"
"$PY" -B "$A/revised/hstool_pc.py" "$MODE"
