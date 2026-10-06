#!/bin/bash -l
#SBATCH --job-name=v27_evidence
#SBATCH --output=logs/01_evidence_%x_%j.out
#SBATCH --error=logs/01_evidence_%x_%j.err
#SBATCH --time=24:00:00
#SBATCH --mem=500G
# NO --partition DIRECTIVE: directives are parsed before the shell runs, so they
# cannot read ${FL_PARTITION}. It is passed at submit time -- see submit_cohort.sh
# and docs/SITE_CONFIG.md. A literal here would hard-code one cluster.
#SBATCH --nodes=1
#SBATCH --cpus-per-task=32
# ============================================================================
# Stage 1 — evidence filter: keep only lake sequences containing an observed
# de novo peptide, pooled across the whole cohort.
#
# Parameters are the canonical ones and are NOT tuned here: --top-percent ${TOP_PERCENT} (default 0.30),
# --min-length ${MIN_LENGTH} (default 9), --max-length 50, --normalize-il true.
# Changing any of them changes the searchable database, so they stay pinned to the
# published recipe. The one exception is the minimum length, which is a declared arm
# variable (FL_MIN_LENGTH, _config.sh) so that stages 1, 3 and 4 always agree on it.
# ============================================================================
set -uo pipefail
# ONE COPY OF THE CONFIG, AND IT IS THE TRACKED ONE.
#
# This line used to be an absolute path into ${FL_CHAIN_DIR:-${REPO}/canonical/publication/chain}/ -- a SECOND, untracked copy
# of the chain. That is the defect the 2026-08-03 archive names as the root cause of the
# v25/v26 churn: edits went to one copy and were hand-copied to the other, and eleven of
# fifteen scripts had drifted between them.
#
# It is also how renaming 00_paths.sh -> _config.sh broke every stage: the old tree still
# holds 00_paths.sh, so this absolute path resolved to nothing and the chain could not start.
#
# The default now points INTO CANONICAL. Override only to run a different checkout:
#   sbatch --export=ALL,FL_CHAIN_DIR=/path/to/checkout/publication/chain ...
# 🔴 UNDER sbatch, ${BASH_SOURCE[0]} IS NOT IN THE CHECKOUT.
# SLURM copies the batch script to /var/spool/slurmd/job<N>/slurm_script and runs THAT,
# so the dirname below resolved to the spool directory and the stage died in five seconds
# with "FATAL: no _config.sh at /var/spool/slurmd/job5840734" (found 2026-09-02, stage 1
# of the CAMPI augmented-lake arm). It only ever worked when the submitting shell had
# already exported FL_CHAIN_DIR; a bare `sbatch 01_evidence.sh` from the chain directory
# itself could not start. Fall back to the submit directory, and accept a candidate ONLY
# if _config.sh is readable there -- a file no stage creates, so the check cannot
# manufacture its own success the way a `mkdir -p` test can.
if [ -z "${FL_CHAIN_DIR:-}" ]; then
    for _fl_d in "$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)" "${SLURM_SUBMIT_DIR:-}"; do
        [ -n "${_fl_d}" ] && [ -r "${_fl_d}/_config.sh" ] && { FL_CHAIN_DIR="${_fl_d}"; break; }
    done
    unset _fl_d
fi
FL_CHAIN_DIR=${FL_CHAIN_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" 2>/dev/null && pwd)}
[ -r "${FL_CHAIN_DIR}/_config.sh" ] || { echo "FATAL: no _config.sh at ${FL_CHAIN_DIR}" >&2; exit 1; }
source "${FL_CHAIN_DIR}/_config.sh"

OUT_DIR=${ROOT}/1_evidence
POOL_DIR=${OUT_DIR}/_pooled_csvs
rm -rf "${POOL_DIR}"; mkdir -p "${POOL_DIR}"

[ -s "${LAKE}" ] || fail "lake missing: ${LAKE}"
[ -x "${FASTA_EXTRACTOR}" ] || fail "fasta_extractor missing: ${FASTA_EXTRACTOR}"

N_EXPECTED=$(n_samples)
N_POOLED=0
for i in $(seq 0 $(( N_EXPECTED - 1 ))); do
    ROW=$(sample_row "$i")
    SAMPLE=$(echo "${ROW}" | cut -f1)
    CSV=$(echo "${ROW}" | cut -f2)
    [ -n "${SAMPLE}" ] && [ -e "${CSV}" ] || fail "bad row ${i}: ${ROW}"
    ln -sf "${CSV}" "${POOL_DIR}/${SAMPLE}_predictions.csv"
    N_POOLED=$((N_POOLED + 1))
done
log "pooled ${N_POOLED}/${N_EXPECTED} prediction CSVs"
# Exact, not a floor. The canonical script used `-ge 300`, which would have
# accepted 300 of 306 and quietly changed the cohort.
[ "${N_POOLED}" -eq "${N_EXPECTED}" ] || fail "pooled ${N_POOLED} != expected ${N_EXPECTED}"

log "stage 1 evidence — lake=${LAKE} ($(count_headers "${LAKE}") proteins)"
# Provenance guard. Refuses to run if this binary's crate has uncommitted
# changes or if the binary is older than its source, and records commit +
# binary sha256 beside the output. See assert_code_committed in _config.sh.
STAGE_NAME=02_evidence assert_code_committed "$(dirname "${EV_OUT:-${ROOT}/1_evidence}")" "${FASTA_EXTRACTOR}"

"${FASTA_EXTRACTOR}" \
    --peptides-dir "${POOL_DIR}" \
    --database     "${LAKE}" \
    --output       "${OUT_DIR}/evidence_lake.fasta" \
    --output-hits  "${OUT_DIR}/evidence_hits.tsv" \
    --top-percent "${TOP_PERCENT}" --min-length "${MIN_LENGTH}" --max-length 50 \
    --normalize-il true \
    -t 32 \
    || fail "fasta_extractor failed"

N=$(count_headers "${OUT_DIR}/evidence_lake.fasta")
[ "${N}" -gt 0 ] || fail "evidence lake is empty"
log "evidence lake: ${N} proteins"
echo "${N}" > "${OUT_DIR}/.n_evidence"
