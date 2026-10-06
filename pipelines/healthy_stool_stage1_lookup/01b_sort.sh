#!/bin/bash -l
#SBATCH --job-name=v27_sort
#SBATCH --output=logs/01b_sort_%x_%j.out
#SBATCH --error=logs/01b_sort_%x_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=96G
# NO --partition DIRECTIVE: directives are parsed before the shell runs, so they
# cannot read ${FL_PARTITION}. It is passed at submit time -- see submit_cohort.sh
# and docs/SITE_CONFIG.md. A literal here would hard-code one cluster.
#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
# ============================================================================
# Stage 1.5 — sort the clustering input by accession. NEW IN v26.
#
# This stage did not exist in the published chain, and its absence is defect D23.
#
# lake_builder writes records in the order it reads sources, so record order
# encodes which catalogue was named first. Stage 2 runs mmseqs with --shuffle 0,
# which makes MMseqs2 respect input order, and greedy set-cover resolves tied
# candidate representatives by the lower database key. So whichever source was
# listed first wins representatives systematically.
#
# That is not cosmetic. Identity and coverage are computed AGAINST the
# representative, so changing it changes which sequences join the cluster, hence
# the cluster count and the size of the searchable database -- measured at ~18%.
# And the representative's accession is what db_breakdown credits for the whole
# cluster.
#
# Accessions are COHORT_<sha256hex>, so sorting by accession IS sorting by sequence
# content: a total order no source ordering can influence.
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

IN=${ROOT}/1_evidence/evidence_lake.fasta
OUT=${ROOT}/1_evidence/evidence_lake.sorted.fasta

[ -s "${IN}" ] || fail "evidence lake missing: ${IN}"
[ -x "${FASTA_SORT}" ] || fail "fasta_sort_by_id missing: ${FASTA_SORT}"

N_IN=$(count_headers "${IN}")
log "sorting ${N_IN} records by accession"
# Provenance guard. Refuses to run if this binary's crate has uncommitted
# changes or if the binary is older than its source, and records commit +
# binary sha256 beside the output. See assert_code_committed in _config.sh.
STAGE_NAME=03_sort assert_code_committed "${ROOT}/1_evidence" "${FASTA_SORT}"

"${FASTA_SORT}" --input "${IN}" --output "${OUT}" || fail "fasta_sort_by_id failed"

N_OUT=$(count_headers "${OUT}")
log "sorted: ${N_OUT} records"
# Sorting must not lose or invent a sequence.
[ "${N_OUT}" -eq "${N_IN}" ] || fail "record count changed in sort: ${N_IN} -> ${N_OUT}"

# Verify the ordering here as well as in stage 04. Two independent checks, because
# an unsorted input produces a perfectly plausible clustering that silently favours
# one catalogue -- there is no downstream symptom to notice.
log "verifying ascending unique accession order"
BAD=$(grep '^>' "${OUT}" | sed 's/^>//' | cut -d' ' -f1 | LC_ALL=C sort -c 2>&1 || true)
[ -z "${BAD}" ] || fail "output is NOT sorted: ${BAD}"
N_UNIQ=$(grep '^>' "${OUT}" | sed 's/^>//' | cut -d' ' -f1 | LC_ALL=C sort -u | wc -l)
[ "${N_UNIQ}" -eq "${N_OUT}" ] || fail "accessions not unique: ${N_UNIQ} unique of ${N_OUT}"
log "PASS sorted, unique, ${N_OUT} records"
