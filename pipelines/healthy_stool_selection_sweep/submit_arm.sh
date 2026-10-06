#!/bin/bash
# One HSTOOL sweep arm through the REAL chain driver (gate -> evidence -> persample ->
# parsimony -> search -> aggregate), FL_NO_FUNCTIONAL=1, no entrapment sets.
#
#   submit_arm.sh <ml|tp|kmer> <value> <PART_ARRAYS> <PART_BIG> <THROTTLE>
#     ml 6      -> publication_data/arms/peptide_length_sweep_runs/minlen6/cohorts/HSTOOL
#     tp 1.00   -> publication_data/arms/score_sweep_runs/top1.00/cohorts/HSTOOL
#
# Prints, on the last line, "TERMINAL <aggregate job id>" for the queue starter's afterany.
# Why mkdir first: the 2026-09-03 overnight length sweep died at 22:00 because the arm base
# directory did not exist when the submit log was opened (logs/sweep_starter_L6_5850562.err).
set -euo pipefail
KIND=${1:?ml|tp|kmer}; VAL=${2:?value}; PA=${3:?partition for arrays}; PB=${4:?partition for big jobs}; THR=${5:?throttle}
CANON=${FASTALAKE_DATA_DIR}/canonical
CHAIN=$CANON/publication/chain
LIST58=$CANON/publication/analyses/fig1_revision_2026-09-02/cohorts/hstool_sample_list_ruled58.tsv
# The arm variables must come ONLY from the case below. A starter submitted with --export=ALL
# from a shell that had FL_TOP_PERCENT exported ran length 7 at top 0.90 (2026-09-06 17:35).
unset FL_MIN_LENGTH FL_TOP_PERCENT FL_SAMPLE_LIST_OVERRIDE
case "$KIND" in
  ml) BASE=$CANON/publication_data/arms/peptide_length_sweep_runs; ARM=$BASE/minlen$VAL; export FL_MIN_LENGTH=$VAL ;;
  tp) BASE=$CANON/publication_data/arms/score_sweep_runs;          ARM=$BASE/top$VAL;    export FL_TOP_PERCENT=$VAL ;;
  # k-mer arm: same score cut, but the de novo CSVs are the stage-1.5 corrected ones
  # (kmer/src/stage15.sbatch); the chain itself is untouched.
  kmer) BASE=$CANON/publication_data/arms/score_sweep_runs;        ARM=$BASE/top${VAL}_kmer; export FL_TOP_PERCENT=$VAL
        LIST58=$ARM/1.5_kmer/hstool_sample_list_ruled58_kmer_top$VAL.tsv
        [ -s "$LIST58" ] || { echo "no corrected sample list at $LIST58" >&2; exit 1; } ;;
  *) echo "bad kind $KIND" >&2; exit 1 ;;
esac
mkdir -p "$ARM/cohorts/HSTOOL"
export FL_STAGE0_ROOT=$CANON/publication_data/cohorts FL_SKIP_CLUSTER=1 FL_CHAIN_DIR=$CHAIN
export FL_DATA_ROOT=$ARM FL_SAMPLE_LIST_OVERRIDE=$LIST58 FL_NO_FUNCTIONAL=1
export FL_PARTITION=$PA FL_PARTITION_BIG=$PB THROTTLE=$THR
# Memory: k10 HSTOOL MaxRSS was evidence 83G, persample 15G, parsimony 9.5G, search 3.5G,
# aggregate 6.2G (10peptide_length_runs/resource_usage.tsv). These arms admit more than
# k10 (shorter seeds, larger score fraction), so a margin above x1.05 is deliberate.
export FL_MEM_EVIDENCE=100G FL_MEM_PERSAMPLE=32G FL_MEM_PARSIMONY=24G FL_MEM_SEARCH=6G FL_MEM_AGGREGATE=8G
# Parsimony scales ~1.4 GB per M proteins in the sample lake (measured: 8.5 M -> 12.8 GB at
# top 1.00, 17.1 M -> 24.1 GB at length 7). Length 7 OOMed at 24G; shorter seeds admit more.
# Sage: ~0.7 GB per MB of razor FASTA (top 1.00: 5.8 MB -> 4.1 GB). Length 7's razor DB is
# 15 MB (27,851 proteins) and died at 6G within 14 s; shorter seeds give bigger razor sets.
case "$KIND:$VAL" in
  ml:6) export FL_MEM_PARSIMONY=64G FL_MEM_PERSAMPLE=48G FL_MEM_SEARCH=32G ;;  # evidence lake 43.4 M proteins (12:29)
  ml:7) export FL_MEM_PARSIMONY=32G FL_MEM_SEARCH=16G ;;
  ml:8) export FL_MEM_SEARCH=12G ;;
esac
# 2026-09-06 10:5x: parsimony OOM at 12G on the top-1.00/0.80 arms (8.5 M-protein sample lakes,
# 11.2-12.8 GB measured); ~1.4 GB per M proteins. 24G. Persample 32G for the length-6 arms.
unset SBATCH_BEGIN
LOG=$ARM/HSTOOL.submit.log
rc=0
( cd "$CHAIN" && ./submit_cohort.sh HSTOOL "" " " ) > "$LOG" 2>&1 || rc=$?
JF=$ARM/cohorts/HSTOOL/.jobids
if [ $rc -ne 0 ] || ! grep -q "^08_aggregate [0-9]" "$JF" 2>/dev/null; then
    echo "SUBMISSION FAILED for $KIND=$VAL (rc=$rc); cancelling whatever was queued" >&2
    cat "$LOG" >&2
    grep -oE "[0-9]{6,}" "$LOG" | sort -u | xargs -r scancel
    exit 1
fi
grep -E "^cohort=|^# (FL_|mem)" "$LOG" | sed "s/^/  $KIND=$VAL /"
echo "TERMINAL $(awk '$1=="08_aggregate" {print $2}' "$JF")"
