#!/bin/bash
# Adjacent-swap lookup (originally R8 task C2). Required: FASTALAKE_OUT_DIR (run directory with inputs/, bin/, C2/),
# FASTALAKE_DATA_DIR (lab tree holding the HSTOOL reservoir and A02 inputs), PIPELINE_DIR (bag_mass_lookup_search), FASTALAKE_PYTHON.
#SBATCH -J R8C_C2_swap_full
#SBATCH --partition=p.hpcl8,p.hpcl92,p.hpcl93
#SBATCH --nice=10000
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=80G
#SBATCH --time=03:00:00
# --mem: benchmark A02 5920358 = 69.3 GB at 207k patterns (61.3 GB of it the mapped reservoir). The largest set here,
# SIM_all, has 454k patterns; the automaton is a contiguous NFA (small) and the query-dependent part is hit-driven,
# so the extrapolated peak is ~72-75 GB; 80G covers it. Parsimony on a ~3 GB matches.fasta needs ~10 GB (chain 04 sizing).
set -euo pipefail
R="${FASTALAKE_OUT_DIR:?set FASTALAKE_OUT_DIR}"; SCR="${PIPELINE_DIR:?set PIPELINE_DIR}/adjacent_swap"
PY="${FASTALAKE_PYTHON:-python}"
EXT=$R/bin/fasta_extractor; PARS=$R/bin/parsimony_engine
RES=${FASTALAKE_DATA_DIR}/canonical/publication_data/cohorts/HSTOOL/0_lake/HSTOOL_lake.fasta
Q=$R/C2/queries; O=$R/C2/full
echo "job $SLURM_JOB_ID node $(hostname) start $(date -u +%FT%TZ)  pilot job ${PILOT_JOB:-?}"
echo "== pilot accounting (measured, for the record)"; sacct -j "${PILOT_JOB:-0}" -P --format=JobIDRaw,JobName,State,Elapsed,MaxRSS || true
echo "== preflight"; "$PY" "$SCR/preflight.py"
[ -s "$R/inputs/reservoir.sha256.json" ] || { echo "FATAL: reservoir sha256 not recorded by the pilot"; exit 1; }
[ -s "$Q/QUERIES.json" ] && [ -s "$Q/REAL_all/queries.csv" ] && [ -s "$Q/SIM_all/queries.csv" ] || { echo "FATAL: query sets missing"; exit 1; }
[ -s "$R/C2/pilot/SUMMARY_P_REAL.json" ] || { echo "FATAL: pilot summary missing"; exit 1; }
run() {
    local s=$1; mkdir -p "$O/$s"
    echo "-- extractor $s : $(($(wc -l < "$Q/$s/queries.csv") - 1)) query rows"
    /usr/bin/time -v "$EXT" --peptides-dir "$Q/$s" --database "$RES" --output "$O/$s/matches.fasta" --output-hits "$O/$s/hits.tsv" \
        --min-length 9 --max-length 50 --normalize-il true -t 24 > "$O/$s/extractor.out" 2> "$O/$s/extractor.err"
    grep -E "Maximum resident|Elapsed \(wall" "$O/$s/extractor.err" | sed 's/^/   /'
    cat "$O/$s/"*.stats.json 2>/dev/null | "$PY" -c "import sys,json; d=json.load(sys.stdin); print('   proteins_matched', d['proteins_matched'], 'total_peptides', d['total_peptides'], 'seconds', round(d['total_seconds']))"
}
for s in REAL_orig REAL_all SIM_orig SIM_all; do run $s; done
echo "== parsimony (razor, hash-acc): run A and run B of the authentic arm"
for s in REAL_orig REAL_all; do
    mkdir -p "$O/parsimony_$s"
    /usr/bin/time -v "$PARS" --strategy razor --tiebreak hash-acc --fasta "$O/$s/matches.fasta" --peptides "$Q/$s/queries.csv" \
        --output-dir "$O/parsimony_$s" --sample "$s" --min-length 9 --max-length 50 --normalize-il -t 8 > "$O/parsimony_$s/parsimony.out" 2> "$O/parsimony_$s/parsimony.err"
    grep -E "Maximum resident|Elapsed \(wall" "$O/parsimony_$s/parsimony.err" | sed 's/^/   /'
done
A02REF=$("$PY" -c "import json; print(json.load(open('${FASTALAKE_DATA_DIR}/MM_review_analysis/A02_4_2_lookup_specificity/inputs/REAL_sample_lake.stats.json'))['proteins_matched'])")
echo "== summary: authentic (positive control: run A proteins must equal A02 archived $A02REF)"
"$PY" "$SCR/summarize.py" "$Q/REAL_all" "$O/REAL_orig" "$O/REAL_all" "$(ls "$O"/parsimony_REAL_all/*.fasta | head -1)" "$(ls "$O"/parsimony_REAL_orig/*.fasta | head -1)" "$O/SUMMARY_REAL.json" "$A02REF"
echo "== summary: scrambled negative control (seed 20260913)"
"$PY" "$SCR/summarize.py" "$Q/SIM_all" "$O/SIM_orig" "$O/SIM_all" - - "$O/SUMMARY_SIM.json"
echo "== output table"
"$PY" - "$O" <<'PYEOF'
import json, sys
O = sys.argv[1]
rows = [("arm", "originals", "swap_only_variants", "(a) proteins_by_originals", "(b) additional_by_swaps_only", "b/a", "originals_without_exact_hit", "rescued_by_swap", "(c) swap_only_in_razor_cover", "razor_cover_size_B", "razor_cover_size_A")]
for arm in ("REAL", "SIM"):
    S = json.load(open(f"{O}/SUMMARY_{arm}.json")); c = S.get("c_parsimony_runB", {}); cA = S.get("c_parsimony_runA", {})
    rows.append((arm, S["n_originals"], S["n_swap_only_variants"], S["a_proteins_by_originals"], S["b_additional_proteins_by_swaps_only"], S["b_over_a_fraction"],
                 S["originals_without_exact_hit"], S["originals_without_exact_hit_rescued_by_a_swap"], c.get("swap_only_candidates_in_cover", "-"), c.get("proteins_in_razor_cover", "-"), cA.get("proteins_in_razor_cover", "-")))
with open(f"{O}/RESULT_TABLE.tsv", "w") as f:
    for r in rows: f.write("\t".join(map(str, r)) + "\n")
print(open(f"{O}/RESULT_TABLE.tsv").read())
PYEOF
echo "FULL_OK end $(date -u +%FT%TZ)"
