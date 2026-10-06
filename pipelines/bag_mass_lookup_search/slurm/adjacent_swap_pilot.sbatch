#!/bin/bash
# Adjacent-swap lookup (originally R8 task C2). Required: FASTALAKE_OUT_DIR (run directory with inputs/, bin/, C2/),
# FASTALAKE_DATA_DIR (lab tree holding the HSTOOL reservoir and A02 inputs), PIPELINE_DIR (bag_mass_lookup_search), FASTALAKE_PYTHON.
#SBATCH -J R8C_C2_swap_pilot
#SBATCH --partition=p.hpcl8,p.hpcl92,p.hpcl93
#SBATCH --nice=10000
#SBATCH --nodes=1
#SBATCH --cpus-per-task=24
#SBATCH --mem=73G
#SBATCH --time=01:30:00
# --mem: A02 job 5920358 (same extractor, same 61.3 GB reservoir, 207k pooled patterns) MaxRSS 69,335,758 K = 69.3 GB,
# of which 61.3 GB is the memory-mapped reservoir counted in RSS. 69.3 x 1.05 = 72.8 -> 73G.
set -euo pipefail
R="${FASTALAKE_OUT_DIR:?set FASTALAKE_OUT_DIR}"; SCR="${PIPELINE_DIR:?set PIPELINE_DIR}/adjacent_swap"
PY="${FASTALAKE_PYTHON:-python}"
EXT=$R/bin/fasta_extractor; PARS=$R/bin/parsimony_engine
RES=${FASTALAKE_DATA_DIR}/canonical/publication_data/cohorts/HSTOOL/0_lake/HSTOOL_lake.fasta
Q=$R/C2/queries; O=$R/C2/pilot
echo "job $SLURM_JOB_ID node $(hostname) start $(date -u +%FT%TZ)"
echo "== preflight"; "$PY" "$SCR/preflight.py"
echo "== reservoir sha256 (61.3 GB; recorded once, re-checked by the full job)"
if [ ! -s "$R/inputs/reservoir.sha256.json" ]; then
    H=$(sha256sum "$RES" | cut -d' ' -f1)
    printf '{"path": "%s", "bytes": %s, "sha256": "%s", "computed_by_job": "%s"}\n' "$RES" "$(stat -c %s "$RES")" "$H" "$SLURM_JOB_ID" > "$R/inputs/reservoir.sha256.json"
fi
cat "$R/inputs/reservoir.sha256.json"
echo "== build query sets (pilot 1,000 + full; deterministic seed 20260917)"
"$PY" "$SCR/build_queries.py"
run() {  # run <set name> ; extractor over one CSV, production flags, no --top-percent (queries are the already-retained top-30% set)
    local s=$1; mkdir -p "$O/$s"
    echo "-- extractor $s : $(($(wc -l < "$Q/$s/queries.csv") - 1)) query rows"
    /usr/bin/time -v "$EXT" --peptides-dir "$Q/$s" --database "$RES" --output "$O/$s/matches.fasta" --output-hits "$O/$s/hits.tsv" \
        --min-length 9 --max-length 50 --normalize-il true -t 24 > "$O/$s/extractor.out" 2> "$O/$s/extractor.err"
    grep -E "Maximum resident|Elapsed \(wall" "$O/$s/extractor.err" | sed 's/^/   /'
    ls "$O/$s"; cat "$O/$s/"*.stats.json 2>/dev/null | "$PY" -c "import sys,json; d=json.load(sys.stdin); print('   proteins_matched', d['proteins_matched'], 'total_peptides', d['total_peptides'], 'seconds', round(d['total_seconds']))"
}
for s in P_REAL_orig P_REAL_all P_SIM_orig P_SIM_all; do run $s; done
echo "== parsimony (razor, hash-acc) on the pilot run B candidates, exercising step (c)"
mkdir -p "$O/parsimony_P_REAL_all" "$O/parsimony_P_REAL_orig"
"$PARS" --strategy razor --tiebreak hash-acc --fasta "$O/P_REAL_all/matches.fasta" --peptides "$Q/P_REAL_all/queries.csv" \
    --output-dir "$O/parsimony_P_REAL_all" --sample P_REAL_all --min-length 9 --max-length 50 --normalize-il -t 8 > "$O/parsimony_P_REAL_all/parsimony.out" 2>&1
"$PARS" --strategy razor --tiebreak hash-acc --fasta "$O/P_REAL_orig/matches.fasta" --peptides "$Q/P_REAL_orig/queries.csv" \
    --output-dir "$O/parsimony_P_REAL_orig" --sample P_REAL_orig --min-length 9 --max-length 50 --normalize-il -t 8 > "$O/parsimony_P_REAL_orig/parsimony.out" 2>&1
ls "$O/parsimony_P_REAL_all" "$O/parsimony_P_REAL_orig"
echo "== summary: authentic pilot"
"$PY" "$SCR/summarize.py" "$Q/P_REAL_all" "$O/P_REAL_orig" "$O/P_REAL_all" "$(ls "$O"/parsimony_P_REAL_all/*.fasta | head -1)" "$(ls "$O"/parsimony_P_REAL_orig/*.fasta | head -1)" "$O/SUMMARY_P_REAL.json"
echo "== summary: scrambled negative control (seed 20260913), same expansion"
"$PY" "$SCR/summarize.py" "$Q/P_SIM_all" "$O/P_SIM_orig" "$O/P_SIM_all" - - "$O/SUMMARY_P_SIM.json"
echo "PILOT_OK end $(date -u +%FT%TZ)"
