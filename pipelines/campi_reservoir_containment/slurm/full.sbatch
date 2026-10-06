#!/bin/bash
#SBATCH --job-name=R8B_full
#SBATCH --partition=p.hpcl8,p.hpcl92,p.hpcl93
#SBATCH --nice=10000
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=155G
#SBATCH --time=12:00:00
#SBATCH --output=logs/full_%j.out
#SBATCH --error=logs/full_%j.err
# Full run (dependency afterok on the pilot; --begin 22:00 for p.hpcl8/p.hpcl93 daytime courtesy):
# (1) one streaming pass over the whole 136 GB reservoir: exact-protein presence + tryptic-context
#     peptide presence for all eight proteomes and their shuffled controls, sha256 verified against
#     the manifest recorded by the pilot; (2) one streaming pass over the 81 GB joint-headers TSV for
#     gut / UniProt / both attribution; (3) B. longum peptides through the production extractor on the
#     whole reservoir (substring presence, the production lookup); (4) summary.
# --mem: production extractor on this reservoir measured MaxRSS 145.6 GB (job 5895250; the file is
# memory-mapped) x 1.05 = 153 GB -> 155G. The Python passes need < 3 GB (pilot measures it).
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=16
D=${FASTALAKE_OUT_DIR}
PY=${FASTALAKE_PYTHON:-python3}
MAN=$D/inputs/manifest.json
OUT=$D/results/full
mkdir -p "$OUT/extractor"
j() { "$PY" -c "import json;print(json.load(open('$MAN'))['inputs']['$1']['path'])"; }
LAKE=$(j reservoir_fasta); JH=$(j joint_headers_tsv); REF=$(j reference_fasta); MEM=$(j members_tsv); EXT=$(j extractor_binary)
T0=$(date +%s)
step() { echo "=== $(date -u +%H:%M:%S) [$(( $(date +%s) - T0 ))s] $*"; }

step "preflight"
bash "${PIPELINE_DIR}/preflight.sh"
"$PY" -c "
import json,sys; m=json.load(open('$MAN'))['inputs']
for k in ('reservoir_fasta','joint_headers_tsv'):
    assert m[k]['sha256'], f'{k}: sha256 not recorded by the pilot; refusing to run'
print('manifest: large-input sha256 recorded by pilot job', m['reservoir_fasta']['sha256_job_id'])"
[ -f "$D/results/pilot/COMPLETE.json" ] || { echo "pilot COMPLETE.json missing" >&2; exit 1; }

step "1 protein + tryptic-context pass, whole reservoir"
"$PY" -B "${PIPELINE_DIR}/overlap.py" protein --manifest "$MAN" --reference "$REF" --members "$MEM" --reservoir "$LAKE" \
  --out "$OUT" --expect-entries 375333654

step "2 source pass, whole joint-headers TSV"
"$PY" -B "${PIPELINE_DIR}/overlap.py" sources --manifest "$MAN" --joint-headers "$JH" --out "$OUT"

step "3 digest B. longum (real + shuffled)"
"$PY" -B "${PIPELINE_DIR}/overlap.py" digest --reference "$REF" --members "$MEM" --out "$OUT" --label "Bifidobacterium longum"

step "4 production extractor on the whole reservoir"
TIME=""; [ -x /usr/bin/time ] && TIME="/usr/bin/time -v"
$TIME "$EXT" -p "$OUT/extractor_queries.txt" -d "$LAKE" -o "$OUT/extractor/matched.fasta" -s DB \
  --min-length 9 --max-length 50 --normalize-il true --output-hits "$OUT/extractor/hits.tsv" -t 16 \
  > "$OUT/extractor/extractor.log" 2> "$OUT/extractor/extractor.time"
tail -5 "$OUT/extractor/extractor.log"; grep -E "Maximum resident|Elapsed" "$OUT/extractor/extractor.time" || true

step "5 extractor summary"
"$PY" -B "${PIPELINE_DIR}/overlap.py" extractor-summary --out "$OUT" --hits "$OUT/extractor/hits.tsv" --label "Bifidobacterium longum"
ls -la "$OUT/extractor/" | tee -a "$OUT/extractor/sizes.txt"
rm -f "$OUT/extractor/matched.fasta"

step "done"
"$PY" - "$OUT" "$T0" <<'EOF'
import json, os, sys, time
out, t0 = sys.argv[1], int(sys.argv[2])
json.dump({"stage": "full", "job_id": os.environ.get("SLURM_JOB_ID"), "elapsed_s": int(time.time()) - t0,
           "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()), "host": os.uname().nodename,
           "protein_pass": json.load(open(f"{out}/protein_pass.json")),
           "sources_pass": json.load(open(f"{out}/sources_pass.json")),
           "extractor_summary": json.load(open(f"{out}/extractor_summary.json"))},
          open(f"{out}/COMPLETE.json", "w"), indent=2, sort_keys=True)
print("FULL COMPLETE")
EOF
