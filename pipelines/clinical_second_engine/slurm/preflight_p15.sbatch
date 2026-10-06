#!/bin/bash
#SBATCH --job-name=SE_p15_preflight
#SBATCH --partition=p.hpcl92
#SBATCH --nice=10000
#SBATCH --cpus-per-task=8
#SBATCH --mem=14G
#SBATCH --time=01:30:00
#SBATCH --output=preflight_%j.out
#SBATCH --error=preflight_%j.err
set -uo pipefail
PY=${FASTALAKE_PYTHON:-python3}
S=${FASTALAKE_OUT_DIR}
P=${FASTALAKE_OUT_DIR}/preflight_p15
mk() {  # mk DIR INDEX [tamper]
  mkdir -p $1/inputs; ln -sfn $S/frozen $1/frozen
  $PY -c "
import json,sys
m=json.load(open('$S/inputs/manifest.json')); t=m['tasks'][int(sys.argv[2])]
if len(sys.argv)>3: t['fasta']['sha256']='0'*64
m['tasks']=[t]; json.dump(m,open(sys.argv[1]+'/inputs/manifest.json','w'),indent=1)" $1 $2 ${3:-}
}
# 1) existing task 888 (CAPSCAN FL_public T760_2169, bench job 6064480) reproduces its receipt
mk $P/repro888 888
${PIPELINE_DIR}/task.sh $P/repro888 0 || { echo "REPRO TASK FAILED"; exit 1; }
ID=CAPSCAN__20211017_EXPL0_PeTr_SA_T760_2169__FL_public
ok=1
for f in accepted_5_50.txt accepted_9_50.txt accepted_canonical_peptides.tsv; do
  if cmp -s $S/results/$ID/$f $P/repro888/results/$ID/$f; then echo "REPRO $f identical"; else echo "REPRO $f DIFFERS"; ok=0; fi
done
$PY -c "
import json; a=json.load(open('$S/results/$ID/RECEIPT.json')); b=json.load(open('$P/repro888/results/$ID/RECEIPT.json'))
k=['peptides_q01','target_entries','labels','inputs','frozen_code']; d=[x for x in k if a[x]!=b[x]]
print('RECEIPT fields', 'identical' if not d else 'DIFFER '+str(d), a['peptides_q01']); raise SystemExit(1 if d else 0)" || ok=0
[ $ok = 1 ] || { echo "REPRO FAILED"; exit 1; }
# 2) negative control: same task, recorded FASTA sha256 tampered -> must fail
mk $P/neg_fasta_hash 888 tamper
if ${PIPELINE_DIR}/task.sh $P/neg_fasta_hash 0 > $P/neg_fasta_hash.log 2>&1; then echo "NEGATIVE CONTROL DID NOT FAIL"; exit 1
else echo "NEGATIVE CONTROL FAILED AS REQUIRED: $(grep -m1 -E 'AssertionError' $P/neg_fasta_hash.log)"; fi
# 3) append rule, dry run on the completed pilot participant-15 acquisition; then run that new row end to end
$PY -B ${PIPELINE_DIR}/append_p15.py $S --dry $P/e2e/inputs 20210615_EXPL0_PeTr_SA_T760_1799 || { echo "APPEND DRY FAILED"; exit 1; }
ln -sfn $S/frozen $P/e2e/frozen
${PIPELINE_DIR}/task.sh $P/e2e 1610 || { echo "E2E TASK FAILED"; exit 1; }
cat $P/e2e/results/CAPSCAN__20210615_EXPL0_PeTr_SA_T760_1799__FL_public/RECEIPT.json | $PY -c "import json,sys; r=json.load(sys.stdin); print('E2E RECEIPT', r['status'], r['peptides_q01'])"
# 4) second collector on the e2e run (1 of 332 FL_public receipts present there)
$PY -B ${PIPELINE_DIR}/collect_fl347.py $P/e2e > $P/collect_e2e.log || { echo "COLLECTOR FAILED"; exit 1; }
cat $P/e2e/RESULTS_CAPSCAN_FL347.md
cmp -s $S/inputs/manifest.json $P/../inputs/manifest.json && echo "production manifest untouched (1610 tasks): $($PY -c "import json;print(len(json.load(open('$S/inputs/manifest.json'))['tasks']))")"
echo PREFLIGHT_DONE
