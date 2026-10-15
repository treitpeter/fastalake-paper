#!/usr/bin/env python3
"""Accepted peptides of one MG search with the frozen rule (split_half_repair pipeline.accepted() without the
odd-scan assertion; identical to iteration3/disjoint/worker.py count_any): rank 1, target, peptide_q <= 0.01,
modifications removed, I->L, 5-50 and nested 9-50. Writes ACCEPTED_5_50.txt, ACCEPTED_9_50.txt, COUNTS.json.
For the control sample, compares with the frozen MG sets of Review_v16/Confidence (expected_sets MG) and exits 1
unless the 5-50 and 9-50 counts equal the archived counts."""
import csv, json, re, sys, hashlib
from pathlib import Path
out, sample = Path(sys.argv[1]), sys.argv[2]
CONF_MAN = Path('LAB_ROOT/review/v15_integration/ready_integration_2026-09-28/source_additions/Review_v16/Confidence/lineage/manifest.json')
acc = set(); n = 0; nd = 0
with open(out / 'results.sage.tsv') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        n += 1; nd += int(r['label']) == -1
        if int(r['rank']) == 1 and int(r['label']) == 1 and float(r['peptide_q']) <= .01:
            p = re.sub(r'\[[^\]]*\]', '', r['peptide']).replace('I', 'L')
            assert p.isalpha() and p.isupper(), p
            if 5 <= len(p) <= 50:
                acc.add(p)
assert n > 0 and nd > 0
sets = {'5': sorted(acc), '9': sorted(p for p in acc if len(p) >= 9)}
res = dict(sample=sample, psms=n, decoy_psms=nd, peptides_5_50=len(sets['5']), peptides_9_50=len(sets['9']))
for L, s in sets.items():
    p = out / f'ACCEPTED_{L}_50.txt'; p.write_text('\n'.join(s) + '\n'); res[f'accepted_{L}_sha256'] = hashlib.sha256(p.read_bytes()).hexdigest()
man = json.loads(CONF_MAN.read_text())
o = next((x for x in man['ordinary'] if x['cohort'] == 'CAPSCAN' and x['sample'] == sample), None)
ok = True
if o is not None:   # an existing 282 pair: the control
    for L in ['5', '9']:
        exp = set(Path(o['expected_sets']['MG'][L]['path']).read_text().split())
        mine = set(sets[L])
        res[f'archived_{L}_50'] = len(exp); res[f'count_equal_{L}'] = len(exp) == len(mine); res[f'set_equal_{L}'] = exp == mine
        res[f'only_rerun_{L}'] = len(mine - exp); res[f'only_archived_{L}'] = len(exp - mine)
        ok &= len(exp) == len(mine)
    res['control_status'] = 'PASS' if ok else 'FAIL'
json.dump(res, open(out / 'COUNTS.json', 'w'), indent=1); print(json.dumps(res))
sys.exit(0 if ok else 1)
