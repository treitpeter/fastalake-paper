#!/usr/bin/env python3
"""Second-engine CapScan paired ratios at 295: production second_engine_all/results_paired.tsv rows (282 planned)
+ the 13 new pairs (MSFragger FL_public / FL_cohortMG receipts from second_engine_all, MSFragger MG receipts from
mg295/second_engine; Sage ratios from phase2/per_acq13/pairs.tsv), same definitions as code/collect_all.py."""
import csv, json, math, statistics
from pathlib import Path
M = Path(__file__).resolve().parents[1]; S = M.parent / 'second_engine_all'
med = lambda x: statistics.median([v for v in x if v is not None and math.isfinite(v)])
prod = [r for r in csv.DictReader(open(S / 'results_paired.tsv'), delimiter='\t') if r['cohort'] == 'CAPSCAN']
pairs = {(r['sample'], r['test_arm'], int(r['min_length'])): r for r in csv.DictReader(open(M / 'phase2/per_acq13/pairs.tsv'), delimiter='\t')}
new = [a['sample'] for a in csv.DictReader(open(M / 'inputs/ACQ13.tsv'), delimiter='\t')]
rc = lambda p: json.loads(Path(p).read_text())['peptides_q01'] if Path(p).exists() else None
out = []
for fl in ['FL_public', 'FL_cohortMG']:
    rows = [dict(sample=r['sample'], r9=float(r['msf_ratio_9']), r5=float(r['msf_ratio_5']), s9=float(r['sage_ratio_9']), s5=float(r['sage_ratio_5'])) for r in prod if r['fl_arm'] == fl]
    n_old = len(rows); old = dict(msf9=med([r['r9'] for r in rows]), msf5=med([r['r5'] for r in rows]), sage9=med([r['s9'] for r in rows]), sage5=med([r['s5'] for r in rows]), more9=sum(r['r9'] > 1 for r in rows))
    missing = []
    for s in new:
        a = rc(S / f'results/CAPSCAN__{s}__{fl}/RECEIPT.json'); b = rc(M / f'second_engine/results/CAPSCAN__{s}__MG/RECEIPT.json')
        if not a or not b: missing.append(s); continue
        p9, p5 = pairs[(s, fl, 9)], pairs[(s, fl, 5)]
        rows.append(dict(sample=s, r9=a['len9_50'] / b['len9_50'], r5=a['len5_50'] / b['len5_50'], s9=int(p9['n_test']) / int(p9['n_reference']), s5=int(p5['n_test']) / int(p5['n_reference'])))
    out.append(dict(comparison=fl + '/MG', n_old=n_old, n_new=len(rows), missing_new=missing, old=old,
                    new=dict(msf9=med([r['r9'] for r in rows]), msf5=med([r['r5'] for r in rows]), sage9=med([r['s9'] for r in rows]), sage5=med([r['s5'] for r in rows]), more9=sum(r['r9'] > 1 for r in rows))))
json.dump(out, open(M / 'phase2/SECOND_ENGINE_295.json', 'w'), indent=1); print(json.dumps(out, indent=1))
