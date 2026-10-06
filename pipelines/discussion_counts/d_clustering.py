"""Task D: CAMPI clustering cost at 9-50 from the same Sage outputs FINAL/clustering_check/recount.py read."""
import csv, json, statistics
from pathlib import Path
from common import FINAL, LAB, OUT, accepted, table, dump

RUN = LAB / 'review/v15_integration/analyses/clustering_90_95_2026-09-29'
ARMS = ['unclustered', 'cluster_099', 'cluster_097', 'cluster_095', 'cluster_090', 'cluster_075', 'cluster_050']
OLD = {'unclustered', 'cluster_097', 'cluster_099'}
m = json.loads((RUN / 'manifest.json').read_text())
prev = {(r['sample'], r['arm']): int(r['peptides']) for r in csv.DictReader(open(FINAL / 'clustering_check/per_acquisition_recount.tsv'), delimiter='\t')}
rows = []; sets = {}; nrows = 0
for t in m['tasks']:
    S = t['sample']; prior = Path(t['prior_search'])
    for arm in ARMS:
        d = prior / arm if arm in OLD else RUN / 'results' / S / 'searches' / arm
        s5, s9, n = accepted(d / 'results.sage.tsv'); nrows += n
        assert len(s5) == prev[S, arm], (S, arm)
        sets[S, arm] = (s5, s9)
for t in m['tasks']:
    S = t['sample']
    for arm in ARMS:
        r = {'sample': S, 'arm': arm}
        for i, L in enumerate([5, 9]):
            b = len(sets[S, 'unclustered'][i]); v = len(sets[S, arm][i])
            r[f'peptides_{L}_50'] = v; r[f'reduction_{L}_50_pct'] = 100 * (1 - v / b)
        rows.append(r)
table(OUT / 'D_per_acquisition.tsv', rows)
summ = []
for arm in ARMS[1:]:
    R = [r for r in rows if r['arm'] == arm]
    summ.append({'arm': arm, 'n': len(R), 'median_reduction_5_50_pct': statistics.median(r['reduction_5_50_pct'] for r in R),
                 'median_reduction_9_50_pct': statistics.median(r['reduction_9_50_pct'] for r in R),
                 'n_below_baseline_9_50': sum(r['reduction_9_50_pct'] > 0 for r in R)})
table(OUT / 'D_summary.tsv', summ)
dump(OUT / 'D_RESULT.json', {'psm_rows_read': nrows, 'searches': len(rows), 'five_50_equals_clustering_check': True, 'summary': summ})
print('PSM rows read', nrows, 'searches', len(rows)); print(json.dumps(summ, indent=1))
