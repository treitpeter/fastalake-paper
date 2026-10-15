"""Task A: sequence accumulation numbers + per-specimen 9-50 medians + CapScan 347 FL-only extension."""
import csv, hashlib, json, statistics, sys
from pathlib import Path
from common import ACC, A08, OUT, sha, fasta, table, dump

res = {}
man = json.loads((ACC / 'manifest.json').read_text())
# 1. accumulation curve read-outs
curve = []
for c in ['PREDICT', 'CAPSCAN']:
    rows = list(csv.DictReader(open(ACC / 'full' / c / 'accumulation.tsv'), delimiter='\t'))
    print('rows read', c, 'accumulation.tsv', len(rows))
    for arm in ['FL_public', 'MG']:
        r = [x for x in rows if x['arm'] == arm]; n = len(r)
        get = lambda k, col: float(r[k - 1][col])
        d = {'cohort': c, 'arm': arm, 'specimens': n}
        for k in [50, 100, n]:
            d[f'fixed_k{k if k != n else "final"}'] = int(get(k, 'cumulative_sequences_fixed_order'))
            d[f'meanorder_k{k if k != n else "final"}'] = round(get(k, 'mean_cumulative_sequences_random_order'), 1)
        # mean new sequences per added specimen over the last 50 positions of the order-averaged curve
        d['meanorder_increment_last50_per_specimen'] = round((get(n, 'mean_cumulative_sequences_random_order') - get(n - 50, 'mean_cumulative_sequences_random_order')) / 50, 1)
        d['meanorder_increment_first50_per_specimen'] = round(get(50, 'mean_cumulative_sequences_random_order') / 50, 1)
        curve.append(d)
table(OUT / 'A_accumulation_readout.tsv', curve)

# 2. per-specimen peptide IDs, 5-50 and 9-50, from the expected_sets worker.py verified against
spec = {}
nfiles = 0
for r in man['rows']:
    for arm, info in r['arms'].items():
        for L in ['5', '9']:
            e = info['expected_sets'][L]; assert sha(e['path']) == e['sha256'], e['path']
            s = set(Path(e['path']).read_text().splitlines()); nfiles += 1
            spec.setdefault((r['cohort'], arm, r['design']['specimen_id'], L), set()).update(s)
print('expected_sets files read', nfiles)
specrows = []
for (c, arm, sid, L), s in spec.items():
    specrows.append({'cohort': c, 'arm': arm, 'specimen_id': sid, 'len_min': int(L), 'peptide_ids': len(s)})
table(OUT / 'A_per_specimen_ids.tsv', specrows)
# check 5-50 against the existing per_specimen.tsv
chk = 0
for c in ['PREDICT', 'CAPSCAN']:
    for x in csv.DictReader(open(ACC / 'full' / c / 'per_specimen.tsv'), delimiter='\t'):
        assert len(spec[(c, x['arm'], x['specimen_id'], '5')]) == int(x['peptide_ids']); chk += 1
print('5-50 per-specimen counts equal existing per_specimen.tsv:', chk)
med = []
for c in ['PREDICT', 'CAPSCAN']:
    for arm in ['FL_public', 'MG']:
        for L in [5, 9]:
            v = [r['peptide_ids'] for r in specrows if r['cohort'] == c and r['arm'] == arm and r['len_min'] == L]
            med.append({'cohort': c, 'arm': arm, 'window': f'{L}-50', 'specimens': len(v), 'median_peptide_ids_per_specimen': statistics.median(v)})
table(OUT / 'A_per_specimen_medians.tsv', med)

# 3. CapScan FL extension 282 -> 347 (A08_347 razor FASTAs; acquisition level; specimen IDs for the 16 are not in the accumulation manifest)
cap282 = {r['sample'] for r in man['rows'] if r['cohort'] == 'CAPSCAN'}
m347 = json.loads((A08 / 'inputs/manifest.json').read_text()); s347 = [x['sample'] for x in m347['samples']]
assert len(s347) == 347 and cap282 <= set(s347)
u282, u347, u282prod = set(), set(), set(); n9 = {}; recs = 0
for s in s347:
    f = A08 / 'results' / s / 'parsimony' / f'{s}_razor.fasta'
    hs = {hashlib.sha256(q.encode()).digest() for _, q in fasta(f)}; recs += 1
    u347 |= hs
    if s in cap282: u282 |= hs
    n9[s] = sum(1 for _ in open(A08 / 'results' / s / 'length9_50.txt'))
for r in man['rows']:
    if r['cohort'] == 'CAPSCAN':
        u282prod |= {hashlib.sha256(q.encode()).digest() for _, q in fasta(r['arms']['FL_public']['fasta'])}
print('A08_347 razor FASTAs read', recs)
ext = {'FL_union_282_production_fastas': len(u282prod), 'FL_union_282_A08_347_fastas': len(u282), 'FL_union_347_A08_347_fastas': len(u347),
       'added_by_16_participant15': len(u347) - len(u282),
       'median_9_50_per_acquisition_282': statistics.median(n9[s] for s in cap282), 'median_9_50_per_acquisition_347': statistics.median(n9.values()),
       'median_9_50_per_acquisition_new16': statistics.median(n9[s] for s in s347 if s not in cap282)}
res = {'curve': curve, 'medians': med, 'capscan347': ext}
dump(OUT / 'A_RESULT.json', res); print(json.dumps(res, indent=1))
