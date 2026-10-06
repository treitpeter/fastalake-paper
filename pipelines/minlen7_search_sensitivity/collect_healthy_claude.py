"""Healthy-stool own / common-FASTA / joint strategies at 5-50 (published) and fresh >=7.
own: wave-1 HSTOOL FL collector sets (baseline5 / fresh7). common: A06 study-union search (baseline, from tasks.json source) vs wave-2 output.
joint: A12 joint search (baseline) vs wave-2 joint output, split per acquisition by `filename`.
Gate for every set: rank 1, label 1, peptide_q <= 0.01, canonical (A10 peptide.canonical: mods stripped, I->L), length window.
Reproduction gate: medians must equal P48 (common +7.1% count, 88.8% own identities; joint +32.9%, 92.5%) to 1 decimal, else FAIL."""
import json, csv, gzip, os, sys, statistics, collections
from pathlib import Path
F = Path('LAB_ROOT')
sys.path.insert(0, str(F/'MM_review_analysis/A10_4_3_4_5_identity_recovery/frozen')); from peptide import canonical
csv.field_size_limit(1 << 30)
W2 = Path(__file__).resolve().parent.parent; RUN = W2.parent
tasks = json.load(open(W2/'tasks.json'))
def gate(path, lo, split=False):
    out = collections.defaultdict(set) if split else set()
    for r in csv.DictReader(open(path), delimiter='\t'):
        if r['rank'] != '1' or r['label'] != '1' or float(r['peptide_q']) > .01: continue
        s = canonical(r['peptide'])
        if not lo <= len(s) <= 50: continue
        if split: out[Path(r['filename']).stem].add(s)
        else: out.add(s)
    return out
w1 = {r['sample']: r for r in csv.DictReader(open(RUN/'tasks.tsv'), delimiter='\t') if r['cohort'] == 'HSTOOL'}
own = {}
for s, r in w1.items():
    p = RUN/'wave2/peptide_tasks'/f"{int(r['index']):04d}"
    own[s] = {k: set(gzip.open(p/f'{k}.txt.gz', 'rt').read().split()) for k in ('baseline5', 'fresh7')}
common = {}
for t in tasks:
    if t['kind'] != 'common': continue
    rec = json.load(open(W2/'receipts'/f"{t['index']:03d}.json")); assert rec['status'] == 'PASS'
    common[t['sample']] = dict(baseline5=gate(Path(t['source']).parent/'results.sage.tsv', 5), fresh7=gate(Path(t['output'])/'results.sage.tsv', 7))
jt = [t for t in tasks if t['kind'] == 'joint'][0]
assert json.load(open(W2/'receipts'/f"{jt['index']:03d}.json"))['status'] == 'PASS'
jb = gate(Path(jt['source']).parent/'results.sage.tsv', 5, split=True); jf = gate(Path(jt['output'])/'results.sage.tsv', 7, split=True)
names = sorted(own); assert len(names) == 58 and set(common) == set(names), (len(names), len(common))
assert set(jb) >= set(names) and set(jf) >= set(names), 'joint filename stems do not match acquisitions'
rows = []
for n in names:
    for pop in ('baseline5', 'fresh7'):
        o = own[n][pop]
        for strat, x in (('common', common[n][pop]), ('joint', (jb if pop == 'baseline5' else jf)[n])):
            rows.append(dict(sample=n, population=pop, strategy=strat, own=len(o), test=len(x), count_change_pct=100*(len(x)/len(o)-1), own_identity_recovery_pct=100*len(o & x)/len(o)))
out = W2/'healthy_strategies'; out.mkdir(exist_ok=True)
with open(out/'per_acquisition.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, list(rows[0]), delimiter='\t'); w.writeheader(); w.writerows(rows)
summ = {}
for pop in ('baseline5', 'fresh7'):
    for strat in ('common', 'joint'):
        g = [r for r in rows if r['population'] == pop and r['strategy'] == strat]
        summ[f'{strat}_{pop}'] = dict(n=len(g), median_count_change_pct=statistics.median(r['count_change_pct'] for r in g), median_own_identity_recovery_pct=statistics.median(r['own_identity_recovery_pct'] for r in g))
# study-wide distinct repertoire
for pop, lo in (('baseline5', 5), ('fresh7', 7)):
    summ[f'repertoire_{pop}'] = dict(own=len(set().union(*(own[n][pop] for n in names))), common=len(set().union(*(common[n][pop] for n in names))), joint=len(set().union(*((jb if pop=='baseline5' else jf)[n] for n in names))))
ok = (round(summ['common_baseline5']['median_count_change_pct'], 1) == 7.1 and round(summ['common_baseline5']['median_own_identity_recovery_pct'], 1) == 88.8
      and round(summ['joint_baseline5']['median_count_change_pct'], 1) == 32.9 and round(summ['joint_baseline5']['median_own_identity_recovery_pct'], 1) == 92.5)
summ['reproduces_P48'] = ok
json.dump(summ, open(out/'SUMMARY.json', 'w'), indent=1); print(json.dumps(summ, indent=1))
sys.exit(0 if ok else 1)
