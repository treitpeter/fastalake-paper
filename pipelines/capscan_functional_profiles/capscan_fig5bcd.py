"""Fig. 5b-d CapScan statistics: copy of biology_refit/capscan_functions.py (refit mode) and regional_host.py
with the participant count generalized from the hard-coded 14 to the roster.

Unchanged computation: acquisitions with positive OG profile, relative profiles averaged within specimen
then participant-region; Hellinger own-stool test (99,999 permutations, seed 20260910); orthogonal
participant/region partition (manuscript_benchmark/common.py partition with 14 -> n); human-assigned share
from taxonomy_signal / total_LFQ; exact sign-flip over 2^n assignments, BH over four regions.
Added: 2,000-resample participant bootstrap CIs (seed 20261003) as printed in v25.
usage: capscan_fig5bcd.py <profiles_dir> <design.tsv> <out_dir>
"""
import json, sys, hashlib
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse
from scipy.spatial.distance import cdist
from scipy.stats import false_discovery_control

REG = ['duodenum', 'jejunum', 'ileum', 'colon', 'stool']
prof, design_path, out = Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=False)
def frame(folder, key):
    a = json.loads((folder/(key + '_axes.json')).read_text()); v = sparse.load_npz(folder/(key + '.npz'))
    return pd.DataFrame(v.toarray(), index=a['samples'], columns=a['features'])
def partition(h, n):
    cube = np.asarray(h).reshape(n, 5, -1); rows = []
    for label, xx in [('All five regions', cube), ('Capsules only', cube[:, :4])]:
        nn, r, k = xx.shape; grand = xx.mean((0, 1)); pm = xx.mean(1) - grand; rm = xx.mean(0) - grand
        res = xx - grand - pm[:, None] - rm[None]; ss = ((xx - grand)**2).sum()
        d = {'Participant': r*(pm**2).sum()/ss, 'Residual': (res**2).sum()/ss}
        if r == 5:
            means = xx.mean(0); contrast = np.array([1, 1, 1, 1, -4])/np.sqrt(20)
            d['Capsules versus stool'] = nn*((contrast[:, None]*means).sum(0)**2).sum()/ss
            d['Between capsule regions'] = nn*((means[:4] - means[:4].mean(0))**2).sum()/ss
        else: d.update({'Capsules versus stool': 0., 'Between capsule regions': nn*(rm**2).sum()/ss})
        assert np.isclose(sum(d.values()), 1, atol=1e-10)
        rows.extend(dict(population=label, component=k, fraction=float(v)) for k, v in d.items())
    return pd.DataFrame(rows)
def boot(x, seed=20261003, B=2000):
    x = np.asarray(x); rng = np.random.default_rng(seed); idx = rng.integers(0, len(x), (B, len(x)))
    m = x[idx].mean(1); return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))

design = pd.read_csv(design_path, sep='\t').set_index('sample_id')
og = frame(prof, 'OG_other_unresolved'); assert set(design.index) == set(og.index), 'profile/design mismatch'
roster = design.index[og.loc[design.index].sum(1) > 0]
# production: participants 1-14 (range(1,15)) = those observed in all five groups; generalized here
obs = design.loc[roster].groupby('participant').region.nunique()
parts = sorted(obs.index[obs.eq(5)]); n = len(parts)
ix = pd.MultiIndex.from_product([parts, REG], names=['participant', 'region'])
d = design.loc[roster]
def cells(v):
    rel = v.div(v.sum(1), axis=0); return rel.groupby([d.participant, d.region, d.specimen_id]).mean().groupby(level=[0, 1]).mean().loc[ix]
cc = cells(og.loc[roster]); x = cc.to_numpy()
perms = np.argsort(np.random.default_rng(20260910).random((99999, n)), axis=1)
rows, don, parts_out = [], [], []
for prev in [0., .1, .2]:
    keep = (x > 0).mean(0) >= prev; v = x[:, keep]; h = np.sqrt(v/v.sum(1, keepdims=True)); cube = h.reshape(n, 5, -1)
    dd = np.stack([cdist(cube[i, :4], cube[:, 4]).mean(0) for i in range(n)]); own = np.diag(dd); other = (dd.sum(1) - own)/(n - 1)
    obs = float((other - own).mean()); pdg = dd[np.arange(n)[None, :], perms].mean(1); null = (dd.sum(1).mean() - n*pdg)/(n - 1)
    p = (1 + int((abs(null) >= abs(obs) - 1e-14).sum()))/100000
    lo, hi = boot(other - own)
    rows.append(dict(prevalence=prev, participants=n, features=int(keep.sum()), positive_features=int((v > 0).any(0).sum()), own_closer=int((own < other).sum()),
                     mean_own=float(own.mean()), mean_other=float(other.mean()), mean_advantage=obs, ci_low=lo, ci_high=hi, p=p))
    don.extend(dict(prevalence=prev, participant=parts[i], own_stool_distance=own[i], other_stool_mean_distance=other[i]) for i in range(n))
    pt = partition(h, n); pt['prevalence'] = prev; parts_out.append(pt)
pd.DataFrame(rows).to_csv(out/'similarity_summary.tsv', sep='\t', index=False)
pd.DataFrame(don).to_csv(out/'own_stool_participants.tsv', sep='\t', index=False)
pd.concat(parts_out).to_csv(out/'orthogonal_partition.tsv', sep='\t', index=False)
# Fig 5b: human-assigned share (capscan_functions refit: taxonomy_signal / total_LFQ, same averaging)
tx = pd.read_csv(prof/'taxonomy_signal.tsv', sep='\t').set_index('sample')
tf = tx.loc[roster].drop(columns='total_LFQ').div(tx.loc[roster, 'total_LFQ'], axis=0)
tc = tf.groupby([d.participant, d.region, d.specimen_id]).mean().groupby(level=[0, 1]).mean().loc[ix]
tc.to_csv(out/'origin_fraction_cells.tsv', sep='\t')
w = tc['Human-assigned'].unstack().loc[parts, REG]*100
signs = 2*((np.arange(2**n, dtype=np.uint64)[:, None] >> np.arange(n, dtype=np.uint64)) & 1).astype(np.int8) - 1
tests = []
for reg in REG[:4]:
    delta = (w[reg] - w.stool).to_numpy(); p = float(np.mean(abs(signs @ delta) >= abs(delta.sum()) - 1e-12)); lo, hi = boot(delta)
    tests.append(dict(region=reg, n=n, mean_region=float(w[reg].mean()), mean_stool=float(w.stool.mean()), mean_difference=float(delta.mean()), ci_low=lo, ci_high=hi, p=p))
for r, q in zip(tests, false_discovery_control([t['p'] for t in tests])): r['q'] = float(q)
pd.DataFrame(tests).to_csv(out/'four_region_tests.tsv', sep='\t', index=False)
w.reset_index().to_csv(out/'patient_region_human_percent.tsv', sep='\t', index=False)
json.dump(dict(status='PASS', participants=n, acquisitions=int(len(roster)), design_rows=len(design), cells=len(ix),
               script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), open(out/'COMPLETE.json', 'w'), indent=2)
print(pd.DataFrame(rows).to_string(index=False)); print(pd.DataFrame(tests).to_string(index=False))
print(pd.concat(parts_out).query('prevalence==0').to_string(index=False))
