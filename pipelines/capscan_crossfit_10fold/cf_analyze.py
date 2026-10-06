"""Summarize C1 cross-fit per_acq/*.json with the iteration4 analyze_final.py statistics
(median; participant-cluster bootstrap 10,000 draws, seed 20261004; NOID_ clusters)."""
import json, glob, os
from pathlib import Path
import numpy as np, pandas as pd
H = Path(os.environ.get('FASTALAKE_OUT_DIR', Path(__file__).resolve().parent)); B = 10000; SEED = 20261004
recs = []
for f in sorted(glob.glob(str(H / 'per_acq/*.json'))):
    j = json.loads(open(f).read())
    sid = j['subject_id'] if isinstance(j['subject_id'], str) and j['subject_id'] not in ('', 'nan') else 'NOID_' + j['sample']
    d = dict(i=j['i'], cohort=j['cohort'], sample=j['sample'], group=j['group'], subject_id=sid,
             baseline_equal_all=all(j[a]['baseline_equals_frozen'] for a in ['FL', 'FULL', 'MG']),
             fallback=';'.join(a for a in ['FL', 'FULL', 'MG'] if j[a]['rescoring'] != 'lda'))
    g = lambda a, b: 100 * (a / b - 1) if b > 0 else np.nan
    for L in ['5', '9']:
        k = f'peptides_{L}_50'
        for opt in ['baseline', 'opt_a', 'opt_b']:
            for arm in ['FL', 'FULL', 'MG']:
                d[f'{arm}_{opt}_{L}'] = j[arm][opt][k]
            d[f'gain_cf_{opt}_{L}'] = g(d[f'FL_{opt}_{L}'], d[f'MG_{opt}_{L}'])
            d[f'cf_vs_full_{opt}_{L}'] = g(d[f'FL_{opt}_{L}'], d[f'FULL_{opt}_{L}'])
    recs.append(d)
a = pd.DataFrame(recs)
assert len(a) == 578 and a.groupby('cohort').size().to_dict() == {'PREDICT': 296, 'CAPSCAN': 282}, a.groupby('cohort').size()
a.to_csv(H / 'CF_PEP_PER_ACQUISITION.tsv', sep='\t', index=False)
COLS = [f'{c}_{o}' for c in ['gain_cf', 'cf_vs_full'] for o in ['baseline', 'opt_a', 'opt_b']]
def summarize(strata, rng):
    """Median gain per stratum with participant-cluster bootstrap 95% interval."""
    out = []
    for name, m in strata:
        g = a[m].reset_index(drop=True); subj = g.subject_id.values; us = np.unique(subj)
        members = [np.flatnonzero(subj == s) for s in us]
        draws = [np.concatenate([members[k] for k in rng.integers(0, len(us), len(us))]) for _ in range(B)]
        for L in ['5', '9']:
            for col in COLS:
                v = g[f'{col}_{L}'].values.astype(float); bs = np.array([np.nanmedian(v[ix]) for ix in draws])
                out.append(dict(stratum=name, min_length=L, comparison=col, n_acquisitions=len(g), n_participants=len(us),
                                median=round(float(np.nanmedian(v)), 2), ci95_low=round(float(np.percentile(bs, 2.5)), 2),
                                ci95_high=round(float(np.percentile(bs, 97.5)), 2), n_pos=int((v > 0).sum()), n_neg=int((v < 0).sum()),
                                n_zero=int((v == 0).sum()), n_nan=int(np.isnan(v).sum())))
    return pd.DataFrame(out)
cap = a.cohort == 'CAPSCAN'
s = summarize([('PREDICT', a.cohort == 'PREDICT'), ('CAPSCAN_all', cap), ('CAPSCAN_stool', cap & (a.group == 'stool')),
               ('CAPSCAN_intestinal', cap & (a.group != 'stool'))], np.random.default_rng(SEED))
s.to_csv(H / 'CF_PEP_SUMMARY.tsv', sep='\t', index=False)
reg = summarize([(f'CAPSCAN_{g}', cap & (a.group == g)) for g in ['duodenum', 'jejunum', 'ileum', 'colon', 'stool']], np.random.default_rng(SEED))
reg.to_csv(H / 'CF_PEP_SUMMARY_BY_REGION.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print(s.pivot_table(index=['stratum', 'min_length'], columns='comparison', values='median', sort=False).to_string())
print('baseline != frozen acquisitions:', int((~a.baseline_equal_all).sum()), ' fallbacks:', int(a.fallback.ne('').sum()))
