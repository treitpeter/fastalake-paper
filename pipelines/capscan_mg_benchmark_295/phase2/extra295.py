#!/usr/bin/env python3
"""(1) same-spectrum mixed model at 282 and 295 with Fig4 panels.stat_mixed_model(view='same');
(2) crossfit_summary.tsv / crossfit_summary_by_region.tsv regenerated at 591 rows with the summarize() of
public347 run/analyze_final.py (copied verbatim below; B=10000, seed 20261004)."""
import sys, json
from pathlib import Path
import numpy as np, pandas as pd
P = Path(__file__).resolve().parent
res = {}
for name in ['repo_base', 'repo']:
    U = P / name / 'figures/Fig4'; sys.path[:0] = [str(U.parents[1]), str(U)]
    for m in [k for k in list(sys.modules) if k in ('panels', 'style', 'style.fl_style')]: del sys.modules[m]
    import os; os.chdir(U)
    import panels
    x = panels.load_crossfit()
    res[name] = {f'{v}_{L}': panels.stat_mixed_model(x, view=v, L=L) for v in ('crossfit', 'same') for L in ('9', '5')}
    sys.path[:] = [p for p in sys.path if not p.startswith(str(P / name))]
json.dump(res, open(P / 'MIXED_MODEL_282_vs_295.json', 'w'), indent=1, default=float)
for k in res['repo']:
    a, b = res['repo_base'][k], res['repo'][k]
    print(k, 'n', a['n_acq'], '->', b['n_acq'], '| fold', round(a['fold'], 3), '->', round(b['fold'], 3), '| CI', round(a['fold_ci_low'], 3), round(a['fold_ci_high'], 3), '->', round(b['fold_ci_low'], 3), round(b['fold_ci_high'], 3), '| P', '%.2g' % a['p'], '->', '%.2g' % b['p'], '| region P', round(a['region_p'], 3), '->', round(b['region_p'], 3))

a = pd.read_csv(P / 'repo/source_data/supplementary_methods_6_7/crossfit_per_acquisition.tsv', sep='\t', dtype={'subject_id': str})
assert len(a) == 591
B = 10000; SEED = 20261004
COLS = ['gain_cf_pooled', 'gain_cf_union', 'gain_cf_pert1', 'gain_cf_pert2', 'gain_full_pooled', 'gain_full_union', 'gain_frozen',
        'gain_odd_full', 'gain_odd_self', 'gain_odd_disjoint', 'cf_vs_full_pooled', 'cf_vs_full_union', 'calib_full_minus_frozen', 'mg_pooled_vs_frozen']
def summarize(strata, rng):
    out = []
    for name, m in strata:
        g = a[m].reset_index(drop=True)
        subj = g.subject_id.values; us = np.unique(subj)
        members = [np.flatnonzero(subj == s) for s in us]
        draws = [np.concatenate([members[k] for k in rng.integers(0, len(us), len(us))]) for _ in range(B)]
        for L in ['5', '9']:
            for col in COLS:
                v = g[f'{col}_{L}'].values.astype(float)
                bs = np.array([np.nanmedian(v[ix]) for ix in draws])
                out.append(dict(stratum=name, min_length=L, comparison=col, n_acquisitions=len(g), n_participants=len(us),
                                median=round(float(np.nanmedian(v)), 2), ci95_low=round(float(np.percentile(bs, 2.5)), 2),
                                ci95_high=round(float(np.percentile(bs, 97.5)), 2), median_exact=float(np.nanmedian(v)),
                                ci95_low_exact=float(np.percentile(bs, 2.5)), ci95_high_exact=float(np.percentile(bs, 97.5)),
                                n_pos=int((v > 0).sum()), n_neg=int((v < 0).sum()), n_zero=int((v == 0).sum()), n_nan=int(np.isnan(v).sum())))
    return pd.DataFrame(out)
cap = a.cohort == 'CAPSCAN'
s = summarize([('PREDICT', a.cohort == 'PREDICT'), ('CAPSCAN_all', cap), ('CAPSCAN_stool', cap & (a.group == 'stool')),
               ('CAPSCAN_intestinal', cap & (a.group != 'stool'))], np.random.default_rng(SEED))
reg = summarize([(f'CAPSCAN_{g}', cap & (a.group == g)) for g in ['duodenum', 'jejunum', 'ileum', 'colon', 'stool']], np.random.default_rng(SEED))
s.to_csv(P / 'crossfit_summary_591.tsv', sep='\t', index=False); reg.to_csv(P / 'crossfit_summary_by_region_591.tsv', sep='\t', index=False)
print(s[(s.min_length == '9') & s.comparison.isin(['gain_cf_pooled', 'gain_frozen', 'gain_odd_disjoint', 'gain_odd_self'])][['stratum', 'comparison', 'n_acquisitions', 'median', 'ci95_low', 'ci95_high', 'n_pos']].to_string())
print(reg[(reg.min_length == '9') & reg.comparison.isin(['gain_cf_pooled', 'gain_frozen'])][['stratum', 'comparison', 'n_acquisitions', 'median', 'ci95_low', 'ci95_high']].to_string())
