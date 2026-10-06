"""D11 cross-fit, final analysis (578 acquisitions). Same definitions as design/analyze.py, plus round-3 context
and solver-perturbation sensitivity. Writes run/PER_ACQUISITION.tsv, run/SUMMARY.tsv, run/SUMMARY_CAPSCAN_BY_REGION.tsv.

Per acquisition and length view L in {5, 9} (5-50, nested 9-50):
  gain_cf_pooled    100*(FL_pooled/MG_pooled-1)   PRIMARY: 10-fold cross-fit FL vs MG, one LDA refit + picked-peptide q per arm over pooled fold PSMs
  gain_cf_union     100*(FL_union/MG_union-1)     sensitivity: union of per-fold accepted sets (each fold's own Sage peptide_q)
  gain_cf_pert1/2   same as primary with both arms' pooled LDA inputs perturbed by 1 ulp (seed 1/2): solver sensitivity
  gain_full_pooled  100*(FULL_pooled/MG_pooled-1) archived all-prediction FASTA, same folds and pooling (same-spectrum reference)
  gain_full_union   same, union
  cf_vs_full_pooled 100*(FL_pooled/FULL_pooled-1) evidence + same-spectrum loss of cross-fitting
  calib_full_minus_frozen  gain_full_pooled - gain_frozen (pooling calibration)
  gain_frozen       frozen manuscript full-spectra gain (roster)
  round-3 (iteration3/disjoint, odd-scan spectra, native Sage q): gain_odd_full, gain_odd_self, gain_odd_disjoint
Median per stratum; 95% CI by participant cluster bootstrap, 10,000 draws, seed 20261004 (one RNG stream, strata in
fixed order, the same resampled index sets used for every column within a stratum). Acquisitions without a participant ID
are their own cluster.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

RUN = Path(__file__).resolve().parent; CF = RUN.parent; IT3 = Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03') / 'iteration3/disjoint'  # public347
B = 10000; SEED = 20261004
ro = pd.read_csv(IT3 / 'roster.tsv', sep='\t', dtype={'subject_id': str})
r3 = pd.read_csv(IT3 / 'PER_ACQUISITION.tsv', sep='\t', dtype={'subject_id': str}).set_index('i')
recs = []
for r in ro.itertuples():
    j = json.loads((CF / f'outputs/{r.cohort}/{r.sample}/COMPLETE.json').read_text())
    d = dict(i=r.i, cohort=r.cohort, sample=r.sample, group=r.group,
             subject_id=r.subject_id if isinstance(r.subject_id, str) and r.subject_id not in ('', 'nan') else 'NOID_' + r.sample,
             fold_ms2_total=sum(f['fold_ms2'] for f in j['folds'].values()),
             fl_fasta_sequences_median=float(np.median([f['fl_fasta_sequences'] for f in j['folds'].values()])),
             lda_fallback=';'.join(a for a in ['FL', 'FULL', 'MG'] if j[f'{a}_pooled']['rescoring'] != 'lda'))
    for L in ['5', '9']:
        k = f'peptides_{L}_50'
        for arm in ['FL', 'FULL', 'MG']:
            for how in ['pooled', 'union']:
                d[f'{arm}_{how}_{L}'] = j[f'{arm}_{how}'][k]
            for s in (0, 1):
                d[f'{arm}_pert{s + 1}_{L}'] = j[f'{arm}_pooled']['solver_perturbation'][s][k]
            d[f'{arm}_pooled_naive_{L}'] = j[f'{arm}_pooled_naive_diagnostic'][k]
        d[f'frozen_fl_{L}'] = int(getattr(r, f'frozen_n_fl_{L}')); d[f'frozen_mg_{L}'] = int(getattr(r, f'frozen_n_mg_{L}'))
        g = lambda a, b: 100 * (a / b - 1) if b > 0 else np.nan
        d[f'gain_cf_pooled_{L}'] = g(d[f'FL_pooled_{L}'], d[f'MG_pooled_{L}'])
        d[f'gain_cf_union_{L}'] = g(d[f'FL_union_{L}'], d[f'MG_union_{L}'])
        d[f'gain_cf_pert1_{L}'] = g(d[f'FL_pert1_{L}'], d[f'MG_pert1_{L}'])
        d[f'gain_cf_pert2_{L}'] = g(d[f'FL_pert2_{L}'], d[f'MG_pert2_{L}'])
        d[f'gain_full_pooled_{L}'] = g(d[f'FULL_pooled_{L}'], d[f'MG_pooled_{L}'])
        d[f'gain_full_union_{L}'] = g(d[f'FULL_union_{L}'], d[f'MG_union_{L}'])
        d[f'gain_frozen_{L}'] = g(d[f'frozen_fl_{L}'], d[f'frozen_mg_{L}'])
        d[f'cf_vs_full_pooled_{L}'] = g(d[f'FL_pooled_{L}'], d[f'FULL_pooled_{L}'])
        d[f'cf_vs_full_union_{L}'] = g(d[f'FL_union_{L}'], d[f'FULL_union_{L}'])
        d[f'calib_full_minus_frozen_{L}'] = d[f'gain_full_pooled_{L}'] - d[f'gain_frozen_{L}']
        d[f'mg_pooled_vs_frozen_{L}'] = g(d[f'MG_pooled_{L}'], d[f'frozen_mg_{L}'])
        x = r3.loc[r.i]; assert x['sample'] == r.sample
        for c in ['gain_odd_full', 'gain_odd_self', 'gain_odd_disjoint']:
            d[f'{c}_{L}'] = x[f'{c}_{L}']
    recs.append(d)
a = pd.DataFrame(recs)
assert len(a) == 578 and a.groupby('cohort').size().to_dict() == {'PREDICT': 296, 'CAPSCAN': 282}
a.to_csv(RUN / 'PER_ACQUISITION.tsv', sep='\t', index=False)

COLS = ['gain_cf_pooled', 'gain_cf_union', 'gain_cf_pert1', 'gain_cf_pert2', 'gain_full_pooled', 'gain_full_union', 'gain_frozen',
        'gain_odd_full', 'gain_odd_self', 'gain_odd_disjoint', 'cf_vs_full_pooled', 'cf_vs_full_union', 'calib_full_minus_frozen',
        'mg_pooled_vs_frozen']


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
                                n_pos=int((v > 0).sum()), n_neg=int((v < 0).sum()),
                                n_zero=int((v == 0).sum()), n_nan=int(np.isnan(v).sum())))
    return pd.DataFrame(out)


rng = np.random.default_rng(SEED)
cap = a.cohort == 'CAPSCAN'
s = summarize([('PREDICT', a.cohort == 'PREDICT'), ('CAPSCAN_all', cap), ('CAPSCAN_stool', cap & (a.group == 'stool')),
               ('CAPSCAN_intestinal', cap & (a.group != 'stool'))], rng)
s.to_csv(RUN / 'SUMMARY.tsv', sep='\t', index=False)
reg = summarize([(f'CAPSCAN_{g}', cap & (a.group == g)) for g in ['duodenum', 'jejunum', 'ileum', 'colon', 'stool']], np.random.default_rng(SEED))
reg.to_csv(RUN / 'SUMMARY_CAPSCAN_BY_REGION.tsv', sep='\t', index=False)
pd.set_option('display.width', 250)
print(s.pivot_table(index=['stratum', 'min_length'], columns='comparison', values='median', sort=False).to_string())
print(reg[reg.comparison.isin(['gain_cf_pooled', 'gain_cf_union', 'gain_odd_disjoint'])].to_string())
print('lda fallbacks:', a.lda_fallback.ne('').sum())
