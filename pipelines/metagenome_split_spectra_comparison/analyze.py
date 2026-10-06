"""Equal-terms comparison on the odd-scan half: FastaLake (FULL / HALF_SELF / HALF_DISJOINT selection)
versus the metagenome database (MG_ODD), plus the MG full-spectra reproduction check.

Statistic: median over acquisitions of 100*(n_FL/n_MG - 1), as for the manuscript gains.
CI: participant-level (cluster) bootstrap, participants resampled with replacement within the
stratum, all acquisitions of a drawn participant kept; 10,000 draws, seed 20261003; percentile 95%.
Acquisitions without a participant ID are their own cluster.
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd

W = Path(__file__).resolve().parent
AUD = Path('LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/split_half_repair/mzml_recovery_2026-09-27/audit/FULL_SOURCE_DATA.tsv')
B = 10000; SEED = 20261003
require_all = '--partial' not in sys.argv

ro = pd.read_csv(W / 'roster.tsv', sep='\t', dtype={'subject_id': str})
fl = pd.read_csv(AUD, sep='\t')
fl = fl.pivot_table(index=['cohort', 'sample'], columns='arm', values=['peptides_5_50', 'peptides_9_50'])
recs = []
for r in ro.itertuples():
    p = W / f'outputs/{r.cohort}/{r.sample}/COMPLETE.json'
    if not p.exists():
        continue
    j = json.loads(p.read_text())
    d = dict(i=r.i, cohort=r.cohort, sample=r.sample, group=r.group,
             subject_id=r.subject_id if isinstance(r.subject_id, str) and r.subject_id not in ('', 'nan') else 'NOID_' + r.sample,
             odd_mzml_identical=j['odd_mzml_identical'])
    for L in ['5', '9']:
        d[f'mg_odd_{L}'] = j['MG_ODD'][f'peptides_{L}_50']
        d[f'mg_full_{L}'] = j['MG_FULL'][f'peptides_{L}_50']
        d[f'mg_full_equals_frozen_{L}'] = j[f'mg_full_equals_frozen_{L}_50']
        d[f'frozen_mg_{L}'] = int(getattr(r, f'frozen_n_mg_{L}'))
        d[f'frozen_fl_{L}'] = int(getattr(r, f'frozen_n_fl_{L}'))
        for arm, k in [('FULL', 'fl_odd_full'), ('HALF_SELF', 'fl_odd_self'), ('HALF_DISJOINT', 'fl_odd_disjoint')]:
            d[f'{k}_{L}'] = int(fl.loc[(r.cohort, r.sample), (f'peptides_{L}_50', arm)])
    recs.append(d)
a = pd.DataFrame(recs)
print('completed', len(a), 'of', len(ro))
if require_all:
    assert len(a) == len(ro) == 578, 'incomplete roster'
for L in ['5', '9']:
    a[f'gain_frozen_fullspectra_{L}'] = 100 * (a[f'frozen_fl_{L}'] / a[f'frozen_mg_{L}'] - 1)
    for k in ['full', 'self', 'disjoint']:
        a[f'gain_odd_{k}_{L}'] = 100 * (a[f'fl_odd_{k}_{L}'] / a[f'mg_odd_{L}'] - 1)
a.to_csv(W / 'PER_ACQUISITION.tsv', sep='\t', index=False)

rng = np.random.default_rng(SEED)
strata = [('PREDICT', a.cohort == 'PREDICT'), ('CAPSCAN_all', a.cohort == 'CAPSCAN'),
          ('CAPSCAN_stool', (a.cohort == 'CAPSCAN') & (a.group == 'stool')),
          ('CAPSCAN_intestinal', (a.cohort == 'CAPSCAN') & (a.group != 'stool'))]
out = []
for name, m in strata:
    g = a[m].reset_index(drop=True)
    if g.empty:
        continue
    subj = g.subject_id.values; us = np.unique(subj)
    members = [np.flatnonzero(subj == s) for s in us]
    draws = [np.concatenate([members[k] for k in rng.integers(0, len(us), len(us))]) for _ in range(B)]
    for L in ['5', '9']:
        for col in ['gain_frozen_fullspectra', 'gain_odd_full', 'gain_odd_self', 'gain_odd_disjoint']:
            v = g[f'{col}_{L}'].values
            bs = np.array([np.median(v[ix]) for ix in draws])
            out.append(dict(stratum=name, min_length=L, comparison=col, n_acquisitions=len(g), n_participants=len(us),
                            median_gain_pct=round(float(np.median(v)), 2),
                            ci95_low=round(float(np.percentile(bs, 2.5)), 2), ci95_high=round(float(np.percentile(bs, 97.5)), 2),
                            n_fl_gt_mg=int((v > 0).sum()), n_fl_lt_mg=int((v < 0).sum()), n_equal=int((v == 0).sum())))
        # paired contrast: disjoint gain minus same-half (self) gain
        v = (g[f'gain_odd_disjoint_{L}'] - g[f'gain_odd_full_{L}']).values
        bs = np.array([np.median(v[ix]) for ix in draws])
        out.append(dict(stratum=name, min_length=L, comparison='disjoint_minus_full_selection_gain_points', n_acquisitions=len(g),
                        n_participants=len(us), median_gain_pct=round(float(np.median(v)), 2),
                        ci95_low=round(float(np.percentile(bs, 2.5)), 2), ci95_high=round(float(np.percentile(bs, 97.5)), 2),
                        n_fl_gt_mg=int((v > 0).sum()), n_fl_lt_mg=int((v < 0).sum()), n_equal=int((v == 0).sum())))
s = pd.DataFrame(out)
s.to_csv(W / 'SUMMARY.tsv', sep='\t', index=False)
print(s.to_string())
san = dict(n=len(a), odd_mzml_identical=int(a.odd_mzml_identical.sum()),
           mg_full_equals_frozen_5=int(a.mg_full_equals_frozen_5.sum()), mg_full_equals_frozen_9=int(a.mg_full_equals_frozen_9.sum()),
           mg_full_count_equals_frozen_5=int((a.mg_full_5 == a.frozen_mg_5).sum()),
           mismatched=a.loc[~(a.mg_full_equals_frozen_5 & a.mg_full_equals_frozen_9), ['cohort', 'sample', 'mg_full_5', 'frozen_mg_5', 'mg_full_9', 'frozen_mg_9']].to_dict('records'),
           bootstrap=dict(B=B, seed=SEED, unit='participant (subject_id)'))
json.dump(san, open(W / 'SANITY.json', 'w'), indent=1)
print(json.dumps(san, indent=1)[:3000])
