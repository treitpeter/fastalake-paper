#!/usr/bin/env python3
"""Check current held-out and error-control summaries from acquisition-level inputs.

The default checks all medians, counts and paired-estimator denominators.
--bootstrap also repeats the 10,000 participant-cluster draws for group-fold and
entrapment intervals, in the original stratum/cluster order (seed 20261004).
Searches and peptide assignments are upstream inputs, not repeated here.
"""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

SD = Path(__file__).resolve().parents[1] / 'source_data/v42_updates'


def read(name):
    return pd.read_csv(SD / name, sep='\t', float_precision='round_trip')


def equal(actual, expected, label):
    assert np.isclose(actual, expected, rtol=1e-12, atol=0, equal_nan=True), (label, actual, expected)


def strata(a):
    cap = a.cohort.eq('CAPSCAN')
    return [('PREDICT', a.cohort.eq('PREDICT')), ('CAPSCAN_all', cap),
            ('CAPSCAN_stool', cap & a.group.eq('stool')),
            ('CAPSCAN_intestinal', cap & ~a.group.eq('stool'))]


def validate_roster(a):
    assert len(a) == 591 and not a.duplicated(['cohort', 'sample']).any()
    assert a.groupby('cohort').size().to_dict() == {'CAPSCAN': 295, 'PREDICT': 296}
    assert a.loc[a.cohort.eq('PREDICT'), 'sample'].str.fullmatch(r'PREDICT_A\d{3}').all()
    assert a.loc[a.cohort.eq('PREDICT'), 'group'].isna().all()


def draws_for(g, rng):
    members = [np.flatnonzero(g.cluster_order.to_numpy() == s) for s in np.unique(g.cluster_order)]
    return [np.concatenate([members[k] for k in rng.integers(0, len(members), len(members))]) for _ in range(10000)]


def summary_checks(v, row, prefix, draws, scale=1):
    equal(scale * np.nanmedian(v), row['median' + prefix], 'median')
    if draws is not None:
        boot = np.array([np.nanmedian(v[ix]) for ix in draws])
        for name, q in [('ci95_low', 2.5), ('ci95_high', 97.5)]:
            equal(scale * np.nanpercentile(boot, q), row[name + prefix], name)
    if 'n_nan' in row and pd.notna(row['n_nan']):
        assert int(np.isnan(v).sum()) == row['n_nan']


def check_group_folds(bootstrap=False):
    a, s = read('group_folds_per_acquisition.tsv'), read('group_folds_summary.tsv')
    validate_roster(a)
    for design in ('group', 'native'):
        for length in (5, 9):
            computed = 100 * (a[f'FL_{design}_{length}'] / a[f'MG_{design}_{length}'].replace(0, np.nan) - 1)
            np.testing.assert_allclose(computed, a[f'gain_{design}_{length}'], rtol=1e-12, atol=0, equal_nan=True)
    rng = np.random.default_rng(20261004)
    for name, mask in strata(a):
        g = a[mask].reset_index(drop=True)
        draws = draws_for(g, rng) if bootstrap else None
        for row in s[s.stratum.eq(name)].to_dict('records'):
            assert row['n'] == len(g)
            q = row['quantity']
            if q.startswith('n_acq_FL_more_than_MG_'):
                design, length = q.removeprefix('n_acq_FL_more_than_MG_').rsplit('_', 1)
                assert row['median'] == int((g[f'FL_{design}_{length}'] > g[f'MG_{design}_{length}']).sum())
                continue
            assert row['n_participants'] == g.cluster_order.nunique()
            v = g[q].to_numpy(float)
            summary_checks(v, row, '', draws)
            for field, count in [('n_pos', (v > 0).sum()), ('n_zero', (v == 0).sum()), ('n_neg', (v < 0).sum())]:
                assert row[field] == count
    return len(s)


def check_entrapment(bootstrap=False):
    a, s = read('entrapment_per_acquisition.tsv'), read('entrapment_summary.tsv')
    validate_roster(a)
    for arm in ('FL', 'MG', 'FULL', 'GAINED'):
        for gate in ('q', 'q_pep'):
            for length in ('5_50', '5_8', '9_50'):
                pre = f'{arm}_{gate}_{length}'
                n = a[f'{pre}_n_accepted'].replace(0, np.nan)
                for estimator, num in [('paired_estimate', a[f'{pre}_paired_numerator']),
                                       ('diagnostic_with_unresolved', a[f'{pre}_paired_numerator'] + a[f'{pre}_n_entrapment_unresolved'])]:
                    np.testing.assert_allclose(num / n, a[f'{pre}_{estimator}'], rtol=1e-12, atol=0, equal_nan=True)
    rng = np.random.default_rng(20261004)
    for name, mask in strata(a):
        g = a[mask].reset_index(drop=True)
        draws = draws_for(g, rng) if bootstrap else None
        for row in s[s.stratum.eq(name)].to_dict('records'):
            assert row['n_acquisitions'] == len(g) and row['n_participants'] == g.cluster_order.nunique()
            arm, gate, length, est = (row[k] for k in ('arm', 'gate', 'length', 'estimator'))
            pre = f'{arm}_{gate}_{length}'
            scale = 100
            if arm == 'FL_minus_MG':
                v = g[f'FL_{gate}_{length}_paired_estimate'] - g[f'MG_{gate}_{length}_paired_estimate']
            elif arm == 'FL_vs_MG_target_gain_with_entrapment':
                v = 100 * (g[f'FL_{gate}_{length}_n_target'] / g[f'MG_{gate}_{length}_n_target'].replace(0, np.nan) - 1)
                scale = 1
            else:
                v = g[f'{pre}_{est}']
                if est in ('paired_estimate', 'diagnostic_with_unresolved'):
                    numerator = g[f'{pre}_paired_numerator'].sum()
                    if est == 'diagnostic_with_unresolved':
                        numerator += g[f'{pre}_n_entrapment_unresolved'].sum()
                    equal(100 * numerator / g[f'{pre}_n_accepted'].sum(), row['pooled_sum_pct'], 'pooled estimate')
            summary_checks(v.to_numpy(float), row, '_pct', draws, scale)
    return len(s)


def check_native_crossfit():
    a = read('crossfit_per_acquisition.tsv')
    validate_roster(a)
    n = 0
    for file in ('crossfit_summary.tsv', 'crossfit_summary_by_region.tsv'):
        for row in read(file).to_dict('records'):
            name = row['stratum']
            mask = dict(strata(a)).get(name)
            if mask is None:
                assert name.startswith('CAPSCAN_')
                mask = a.cohort.eq('CAPSCAN') & a.group.eq(name.removeprefix('CAPSCAN_'))
            g = a[mask]
            assert row['n_acquisitions'] == len(g) and row['n_participants'] == g.cluster_order.nunique()
            v = g[f"{row['comparison']}_{row['min_length']}"].to_numpy(float)
            equal(np.nanmedian(v), row['median_exact'], name)
            n += 1
    return n


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--bootstrap', action='store_true')
    args = parser.parse_args()
    print('Group-fold summaries:', check_group_folds(args.bootstrap))
    print('Entrapment summaries:', check_entrapment(args.bootstrap))
    print('Current native cross-fit summaries:', check_native_crossfit())
    print('PASS' + (' including participant-bootstrap intervals' if args.bootstrap else ' (intervals retained from upstream participant bootstrap)'))


if __name__ == '__main__':
    main()
