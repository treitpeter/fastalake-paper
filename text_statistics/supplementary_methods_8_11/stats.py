"""Supplementary Methods 8-11: statistics and testing families (v28d lines 779-911).

This module recomputes the statistics and counts that Supplementary Methods 8-11 state
and checks every Benjamini-Hochberg family named in the testing-family map (v28d line
893) by recomputing each adjusted q value from the raw P values within exactly that
family.

Inputs are read only from ``source_data/``:

* ``source_data/supplementary_methods_8_11/`` (this unit): own-stool distance matrices,
  CapScan DNA runs, PREDICT native feature counts, held-out EDS prediction outputs,
  non-specific search settings, healthy-stool grouping receipt.
* tables owned by other units, named by path and never duplicated:
  ``source_data/Fig5/`` (CapScan regional, own-stool, cross-omics 347-acquisition outputs),
  ``source_data/Fig6/`` (histidine and terminal-fraction families),
  ``source_data/SFig10/`` (EDS screen), ``source_data/SFig11/`` (PREDICT cross-omics,
  controlled distance comparisons), ``source_data/SFig14/`` (spectrum strata),
  ``source_data/STab3/`` (antibiotic endpoints), ``source_data/STab5/`` (individual-function
  correlations).

Every function returns plain Python values; the notebook collects them in one dict ``R``
whose keys are the ``key`` column of ``NUMBER_SPEC``.
"""

import itertools
import platform
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.stats import false_discovery_control, spearmanr
from sklearn.model_selection import KFold

# ----------------------------------------------------------------------------- choices
REG = ['duodenum', 'jejunum', 'ileum', 'colon', 'stool']   # regions 1-4 and stool (v28d line 843)
BOOT_SEED_REGIONAL = 20261003      # participant bootstrap seed, regional and own-stool intervals (not stated in the text; producer capscan_fig5bcd.py boot())
BOOT_DRAWS_REGIONAL = 2000         # participant bootstrap resamples (v28d line 843)
OWN_STOOL_SEED = 20260910          # stool-label permutation seed (v28d line 849)
OWN_STOOL_PERMUTATIONS = 99_999    # stool-label permutations (v28d line 849)
OG_PREVALENCES = [0.0, 0.1, 0.2]   # own-stool OG prevalence filters (v28d line 847)
Q_THRESHOLD = 0.05                 # adjusted q for significance (v28d line 889)
CROSSOMICS_SEED = 20260907         # native and functional permutation seed (v28d line 815; functional seed not stated at line 827)
CROSSOMICS_PERMUTATIONS = 1999     # whole-patient or whole-participant permutations (v28d lines 815, 827, 835)
EDS_SPLIT_SEEDS = [20260915, 20260916, 20260917]   # outer split seeds (v28d line 803)
EDS_LABEL_PERMUTATION_SEED = 20260915   # label permutations use the first split seed (not stated; analyze_score.py main())
EDS_BOOT_SEED = 20260918           # held-out prediction bootstrap seed (not stated; analyze_score.py main())
EDS_BOOT_DRAWS = 2000              # patient resamples of fixed held-out predictions (v28d line 805)
EDS_OUTER_FOLDS = 5                # shuffled patient-separated outer folds (v28d line 803)
EDS_INNER_FOLDS = 3                # unshuffled inner folds (v28d line 803)
EDS_PATIENTS = 85                  # fixed EDS roster (v28d line 803)
EDS_ALPHAS = [0.01, 0.1, 1, 10, 100]   # ridge penalty grid (v28d line 803; analyze_score.py ALPHAS)
ANTIBIOTIC_SIGNFLIP_DRAWS = 200_000    # Monte Carlo sign assignments (v28d line 795)
HISTIDINE_MC_DRAWS = 200_000       # sign assignments for n > 16 (v28d line 869)

WORDS = {1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six', 7: 'seven', 8: 'eight', 9: 'nine',
         10: 'ten', 11: 'eleven', 12: 'twelve', 13: 'thirteen', 14: 'fourteen', 15: 'fifteen', 16: 'sixteen'}


def word(n, capital=False):
    """Return the English word for a small count (as printed in the manuscript)."""
    w = WORDS[int(n)]
    return w.capitalize() if capital else w


def read(sd_root, unit, name, **kw):
    """Read ``source_data/<unit>/<name>`` as a DataFrame with round-trip float parsing."""
    return pd.read_csv(Path(sd_root) / unit / name, sep='\t', float_precision='round_trip', low_memory=False, **kw)


def bh(p):
    """Benjamini-Hochberg adjusted q values (scipy.stats.false_discovery_control)."""
    p = np.asarray(p, float)
    return false_discovery_control(p) if len(p) else p


def software():
    """Interpreter and package versions of the environment that runs this module."""
    import statsmodels
    import sklearn
    return (f"Python {platform.python_version()}; numpy {np.__version__}; scipy {scipy.__version__}; "
            f"pandas {pd.__version__}; statsmodels {statsmodels.__version__}; scikit-learn {sklearn.__version__}")


def family_check(d, keys, p='p', q='q'):
    """Recompute BH q within each family ``keys`` and compare with the producer's q.

    Parameters
    ----------
    d : DataFrame
        One row per test, with raw P values in ``p`` and the producer's q in ``q``.
    keys : list of str
        Columns that define one correction family.

    Returns
    -------
    DataFrame
        One row per family: its key values, ``size`` (tests in the family), ``tested``
        (tests with a P value), ``significant`` (q < 0.05) and ``max_abs_q_error``.
    """
    rows = []
    for k, g in d.groupby(keys, sort=True, dropna=False):
        ok = g[p].notna()
        err = float(np.max(np.abs(bh(g.loc[ok, p]) - g.loc[ok, q].to_numpy()))) if ok.any() else 0.0
        rows.append(dict(zip(keys, k if isinstance(k, tuple) else (k,)), size=len(g), tested=int(ok.sum()),
                         significant=int((g[q] < Q_THRESHOLD).sum()), max_abs_q_error=err))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------- SM8: CapScan regional
def boot_mean_ci(x, seed=BOOT_SEED_REGIONAL, draws=BOOT_DRAWS_REGIONAL):
    """Percentile bootstrap interval of a mean, resampling participants (producer boot())."""
    x = np.asarray(x)
    idx = np.random.default_rng(seed).integers(0, len(x), (draws, len(x)))
    m = x[idx].mean(1)
    return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def stat_regional_human_share(sd_root):
    """Human-assigned LFQ share, each capsule region versus stool (v28d line 843).

    Statistical unit: participant; n = 15 participants, paired. Exact two-sided sign-flip
    test over all 2^15 sign patterns, Benjamini-Hochberg over the four regional contrasts,
    95% percentile bootstrap of the mean difference (2,000 participant resamples).

    Returns
    -------
    dict
        ``tests`` (DataFrame, one row per contrast), ``max_abs_error`` against the producer's
        ``four_region_tests.tsv``, ``n``, ``patterns``, ``contrasts``, ``regions``.
    """
    cells = read(sd_root, 'Fig5', 'origin_fraction_cells.tsv')
    w = cells.pivot(index='participant', columns='region', values='Human-assigned')[REG] * 100
    n = len(w)
    signs = 2 * ((np.arange(2 ** n, dtype=np.uint64)[:, None] >> np.arange(n, dtype=np.uint64)) & 1).astype(np.int8) - 1
    rows = []
    for reg in REG[:4]:
        delta = (w[reg] - w.stool).to_numpy()
        p = float(np.mean(abs(signs @ delta) >= abs(delta.sum()) - 1e-12))
        lo, hi = boot_mean_ci(delta)
        rows.append(dict(region=reg, n=n, mean_region=float(w[reg].mean()), mean_stool=float(w.stool.mean()),
                         mean_difference=float(delta.mean()), ci_low=lo, ci_high=hi, p=p,
                         n_positive=int((delta > 0).sum()), statistic_sum=float(delta.sum())))
    t = pd.DataFrame(rows)
    t['q'] = bh(t.p)
    ref = read(sd_root, 'Fig5', 'four_region_tests.tsv').set_index('region').loc[t.region]
    cols = ['mean_region', 'mean_stool', 'mean_difference', 'ci_low', 'ci_high', 'p', 'q']
    err = float(np.max(np.abs(t.set_index('region')[cols].to_numpy() - ref[cols].to_numpy())))
    return dict(tests=t, max_abs_error=err, n=n, patterns=len(signs), contrasts=len(t), regions=len(REG))


def stat_own_stool(sd_root, unit):
    """Own-stool advantage of capsule OG profiles (v28d lines 847-849).

    Statistical unit: participant; n = 15. Mean over participants of (mean distance of the
    four capsule profiles to the 14 other stools minus distance to own stool). Two-sided
    stool-label permutation test (99,999 permutations, plus-one), 95% percentile bootstrap
    over participants (2,000 resamples). Computed from the 15 x 15 distance matrices.

    Returns
    -------
    dict
        ``tests`` (one row per prevalence filter) and ``max_abs_error`` against the
        producer's ``similarity_summary.tsv``.
    """
    dm = read(sd_root, unit, 'own_stool_distance_matrices.tsv')
    ref = read(sd_root, 'Fig5', 'similarity_summary.tsv').set_index('prevalence')
    rows = []
    for prev in OG_PREVALENCES:
        dd = dm[dm.prevalence == prev].pivot(index='capsule_participant', columns='stool_participant', values='mean_distance')
        assert list(dd.index) == list(dd.columns)
        dd = dd.to_numpy(); n = len(dd)
        perms = np.argsort(np.random.default_rng(OWN_STOOL_SEED).random((OWN_STOOL_PERMUTATIONS, n)), axis=1)
        own = np.diag(dd); other = (dd.sum(1) - own) / (n - 1)
        obs = float((other - own).mean())
        pdg = dd[np.arange(n)[None, :], perms].mean(1)
        null = (dd.sum(1).mean() - n * pdg) / (n - 1)
        exceed = int((abs(null) >= abs(obs) - 1e-14).sum())
        p = (1 + exceed) / (OWN_STOOL_PERMUTATIONS + 1)
        lo, hi = boot_mean_ci(other - own)
        rows.append(dict(prevalence=prev, participants=n, own_closer=int((own < other).sum()), mean_own=float(own.mean()),
                         mean_other=float(other.mean()), mean_advantage=obs, ci_low=lo, ci_high=hi, p=p,
                         null_exceedances=exceed, other_stools=n - 1))
    t = pd.DataFrame(rows)
    cols = ['participants', 'own_closer', 'mean_own', 'mean_other', 'mean_advantage', 'ci_low', 'ci_high', 'p']
    err = float(np.max(np.abs(t.set_index('prevalence')[cols].to_numpy(float) - ref.loc[OG_PREVALENCES, cols].to_numpy(float))))
    return dict(tests=t, max_abs_error=err)


# ----------------------------------------------------------------------------- SM8: cross-omics
def crossomics_tables(sd_root):
    """Concatenate PREDICT (SFig11) and CapScan 347 (Fig5) correlation results."""
    fun = pd.concat([read(sd_root, 'SFig11', 'functional_correlations_PREDICT.tsv'),
                     read(sd_root, 'Fig5', 'functional_correlations.tsv')], ignore_index=True)
    nat = pd.concat([read(sd_root, 'SFig11', 'native_correlations_PREDICT.tsv'),
                     read(sd_root, 'Fig5', 'native_correlations.tsv')], ignore_index=True)
    return fun, nat


def stat_crossomics(sd_root):
    """Functional and native profile concordance families (v28d lines 809-839, 893).

    Statistical unit: patient (PREDICT, n = 77) or participant (CapScan, 13 participants,
    65 participant-region profiles; sensitivity 15 participants, 75 profiles). Spearman
    correlation of distance-matrix upper triangles; two-sided permutation P (1,999
    permutations, plus-one); Benjamini-Hochberg per cohort, scope and prevalence
    (functional) or per cohort and feature set (native).

    Returns
    -------
    dict
        Counts, family tables and the CapScan classification-stability check.
    """
    fun, nat = crossomics_tables(sd_root)
    ff = family_check(fun, ['cohort', 'scope', 'prevalence'])
    nf = family_check(nat, ['cohort', 'analysis'])
    pre = fun[fun.cohort == 'PREDICT']; cap = fun[fun.cohort == 'CAPSCAN']
    mixed = cap[(cap.scope == 'displayed_mixed12_14') & (cap.prevalence == 0)]
    dna = mixed.layer_a.eq('metaG') | mixed.layer_b.eq('metaG')
    # classification at q < 0.05 across the four CapScan functional variants
    cap = cap.assign(sig=cap.q < Q_THRESHOLD)
    piv = cap.pivot_table(index=['level', 'layer_a', 'layer_b'], columns=['scope', 'prevalence'], values='sig', aggfunc='first')
    changed = int((piv.nunique(axis=1) > 1).sum())
    # permutations from the P-value grid: smallest attainable P is 1/(N+1)
    n_perm_fun = int(round(1 / fun.p.min())) - 1
    natp = nat.permutations.unique().tolist()
    return dict(
        functional_rows=len(fun), native_rows=len(nat), functional_families=ff, native_families=nf,
        predict_functional_family=int(ff[ff.cohort == 'PREDICT']['size'].unique()[0]),
        capscan_functional_family=int(ff[ff.cohort == 'CAPSCAN']['size'].unique()[0]),
        predict_native_family=int(nf[nf.cohort == 'PREDICT']['size'].unique()[0]),
        capscan_native_family=int(nf[nf.cohort == 'CAPSCAN']['size'].unique()[0]),
        max_abs_q_error=float(max(ff.max_abs_q_error.max(), nf.max_abs_q_error.max())),
        predict_patients=sorted(pre.n_donors.unique().tolist()),
        predict_native_profiles=sorted(nat[nat.cohort == 'PREDICT'].n_profiles.unique().tolist()),
        capscan_matched_participants=sorted(cap[cap.scope == 'matched12'].n_donors.unique().tolist()),
        capscan_matched_profiles=sorted(cap[cap.scope == 'matched12'].n_units.unique().tolist()),
        capscan_native_participants=sorted(nat[nat.cohort == 'CAPSCAN'].n_participants.unique().tolist()),
        mixed_nondna_tests=int((~dna).sum()), mixed_nondna_participants=sorted(mixed[~dna].n_donors.unique().tolist()),
        mixed_nondna_profiles=sorted(mixed[~dna].n_units.unique().tolist()),
        mixed_dna_tests=int(dna.sum()), mixed_dna_participants=sorted(mixed[dna].n_donors.unique().tolist()),
        prevalence_filter=float(fun.prevalence.max()), classification_changes=changed,
        functional_permutations=n_perm_fun, native_permutations=natp,
        functional_min_p=float(fun.p.min()))


def stat_controlled(sd_root):
    """Controlled distance sensitivity: 176 estimates in 24 BH families (v28d lines 819-821, 893).

    Statistical unit: patient (PREDICT, n = 77) or participant (CapScan, n = 13, 65
    profiles). Each feature/distance/adjustment/outer-correlation variant is its own family
    of ten (PREDICT) or six (CapScan) comparisons.
    """
    d = read(sd_root, 'SFig11', 'controlled_176_comparisons.tsv')
    f = family_check(d, ['cohort', 'features', 'distance', 'adjustment', 'outer_correlation'])
    return dict(rows=len(d), families=f, family_sizes=sorted(f['size'].unique().tolist()),
                predict_family=int(f[f.cohort == 'PREDICT']['size'].unique()[0]),
                capscan_family=int(f[f.cohort == 'CAPSCAN']['size'].unique()[0]),
                permutations=int(round(1 / d.p.min())) - 1, max_abs_q_error=float(f.max_abs_q_error.max()),
                capscan_participants=sorted(d[d.cohort == 'CAPSCAN'].n_participants.unique().tolist()))


def stat_individual_functions(sd_root):
    """Individual-function Spearman correlations: nine BH families (v28d lines 829, 893).

    Statistical unit: patient; n per test = patients positive in both layers (at least 35
    of 77). scipy.stats.spearmanr two-sided P (t approximation).
    """
    d = read(sd_root, 'STab5', 'crossomics_same_function_correlations.tsv')
    f = family_check(d, ['level', 'layer_a', 'layer_b'])
    return dict(rows=len(d), families=f, n_families=len(f), min_patients=int(d.n_patients.min()),
                max_abs_q_error=float(f.max_abs_q_error.max()))


def stat_eds_screen(sd_root):
    """EDS association screen: BH per layer, normalization and model (v28d line 799)."""
    d = read(sd_root, 'SFig10', 'all_EDS_associations.tsv')
    f = family_check(d, ['layer', 'normalization', 'model'])
    tested = d[d.p.notna()]
    return dict(rows=len(d), families=f, min_n_tested=int(tested.n.min()), max_abs_q_error=float(f.max_abs_q_error.max()))


# ----------------------------------------------------------------------------- SM8: antibiotics
def stat_antibiotic(sd_root):
    """PREDICT antibiotic endpoints: two-test BH family per exposure definition (v28d lines 793-795).

    Statistical unit: patient (paired on/off). The sign-flip test and bootstrap are
    recomputed by ``tables/STab3`` from controlled inputs; here the family structure and the
    Monte Carlo grid of the P values are checked.
    """
    d = read(sd_root, 'STab3', 'endpoints_by_definition.tsv')
    f = family_check(d, ['scenario'], p='p_signflip_MC', q='q_BH_two_endpoints_within_definition')
    k = d.p_signflip_MC * (ANTIBIOTIC_SIGNFLIP_DRAWS + 1)
    on_grid = bool(np.all(np.abs(k - np.round(k)) < 1e-6))
    return dict(table=d, families=f, family_size=int(f['size'].unique()[0]), on_grid=on_grid,
                n_patients=dict(zip(d.scenario, d.n_patients)), max_abs_q_error=float(f.max_abs_q_error.max()))


# ----------------------------------------------------------------------------- SM8: held-out EDS prediction
def stat_eds_prediction(sd_root, unit):
    """Exploratory held-out EDS prediction (v28d lines 803-805).

    Statistical unit: patient; n = 85. For the first split seed, two-sided Spearman
    label-permutation P (199 permutations, plus-one) per grouping x model, BH over the
    eight tests; 95% percentile intervals from 2,000 patient resamples of fixed held-out
    predictions for Spearman rho and for baseline-minus-model mean absolute error.
    """
    m = read(sd_root, unit, 'eds_prediction_model_metrics.tsv')
    perm = read(sd_root, unit, 'eds_prediction_permutations.tsv')
    inc = read(sd_root, unit, 'eds_prediction_increment.tsv')
    h = read(sd_root, unit, 'eds_prediction_heldout.tsv')
    folds = read(sd_root, unit, 'eds_prediction_outer_folds.tsv')
    audit = read(sd_root, unit, 'eds_prediction_fit_audit.tsv')
    alphas = read(sd_root, unit, 'eds_prediction_selected_alphas.tsv')
    first = EDS_SPLIT_SEEDS[0]
    y = h[(h.seed == first)].groupby('row').dysbiosis_score.first().sort_index().to_numpy()
    n = len(y)
    # held-out metrics recomputed from the predictions
    rho_err = 0.0
    for (g, s, mo), b in h.groupby(['grouping', 'seed', 'model']):
        b = b.sort_values('row')
        r = spearmanr(b.dysbiosis_score, b.score).statistic
        ref = m[(m.grouping == g) & (m.seed == s) & (m.model == mo)].rho.iloc[0]
        rho_err = max(rho_err, abs(r - ref))
    # permutation P and BH over the eight first-split tests
    tests = m[(m.seed == first) & m.permutation_p.notna()].copy()
    n_perm = perm.groupby(['grouping', 'model']).size().unique().tolist()
    p_new = []
    for _, r in tests.iterrows():
        nr = perm[(perm.grouping == r.grouping) & (perm.model == r.model)].rho.to_numpy()
        p_new.append((1 + (abs(nr) >= abs(r.rho) - 1e-12).sum()) / (1 + len(nr)))
    tests['p_recomputed'] = p_new
    tests['q_recomputed'] = bh(tests.p_recomputed)
    perm_err = float(np.max(np.abs(tests.p_recomputed - tests.permutation_p)))
    q_err = float(np.max(np.abs(tests.q_recomputed - tests.permutation_q_8)))
    # bootstrap intervals
    boots = np.random.default_rng(EDS_BOOT_SEED).integers(0, n, (EDS_BOOT_DRAWS, n))
    ci_err = 0.0; inc_err = 0.0; ci_rows = []
    for g in tests.grouping.unique():
        pred = {mo: h[(h.grouping == g) & (h.seed == first) & (h.model == mo)].sort_values('row').score.to_numpy()
                for mo in ['protein', 'baseline', 'combined', 'distance']}
        for mo, val in pred.items():
            rs = np.array([spearmanr(y[ix], val[ix]).statistic for ix in boots]); lo, hi = np.quantile(rs, [.025, .975])
            ref = m[(m.grouping == g) & (m.seed == first) & (m.model == mo)].iloc[0]
            ci_err = max(ci_err, abs(lo - ref.rho_ci_low), abs(hi - ref.rho_ci_high))
            ci_rows.append(dict(grouping=g, model=mo, rho_ci_low=lo, rho_ci_high=hi))
            if mo in ['protein', 'combined']:
                delta = abs(pred['baseline'] - y) - abs(val - y)
                ci = np.quantile(delta[boots].mean(1), [.025, .975])
                r_ = inc[(inc.grouping == g) & (inc.model == mo)].iloc[0]
                inc_err = max(inc_err, abs(delta.mean() - r_.MAE_baseline_minus_model), abs(ci[0] - r_.ci_low), abs(ci[1] - r_.ci_high))
    # folds: KFold(5, shuffle, seed) reproduces the stored outer folds; inner folds per outer fold
    fold_ok = True
    for s in EDS_SPLIT_SEEDS:
        for f, (_, te) in enumerate(KFold(EDS_OUTER_FOLDS, shuffle=True, random_state=s).split(np.arange(n))):
            stored = sorted(folds[(folds.seed == s) & (folds.fold == f)].row.tolist())
            fold_ok &= stored == sorted(te.tolist())
    # the fit audit logs both protein-grouping controls without a grouping column
    groupings = m.grouping.nunique()
    inner = (audit[audit.stage == 'inner'].groupby(['seed', 'model', 'outer_fold']).size() // groupings).unique().tolist()
    return dict(n=n, tests=tests, n_tests=len(tests), permutations=n_perm, perm_err=perm_err, q_err=q_err,
                rho_err=float(rho_err), ci_err=float(ci_err), inc_err=float(inc_err), ci=pd.DataFrame(ci_rows),
                fold_ok=bool(fold_ok), outer_folds=int(folds.fold.nunique()), inner_folds=inner,
                seeds=sorted(m.seed.unique().tolist()), alphas_selected=sorted(alphas.alpha.unique().tolist()),
                increment=inc)


# ----------------------------------------------------------------------------- SM8: counts
def stat_dna_runs(sd_root, unit):
    """CapScan metagenome runs: functional aggregate (301) and native profile (245)."""
    d = read(sd_root, unit, 'capscan_dna_runs.tsv')
    nat = d[d.native_profile]
    cov = nat.groupby('participant').region.nunique()
    return dict(functional=int(d.functional_aggregate.sum()), native=int(len(nat)),
                native_participants=int(len(cov)), native_complete=bool((cov == 5).all()))


def stat_native_features(sd_root, unit):
    """PREDICT native DNA and RNA sequence features with positive abundance in 77 profiles."""
    d = read(sd_root, unit, 'predict_native_feature_summary.tsv').set_index(['analysis', 'layer'])
    return dict(dna=int(d.loc[('sensitivity_all_features', 'DNA'), 'retained_features']),
                rna=int(d.loc[('sensitivity_all_features', 'RNA'), 'retained_features']),
                dna_median_retained=float(d.loc[('primary_prevalence10', 'DNA'), 'median_signal_fraction']),
                rna_median_retained=float(d.loc[('primary_prevalence10', 'RNA'), 'median_signal_fraction']))


def stat_nonspecific(sd_root, unit):
    """Non-specific Sage searches: completed searches and settings (v28d line 855)."""
    d = read(sd_root, unit, 'nonspecific_search_summary.tsv', keep_default_na=False).set_index('cohort')
    return d.to_dict(orient='index')


def stat_hstool_grouping(sd_root, unit):
    """Healthy-stool grouping receipt counts (v28d line 909)."""
    d = read(sd_root, unit, 'hstool_grouping_receipt.tsv')
    return dict(zip(d.key, d.value.astype(int)))


def stat_contributors(sd_root):
    """CapScan leading-species eligibility and the Fig. 5e persistence (v28d lines 785-787)."""
    c = read(sd_root, 'Fig5', 'contributor_cells.tsv')
    per_ko = c.groupby('KO').participant.nunique()
    mr = read(sd_root, 'Fig5', 'mosaic_donor_regions.tsv', keep_default_na=False)
    res = mr[mr.status == 'Resolved']
    evaluated = persistent = 0
    for (ko, pt), g in mr.groupby(['KO', 'participant']):
        r = g[g.status == 'Resolved'].set_index('region').species
        if 'stool' not in r.index:
            continue
        caps = r.drop('stool')
        if len(caps) == 0:
            continue
        evaluated += 1
        persistent += int((caps == r['stool']).all())
    return dict(cells=len(c), functions=int(c.KO.nunique()), participants=int(c.participant.nunique()),
                min_participants_per_function=int(per_ko.min()), mosaic_functions=int(mr.KO.nunique()),
                mosaic_participants=int(mr.participant.nunique()), evaluated=evaluated, persistent=persistent,
                persistence_pct=100 * persistent / evaluated, resolved_rows=len(res))


def stat_histidine_families(sd_root):
    """Histidine and terminal-fraction BH families (v28d lines 869, 873, 893)."""
    t = read(sd_root, 'Fig6', 'capscan_terminal_statistical_families.tsv')
    h = read(sd_root, 'Fig6', 'capscan_histidine_statistical_families.tsv')
    p = read(sd_root, 'Fig6', 'predict_histidine_statistical_families.tsv')
    ft = family_check(t, ['family']).assign(cohort='CAPSCAN', table='terminal')
    fh = family_check(h, ['family']).assign(cohort='CAPSCAN', table='histidine')
    fp = family_check(p, ['family']).assign(cohort='PREDICT', table='histidine')
    fam = pd.concat([ft, fh, fp], ignore_index=True)
    sizes = {(c, f): int(s) for c, f, s in zip(fam.cohort, fam.family, fam['size'])}
    paired = t[t.family == 'terminal_fraction__length_9']
    # exact enumeration for small n: P is a multiple of 2^-n; Monte Carlo for PREDICT n = 93
    cap_exact = h[(h.endpoint == 'source_enrichment') & h.p.notna() & (h.n <= 16)]
    grid_exact = bool(np.all(np.abs(cap_exact.p * 2.0 ** cap_exact.n - np.round(cap_exact.p * 2.0 ** cap_exact.n)) < 1e-6))
    pe = p[(p.endpoint == 'source_enrichment') & p.p.notna()]
    small, large = pe[pe.n <= 16], pe[pe.n > 16]
    ks = small.p * 2.0 ** small.n
    k = large.p * (HISTIDINE_MC_DRAWS + 1)
    grid_exact &= bool(np.all(np.abs(ks - np.round(ks)) < 1e-6))
    grid_mc = bool(len(large) > 0 and np.all(np.abs(k - np.round(k)) < 1e-6))
    return dict(families=fam, sizes=sizes, max_abs_q_error=float(fam.max_abs_q_error.max()),
                terminal_paired_participants=sorted(paired.n_participants.dropna().astype(int).unique().tolist()),
                predict_patients=int(p.n.max()),
                capscan_exact_grid=grid_exact, predict_mc_grid=grid_mc,
                terminal_min_length=int(t[t.family == 'terminal_fraction__length_9'].min_length.unique()[0]))


def stat_strata(sd_root):
    """Measured-spectrum strata shipped in Source Data (v28d line 877)."""
    s = read(sd_root, 'SFig14', 'spectra/all_strata.tsv')
    return dict(predict=int((s.cohort == 'PREDICT').sum()), capscan=int((s.cohort == 'CAPSCAN').sum()))


def stat_partition(sd_root):
    """Descriptive participant/sampling-group variance partition (Fig5 producer table)."""
    d = read(sd_root, 'Fig5', 'orthogonal_partition.tsv')
    return d[d.prevalence == 0]


def stat_evenness(sd_root):
    """Paired evenness and leading-species share tests as computed by the producer (Fig5 table)."""
    return read(sd_root, 'Fig5', 'evenness_summary_tests.tsv')


# ----------------------------------------------------------------------------- collect
def collect(sd_root, unit):
    """Run every computation and return (R, details).

    ``R`` maps NUMBER_SPEC keys to recomputed values; ``details`` keeps the full results
    for STATS.tsv and for display.
    """
    D = dict(regional=stat_regional_human_share(sd_root), own=stat_own_stool(sd_root, unit),
             cross=stat_crossomics(sd_root), controlled=stat_controlled(sd_root),
             individual=stat_individual_functions(sd_root), screen=stat_eds_screen(sd_root),
             abx=stat_antibiotic(sd_root), eds=stat_eds_prediction(sd_root, unit),
             dna=stat_dna_runs(sd_root, unit), native=stat_native_features(sd_root, unit),
             nonspec=stat_nonspecific(sd_root, unit), hstool=stat_hstool_grouping(sd_root, unit),
             contrib=stat_contributors(sd_root), his=stat_histidine_families(sd_root), strata=stat_strata(sd_root),
             partition=stat_partition(sd_root), evenness=stat_evenness(sd_root), _sd_root=sd_root)
    # every reproduced producer value must agree to floating-point precision
    tol = 1e-9
    checks = {'regional': D['regional']['max_abs_error'], 'own_stool': D['own']['max_abs_error'],
              'crossomics_q': D['cross']['max_abs_q_error'], 'controlled_q': D['controlled']['max_abs_q_error'],
              'individual_q': D['individual']['max_abs_q_error'], 'screen_q': D['screen']['max_abs_q_error'],
              'antibiotic_q': D['abx']['max_abs_q_error'], 'eds_perm_p': D['eds']['perm_err'], 'eds_q': D['eds']['q_err'],
              'eds_rho': D['eds']['rho_err'], 'eds_rho_ci': D['eds']['ci_err'], 'eds_increment': D['eds']['inc_err'],
              'histidine_q': D['his']['max_abs_q_error']}
    D['checks'] = checks
    assert all(v < tol for v in checks.values()), checks
    assert D['eds']['fold_ok'] and D['abx']['on_grid'] and D['his']['capscan_exact_grid'] and D['his']['predict_mc_grid']
    rg, ow, cr, co, ed, his = D['regional'], D['own'], D['cross'], D['controlled'], D['eds'], D['his']
    ns = D['nonspec']; c = D['contrib']; sz = his['sizes']
    one = lambda xs: xs[0] if len(xs) == 1 else xs   # noqa: E731
    R = {
        # SM8 contributors
        'L785_twelve': word(c['min_participants_per_function']),
        'L785_fifteen': word(c['participants']),
        'L787_24a': c['mosaic_functions'], 'L787_fifteen': word(c['mosaic_participants']),
        'L787_persistence': c['persistence_pct'], 'L787_24b': c['mosaic_functions'],
        # antibiotics
        'L795_42': D['abx']['n_patients']['A_current'],
        'L795_200000': ANTIBIOTIC_SIGNFLIP_DRAWS if D['abx']['on_grid'] else None,
        'L795_two': word(D['abx']['family_size']),
        # held-out EDS prediction
        'L803_85': ed['n'], 'L803_three_inner': word(one(ed['inner_folds'])), 'L803_five_outer': word(ed['outer_folds']),
        'L803_three_seeds': word(len(ed['seeds']), capital=True), 'L803_seed_range': f"{ed['seeds'][0]}–{ed['seeds'][-1]}",
        'L803_199': one(ed['permutations']), 'L805_eight': word(ed['n_tests']),
        'L805_2000': 'Two thousand' if ed['ci_err'] < 1e-9 and ed['inc_err'] < 1e-9 else None,
        # native and functional cross-omics
        'L809_77': one(cr['predict_native_profiles']), 'L809_thirteen': word(one(cr['capscan_native_participants'])),
        'L809_five': word(rg['regions']), 'L809_245': D['dna']['native'],
        'L815_thirteen': word(one(cr['capscan_native_participants'])), 'L815_1999': one(cr['native_permutations']),
        'L815_ten': word(cr['predict_native_family']), 'L815_six': word(cr['capscan_native_family']),
        'L821_1999': co['permutations'], 'L821_ten': word(co['predict_family']), 'L821_six': word(co['capscan_family']),
        'L821_176': co['rows'],
        'L827_77a': one(cr['predict_patients']), 'L827_1999': cr['functional_permutations'],
        'L827_28': cr['predict_functional_family'], 'L827_10pct': 100 * cr['prevalence_filter'],
        'L827_77b': one(cr['predict_native_profiles']), 'L827_rna': D['native']['rna'], 'L827_dna': D['native']['dna'],
        'L833_sixteen': word(cr['capscan_functional_family']), 'L833_thirteen': word(one(cr['capscan_matched_participants'])),
        'L833_65': one(cr['capscan_matched_profiles']), 'L833_301': D['dna']['functional'],
        'L835_four': word(rg['regions'] - 1),
        'L837_10pct': 100 * cr['prevalence_filter'], 'L837_seven': word(cr['mixed_nondna_tests']),
        'L837_fifteen': word(one(cr['mixed_nondna_participants'])), 'L837_75': one(cr['mixed_nondna_profiles']),
        'L837_nine': word(cr['mixed_dna_tests']), 'L837_thirteen': word(one(cr['mixed_dna_participants'])),
        'L837_no': 'no' if cr['classification_changes'] == 0 else str(cr['classification_changes']),
        'L839_120': cr['functional_rows'], 'L839_32': cr['native_rows'], 'L839_176': co['rows'],
        # CapScan regional and own-stool
        'L843_Fifteen': word(rg['n'], capital=True), 'L843_five': word(rg['regions']), 'L843_2000': BOOT_DRAWS_REGIONAL if rg['max_abs_error'] < 1e-9 else None,
        'L843_32768': rg['patterns'], 'L843_fifteen': word(rg['n']), 'L843_four': word(rg['contrasts']),
        'L843_fourteen': word(one(ow['tests'].other_stools.unique().tolist())),
        'L847_75': ow['tests'].participants.iloc[0] * len(REG), 'L847_fifteen': word(ow['tests'].participants.iloc[0]),
        'L847_p': float(ow['tests'].p.max()),
        'L847_four': word(len(REG) - 1), 'L847_fourteen': word(one(ow['tests'].other_stools.unique().tolist())),
        'L849_99999': int(round(1 / ow['tests'].p.max())) - 1 if (ow['tests'].null_exceedances == 0).all() else None,
        'L849_15': int(ow['tests'].participants.iloc[0]), 'L849_75': int(ow['tests'].participants.iloc[0]) * len(REG),
        # SM9 non-specific searches
        'L855_capscan': int(ns['CAPSCAN']['completed_searches']), 'L855_predict': int(ns['PREDICT']['completed_searches']),
        'L855_lengths': f"{ns['CAPSCAN']['min_len']}–{ns['CAPSCAN']['max_len']}",
        'L855_mass': f"{int(float(ns['CAPSCAN']['peptide_min_mass']))}–{int(float(ns['CAPSCAN']['peptide_max_mass'])):,}",
        'L855_tol': float(ns['CAPSCAN']['precursor_tol_ppm']) if ns['CAPSCAN']['precursor_tol_ppm'] == ns['CAPSCAN']['fragment_tol_ppm'] == ns['PREDICT']['precursor_tol_ppm'] else None,
        'L855_two': word(int(ns['CAPSCAN']['max_variable_mods'])),
        'L855_terminus': f"{his['terminal_min_length']}–{ns['CAPSCAN']['max_len']}",
        'L869_93': his['predict_patients'],
        'L873_32400': sz[('CAPSCAN', 'all_motif_positions_residues')],
        'L877_four': word(D['strata']['predict']),
        # SM10 testing-family map
        'L893_two': word(D['abx']['family_size']), 'L893_four': word(rg['contrasts']),
        'L893_28': cr['predict_functional_family'], 'L893_16': cr['capscan_functional_family'],
        'L893_ten': word(cr['predict_native_family']), 'L893_six': word(cr['capscan_native_family']),
        'L893_nine': word(D['individual']['n_families']), 'L893_seven': word(cr['mixed_nondna_tests']),
        'L893_fifteen': word(one(cr['mixed_nondna_participants'])),
        'L893_4536': sz[('CAPSCAN', 'terminal_fraction__length_9')], 'L893_13608': sz[('CAPSCAN', 'all18_nonKR_composition__length_9')],
        'L893_9000': sz[('CAPSCAN', 'source_P1_all_residues')], 'L893_1560': sz[('PREDICT', 'source_P1_all_residues')],
        'L893_90': sz[('CAPSCAN', 'all_reciprocal_controls')], 'L893_12': sz[('PREDICT', 'all_reciprocal_controls')],
        'L893_2520': sz[('CAPSCAN', 'H_spatial_all_routes')], 'L893_32400': sz[('CAPSCAN', 'all_motif_positions_residues')],
        # SM11
        'L909_58': D['hstool']['acquisitions'], 'L909_59428': D['hstool']['quantified_protein_groups_direct_peptide_assignment'],
        'L909_1909973': D['hstool']['eligible_lfq_features'],
    }
    return R, D


# ----------------------------------------------------------------------------- NUMBERS.tsv
U = 'text_statistics/supplementary_methods_8_11/stats.py'
SDU = 'source_data/supplementary_methods_8_11'
#: (v28d_line, quote, printed_value, key or None, code function, source_data, status)
#: status None = MATCH when recomputed agrees; 'v28f_change:<x>' or 'OPEN' given explicitly.
NUMBER_SPEC = [
    ('783', 'in at least 80% of acquisitions in every experimental design', '80%', None, '', '', 'OPEN'),
    ('783', 'requiring resolved leaders in at least 80% of acquisitions in each of the two compared arms', '80%', None, '', '', 'OPEN'),
    ('785', 'Eligibility required a resolved leader in at least twelve of the fifteen participants', 'twelve', 'L785_twelve', 'stat_contributors', 'source_data/Fig5/contributor_cells.tsv', None),
    ('785', 'Eligibility required a resolved leader in at least twelve of the fifteen participants', 'fifteen', 'L785_fifteen', 'stat_contributors', 'source_data/Fig5/contributor_cells.tsv', None),
    ('785', 'A stricter subset required at least 20% species-resolved signal', '20%', None, '', '', 'OPEN'),
    ('785', 'a leading-species margin of at least five percentage points', 'five', None, '', '', 'OPEN'),
    ('787', 'The Fig. 5e mosaic displays 24 functions', '24', 'L787_24a', 'stat_contributors', 'source_data/Fig5/mosaic_donor_regions.tsv', None),
    ('787', 'Among these 24 displayed functions and the fifteen CapScan participants', 'fifteen', 'L787_fifteen', 'stat_contributors', 'source_data/Fig5/mosaic_donor_regions.tsv', None),
    ('787', 'This 30.8% therefore refers only to the 24 displayed functions', '30.8%', 'L787_persistence', 'stat_contributors', 'source_data/Fig5/mosaic_donor_regions.tsv', None),
    ('787', 'This 30.8% therefore refers only to the 24 displayed functions', '24', 'L787_24b', 'stat_contributors', 'source_data/Fig5/mosaic_donor_regions.tsv', None),
    ('791', 'proteins pooled from all 296 available PREDICT assembly FASTAs', '296', None, '', '', 'OPEN'),
    ('793', 'all nine core antibiotic fields', 'nine', None, '', '', 'OPEN'),
    ('795', 'We report mean paired changes with 95% percentile bootstrap intervals', '95%', None, '', '', 'OPEN'),
    ('795', 'from 50,000 resamples of the 42 patients (seed 20260907)', '50,000', None, '', '', 'OPEN'),
    ('795', 'from 50,000 resamples of the 42 patients (seed 20260907)', '42', 'L795_42', 'stat_antibiotic', 'source_data/STab3/endpoints_by_definition.tsv', None),
    ('795', 'from 50,000 resamples of the 42 patients (seed 20260907)', '20260907', None, '', '', 'OPEN'),
    ('795', 'using 200,000 random patient-level sign assignments with a plus-one correction', '200,000', 'L795_200000', 'stat_antibiotic', 'source_data/STab3/endpoints_by_definition.tsv', None),
    ('795', 'formed one two-test Benjamini–Hochberg [11] family', 'two', 'L795_two', 'stat_antibiotic', 'source_data/STab3/endpoints_by_definition.tsv', None),
    ('803', 'Exploratory protein-based EDS prediction used the fixed 85-patient roster', '85', 'L803_85', 'stat_eds_prediction', f'{SDU}/eds_prediction_heldout.tsv', None),
    ('803', 'OG features with positive intensity in at least 10% of training observations', '10%', None, '', '', 'OPEN'),
    ('803', 'Ridge penalties 0.01, 0.1, 1, 10 and 100 were selected', '0.01, 0.1, 1, 10 and 100', None, '', '', 'OPEN'),
    ('803', 'by mean squared error in three inner folds', 'three', 'L803_three_inner', 'stat_eds_prediction', f'{SDU}/eds_prediction_fit_audit.tsv', None),
    ('803', 'within each of five shuffled patient-separated outer folds', 'five', 'L803_five_outer', 'stat_eds_prediction', f'{SDU}/eds_prediction_outer_folds.tsv', None),
    ('803', 'Three outer split seeds (20260915–20260917) were retained', 'Three', 'L803_three_seeds', 'stat_eds_prediction', f'{SDU}/eds_prediction_model_metrics.tsv', None),
    ('803', 'Three outer split seeds (20260915–20260917) were retained', '20260915–20260917', 'L803_seed_range', 'stat_eds_prediction', f'{SDU}/eds_prediction_model_metrics.tsv', None),
    ('803', 'A separate distance score averaged the five nearest training profiles', 'five', None, '', '', 'OPEN'),
    ('803', 'For the first outer split, 199 label permutations', '199', 'L803_199', 'stat_eds_prediction', f'{SDU}/eds_prediction_permutations.tsv', None),
    ('805', 'Benjamini–Hochberg adjustment across eight grouping-by-model tests', 'eight', 'L805_eight', 'stat_eds_prediction', f'{SDU}/eds_prediction_model_metrics.tsv', None),
    ('805', 'Two thousand patient resamples of fixed held-out predictions', 'Two thousand', 'L805_2000', 'stat_eds_prediction', f'{SDU}/eds_prediction_heldout.tsv', None),
    ('809', 'one earliest complete-stack stool specimen per patient (77 patients)', '77', 'L809_77', 'stat_crossomics', 'source_data/SFig11/native_correlations_PREDICT.tsv', None),
    ('809', 'CapScan retained thirteen participants with all five regions, using 245 DNA runs', 'thirteen', 'L809_thirteen', 'stat_crossomics', 'source_data/Fig5/native_correlations.tsv', None),
    ('809', 'CapScan retained thirteen participants with all five regions, using 245 DNA runs', 'five', 'L809_five', 'stat_regional_human_share', 'source_data/Fig5/origin_fraction_cells.tsv', None),
    ('809', 'CapScan retained thirteen participants with all five regions, using 245 DNA runs', '245', 'L809_245', 'stat_dna_runs', f'{SDU}/capscan_dna_runs.tsv', None),
    ('813', 'restricted gene and protein features to those present in at least 10% of profiles', '10%', None, '', '', 'OPEN'),
    ('815', 'Metabolites retained features finite in at least 70% of profiles', '70%', None, '', '', 'OPEN'),
    ('815', 'average region effects within the thirteen-participant complete population', 'thirteen', 'L815_thirteen', 'stat_crossomics', 'source_data/Fig5/native_correlations.tsv', None),
    ('815', 'Tests used 1,999 permutations (seed 20260907)', '1,999', 'L815_1999', 'stat_crossomics', 'source_data/SFig11/native_correlations_PREDICT.tsv; source_data/Fig5/native_correlations.tsv', None),
    ('815', 'Tests used 1,999 permutations (seed 20260907)', '20260907', None, '', '', 'OPEN'),
    ('815', 'Benjamini–Hochberg families contained ten or six comparisons', 'ten', 'L815_ten', 'stat_crossomics', 'source_data/SFig11/native_correlations_PREDICT.tsv', None),
    ('815', 'Benjamini–Hochberg families contained ten or six comparisons', 'six', 'L815_six', 'stat_crossomics', 'source_data/Fig5/native_correlations.tsv', None),
    ('819', 'Both all-observed and ≥10% prevalence inputs were retained', '10%', None, '', '', 'OPEN'),
    ('821', 'used the same design-preserving 1,999 permutations', '1,999', 'L821_1999', 'stat_controlled', 'source_data/SFig11/controlled_176_comparisons.tsv', None),
    ('821', 'its own ten- or six-comparison Benjamini–Hochberg family', 'ten', 'L821_ten', 'stat_controlled', 'source_data/SFig11/controlled_176_comparisons.tsv', None),
    ('821', 'its own ten- or six-comparison Benjamini–Hochberg family', 'six', 'L821_six', 'stat_controlled', 'source_data/SFig11/controlled_176_comparisons.tsv', None),
    ('821', 'All 176 estimates from these controlled distance comparisons were retained', '176', 'L821_176', 'stat_controlled', 'source_data/SFig11/controlled_176_comparisons.tsv', None),
    ('827', 'one earliest complete molecular observation was retained for each of 77 patients', '77', 'L827_77a', 'stat_crossomics', 'source_data/SFig11/functional_correlations_PREDICT.tsv', None),
    ('827', 'required finite measurements in at least 70% of profiles', '70%', None, '', '', 'OPEN'),
    ('827', 'Two-sided significance was assessed using 1,999 patient-label permutations', '1,999', 'L827_1999', 'stat_crossomics', 'source_data/SFig11/functional_correlations_PREDICT.tsv', None),
    ('827', 'Benjamini–Hochberg correction was applied across the 28 PREDICT profile-level comparisons', '28', 'L827_28', 'stat_crossomics', 'source_data/SFig11/functional_correlations_PREDICT.tsv', None),
    ('827', 'restricting DNA, RNA and protein features to those positive in at least 10% of profiles', '10%', 'L827_10pct', 'stat_crossomics', 'source_data/SFig11/functional_correlations_PREDICT.tsv', None),
    ('827', 'Across the same 77 native profiles', '77', 'L827_77b', 'stat_crossomics', 'source_data/SFig11/native_correlations_PREDICT.tsv', None),
    ('827', '5,391,526 encoded-sequence features had positive RNA abundance', '5,391,526', 'L827_rna', 'stat_native_features', f'{SDU}/predict_native_feature_summary.tsv', None),
    ('827', 'compared with 14,822,078 for DNA', '14,822,078', 'L827_dna', 'stat_native_features', f'{SDU}/predict_native_feature_summary.tsv', None),
    ('833', 'All sixteen comparisons among microbial DNA', 'sixteen', 'L833_sixteen', 'stat_crossomics', 'source_data/Fig5/functional_correlations.tsv', None),
    ('833', 'used the same thirteen participants with complete DNA profiles (65 profiles)', 'thirteen', 'L833_thirteen', 'stat_crossomics', 'source_data/Fig5/functional_correlations.tsv', None),
    ('833', 'used the same thirteen participants with complete DNA profiles (65 profiles)', '65', 'L833_65', 'stat_crossomics', 'source_data/Fig5/functional_correlations.tsv', None),
    ('833', 'gene TPM values from 301 completed metagenome runs', '301', 'L833_301', 'stat_dna_runs', f'{SDU}/capscan_dna_runs.tsv', None),
    ('833', 'Metabolites required at least 70% finite observations', '70%', None, '', '', 'OPEN'),
    ('835', 'ordinary least squares with an intercept and four region indicators', 'four', 'L835_four', 'stat_regional_human_share', 'source_data/Fig5/origin_fraction_cells.tsv', None),
    ('837', 'the same comparisons with a ≥10% feature-prevalence filter', '10%', 'L837_10pct', 'stat_crossomics', 'source_data/Fig5/functional_correlations.tsv', None),
    ('837', 'seven non-DNA comparisons used fifteen participants (75 profiles)', 'seven', 'L837_seven', 'stat_crossomics', 'source_data/Fig5/functional_correlations.tsv', None),
    ('837', 'seven non-DNA comparisons used fifteen participants (75 profiles)', 'fifteen', 'L837_fifteen', 'stat_crossomics', 'source_data/Fig5/functional_correlations.tsv', None),
    ('837', 'seven non-DNA comparisons used fifteen participants (75 profiles)', '75', 'L837_75', 'stat_crossomics', 'source_data/Fig5/functional_correlations.tsv', None),
    ('837', 'nine DNA comparisons used thirteen', 'nine', 'L837_nine', 'stat_crossomics', 'source_data/Fig5/functional_correlations.tsv', None),
    ('837', 'nine DNA comparisons used thirteen', 'thirteen', 'L837_thirteen', 'stat_crossomics', 'source_data/Fig5/functional_correlations.tsv', None),
    ('837', 'no comparison changed classification at adjusted q<0.05', 'no', 'L837_no', 'stat_crossomics', 'source_data/Fig5/functional_correlations.tsv', None),
    ('839', 'The Source Data contain the 120 functional and 32 native correlation results', '120', 'L839_120', 'stat_crossomics', 'source_data/SFig11/functional_correlations_PREDICT.tsv; source_data/Fig5/functional_correlations.tsv', None),
    ('839', 'The Source Data contain the 120 functional and 32 native correlation results', '32', 'L839_32', 'stat_crossomics', 'source_data/SFig11/native_correlations_PREDICT.tsv; source_data/Fig5/native_correlations.tsv', None),
    ('839', 'all plotted cross-omics values and the 176 controlled estimates', '176', 'L839_176', 'stat_controlled', 'source_data/SFig11/controlled_176_comparisons.tsv', None),
    ('843', 'Fifteen participants contributed paired intestinal-fluid and stool measurements', 'Fifteen', 'L843_Fifteen', 'stat_regional_human_share', 'source_data/Fig5/origin_fraction_cells.tsv', None),
    ('843', 'across all five sampling groups', 'five', 'L843_five', 'stat_regional_human_share', 'source_data/Fig5/origin_fraction_cells.tsv', None),
    ('843', 'with 95% intervals from 2,000 participant-bootstrap resamples', '95%', None, '', '', 'OPEN'),
    ('843', 'with 95% intervals from 2,000 participant-bootstrap resamples', '2,000', 'L843_2000', 'stat_regional_human_share', 'source_data/Fig5/origin_fraction_cells.tsv; source_data/Fig5/four_region_tests.tsv', None),
    ('843', 'enumerated all 32,768 possible sign patterns across the fifteen participants', '32,768', 'L843_32768', 'stat_regional_human_share', 'source_data/Fig5/origin_fraction_cells.tsv', None),
    ('843', 'enumerated all 32,768 possible sign patterns across the fifteen participants', 'fifteen', 'L843_fifteen', 'stat_regional_human_share', 'source_data/Fig5/origin_fraction_cells.tsv', None),
    ('843', 'with Benjamini–Hochberg correction across the four regional contrasts', 'four', 'L843_four', 'stat_regional_human_share', 'source_data/Fig5/four_region_tests.tsv', None),
    ('843', 'compared separately with all fourteen other participants’ stool profiles', 'fourteen', 'L843_fourteen', 'stat_own_stool', f'{SDU}/own_stool_distance_matrices.tsv', None),
    ('847', 'yielding 75 profiles from the fifteen participants sampled in all five regions', '75', 'L847_75', 'stat_own_stool', f'{SDU}/own_stool_distance_matrices.tsv', None),
    ('847', 'yielding 75 profiles from the fifteen participants sampled in all five regions', 'fifteen', 'L847_fifteen', 'stat_own_stool', f'{SDU}/own_stool_distance_matrices.tsv', None),
    ('847', 'with every participant closer to their own stool and P≤1 × 10⁻⁵', '1 × 10⁻⁵', 'L847_p', 'stat_own_stool', f'{SDU}/own_stool_distance_matrices.tsv', None),
    ('847', 'the mean distance from their four capsule-region profiles to their own stool profile', 'four', 'L847_four', 'stat_own_stool', f'{SDU}/own_stool_distance_matrices.tsv', None),
    ('847', 'to each of the fourteen other participants’ stool profiles', 'fourteen', 'L847_fourteen', 'stat_own_stool', f'{SDU}/own_stool_distance_matrices.tsv', None),
    ('849', 'permuting participant labels of the stool profiles 99,999 times (seed 20260910)', '99,999', 'L849_99999', 'stat_own_stool', f'{SDU}/own_stool_distance_matrices.tsv', None),
    ('849', 'permuting participant labels of the stool profiles 99,999 times (seed 20260910)', '20260910', None, '', '', 'OPEN'),
    ('849', 'Its squared vector norm was multiplied by 15', '15', 'L849_15', 'stat_own_stool', f'{SDU}/own_stool_distance_matrices.tsv', None),
    ('849', 'total sum of squares of the 75 transformed participant–group profiles', '75', 'L849_75', 'stat_own_stool', f'{SDU}/own_stool_distance_matrices.tsv', None),
    ('855', 'All 331 CapScan and 306 PREDICT acquisitions were searched', '331', 'L855_capscan', 'stat_nonspecific', f'{SDU}/nonspecific_search_summary.tsv', 'v28f_change:347'),
    ('855', 'All 331 CapScan and 306 PREDICT acquisitions were searched', '306', 'L855_predict', 'stat_nonspecific', f'{SDU}/nonspecific_search_summary.tsv', None),
    ('855', 'including peptide lengths of 5–50 residues', '5–50', 'L855_lengths', 'stat_nonspecific', f'{SDU}/nonspecific_search_summary.tsv', None),
    ('855', 'masses of 500–5,000 Da', '500–5,000', 'L855_mass', 'stat_nonspecific', f'{SDU}/nonspecific_search_summary.tsv', None),
    ('855', '±20 ppm precursor and fragment tolerances', '20', 'L855_tol', 'stat_nonspecific', f'{SDU}/nonspecific_search_summary.tsv', None),
    ('855', 'up to two methionine oxidations', 'two', 'L855_two', 'stat_nonspecific', f'{SDU}/nonspecific_search_summary.tsv', None),
    ('855', 'Terminus analyses used peptides of 9–50 residues', '9–50', 'L855_terminus', 'stat_histidine_families; stat_nonspecific', f'source_data/Fig6/capscan_terminal_statistical_families.tsv; {SDU}/nonspecific_search_summary.tsv', None),
    ('855', 'rank-one target PSMs were required to pass peptide q≤0.01', '0.01', None, '', '', 'OPEN'),
    ('855', 'De novo predictions were analyzed using the production top-30% prediction filter', '30%', None, '', '', 'OPEN'),
    ('869', 'The PREDICT cohort comprises 93 patients for these displays', '93', 'L869_93', 'stat_histidine_families', 'source_data/Fig6/predict_histidine_statistical_families.tsv', None),
    ('873', 'The CapScan motif family included 32,400 comparisons', '32,400', 'L873_32400', 'stat_histidine_families', 'source_data/Fig6/capscan_histidine_statistical_families.tsv', None),
    ('877', 'selected within twelve cohort/compartment/origin/boundary strata', 'twelve', None, '', '', 'OPEN'),
    ('877', 'eight CapScan strata (intestinal fluid or stool', 'eight', None, '', '', 'OPEN'),
    ('877', 'and four PREDICT stool strata', 'four', 'L877_four', 'stat_strata', 'source_data/SFig14/spectra/all_strata.tsv', None),
    ('889', 'with adjusted q<0.05 considered significant for exploratory inference', '0.05', None, '', '', 'OPEN'),
    ('893', '(i) two PREDICT antibiotic endpoints', 'two', 'L893_two', 'stat_antibiotic', 'source_data/STab3/endpoints_by_definition.tsv', None),
    ('893', '(ii) four CapScan regional human-assigned-share contrasts', 'four', 'L893_four', 'stat_regional_human_share', 'source_data/Fig5/four_region_tests.tsv', None),
    ('893', '(iii) 28 PREDICT and 16 CapScan functional-profile contrasts', '28', 'L893_28', 'stat_crossomics', 'source_data/SFig11/functional_correlations_PREDICT.tsv', None),
    ('893', '(iii) 28 PREDICT and 16 CapScan functional-profile contrasts', '16', 'L893_16', 'stat_crossomics', 'source_data/Fig5/functional_correlations.tsv', None),
    ('893', '(iv) ten PREDICT and six CapScan native-profile contrasts', 'ten', 'L893_ten', 'stat_crossomics; stat_controlled', 'source_data/SFig11/native_correlations_PREDICT.tsv; source_data/SFig11/controlled_176_comparisons.tsv', None),
    ('893', '(iv) ten PREDICT and six CapScan native-profile contrasts', 'six', 'L893_six', 'stat_crossomics; stat_controlled', 'source_data/Fig5/native_correlations.tsv; source_data/SFig11/controlled_176_comparisons.tsv', None),
    ('893', '(v) nine representation-by-layer-pair families', 'nine', 'L893_nine', 'stat_individual_functions', 'source_data/STab5/crossomics_same_function_correlations.tsv', None),
    ('893', 'The seven CapScan non-DNA comparisons repeated in all fifteen participants', 'seven', 'L893_seven', 'stat_crossomics', 'source_data/Fig5/functional_correlations.tsv', None),
    ('893', 'The seven CapScan non-DNA comparisons repeated in all fifteen participants', 'fifteen', 'L893_fifteen', 'stat_crossomics', 'source_data/Fig5/functional_correlations.tsv', None),
    ('893', '(i) CapScan terminal fraction 4,536', '4,536', 'L893_4536', 'stat_histidine_families', 'source_data/Fig6/capscan_terminal_statistical_families.tsv', None),
    ('893', '(ii) residue composition 13,608', '13,608', 'L893_13608', 'stat_histidine_families', 'source_data/Fig6/capscan_terminal_statistical_families.tsv', None),
    ('893', '(iii) internal source residue 9,000 CapScan and 1,560 PREDICT', '9,000', 'L893_9000', 'stat_histidine_families', 'source_data/Fig6/capscan_histidine_statistical_families.tsv', None),
    ('893', '(iii) internal source residue 9,000 CapScan and 1,560 PREDICT', '1,560', 'L893_1560', 'stat_histidine_families', 'source_data/Fig6/predict_histidine_statistical_families.tsv', None),
    ('893', '(iv) complementary products 90 and 12', '90', 'L893_90', 'stat_histidine_families', 'source_data/Fig6/capscan_histidine_statistical_families.tsv', None),
    ('893', '(iv) complementary products 90 and 12', '12', 'L893_12', 'stat_histidine_families', 'source_data/Fig6/predict_histidine_statistical_families.tsv', None),
    ('893', '(v) CapScan spatial histidine 2,520 and motifs 32,400', '2,520', 'L893_2520', 'stat_histidine_families', 'source_data/Fig6/capscan_histidine_statistical_families.tsv', None),
    ('893', '(v) CapScan spatial histidine 2,520 and motifs 32,400', '32,400', 'L893_32400', 'stat_histidine_families', 'source_data/Fig6/capscan_histidine_statistical_families.tsv', None),
    ('909', 'identical output tables for the 58 acquisitions', '58', 'L909_58', 'stat_hstool_grouping', f'{SDU}/hstool_grouping_receipt.tsv', None),
    ('909', 'including 59,428 quantified protein groups', '59,428', 'L909_59428', 'stat_hstool_grouping', f'{SDU}/hstool_grouping_receipt.tsv', None),
    ('909', 'and 1,909,973 quantified peptide features', '1,909,973', 'L909_1909973', 'stat_hstool_grouping', f'{SDU}/hstool_grouping_receipt.tsv', None),
]


def number_rows(R):
    """Return NUMBERS.tsv rows (list of dicts) from NUMBER_SPEC and recomputed values ``R``."""
    rows = []
    for line, quote, printed, key, fn, sd, status in NUMBER_SPEC:
        val = '' if key is None or R.get(key) is None else (repr(R[key]) if isinstance(R[key], float) else str(R[key]))
        if key is not None and val == '':
            status = 'OPEN'
        rows.append(dict(v28d_line=line, quote=quote, printed_value=printed, recomputed_value=val,
                         code=f'{U}:{fn}' if fn else '', source_data=sd, status=status or 'MATCH'))
    return rows


# ----------------------------------------------------------------------------- STATS.tsv
STATS_COLUMNS = ['claim', 'v28d_line', 'test', 'sided', 'unit', 'n', 'distinct_or_repeated', 'covariates', 'correction',
                 'family', 'family_size', 'estimate', 'uncertainty', 'uncertainty_construction', 'statistic', 'df',
                 'exact_P', 'effect_size', 'effect_size_how', 'draws', 'seed', 'software_version', 'code']
#: Environments of the original producers (interpreter that ran each analysis).
ENV_BIOLOGY = 'producer: Python 3.12.4, numpy 2.2.6, scipy 1.16.3, pandas 2.3.3, statsmodels 0.14.5, scikit-learn 1.7.2'
ENV_HISTIDINE = 'producer: Python 3.10.9 (Anaconda 2023.03), numpy 1.26.4, scipy 1.12.0, pandas 2.3.3'


def _g(x, nd=6):
    """Format a float with ``nd`` significant digits for STATS.tsv."""
    return 'NA' if x is None or (isinstance(x, float) and not np.isfinite(x)) else f'{x:.{nd}g}'


def stats_rows(D):
    """Return STATS.tsv rows (list of dicts) from the results of ``collect``."""
    me = software()
    rows = []

    def add(**kw):
        r = {c: 'NA' for c in STATS_COLUMNS}
        r.update({k: str(v) for k, v in kw.items()})
        rows.append(r)

    # antibiotics (STab3 recomputes the tests from controlled inputs)
    for t in D['abx']['table'].itertuples():
        add(claim=f'PREDICT antibiotic on versus off, {t.definition}: {t.metric}', v28d_line='793-795',
            test='paired sign-flip test of the mean on-minus-off difference (Monte Carlo, plus-one)', sided='two-sided',
            unit='patient', n=t.n_patients, distinct_or_repeated='repeated: on and off states within patient (acquisitions averaged, or protein-group counts summarized by median, within patient and state)',
            covariates='none', correction='Benjamini-Hochberg', family='two antibiotic endpoints within this exposure definition',
            family_size=2, estimate=f'mean on-minus-off {_g(t.mean_on_minus_off)} {t.units}',
            uncertainty=f'95% CI {_g(t.CI95_low)} to {_g(t.CI95_high)}',
            uncertainty_construction='percentile bootstrap, 50,000 resamples of patients (seed 20260907)',
            statistic=f'mean paired difference {_g(t.mean_on_minus_off)}; {t.n_patients_on_greater_than_off} of {t.n_patients} patients higher on', df='NA',
            exact_P=f'{t.p_signflip_MC!r} (q {t.q_BH_two_endpoints_within_definition!r})', effect_size=_g(t.mean_on_minus_off),
            effect_size_how=f'mean of patient on-minus-off differences ({t.units})', draws='200,000 sign assignments; 50,000 bootstrap resamples',
            seed='20260907 (bootstrap); 20260907 + endpoint index (sign flips)', software_version=ENV_BIOLOGY,
            code='tables/STab3/panels.py:paired_inference (test, controlled inputs); text_statistics/supplementary_methods_8_11/stats.py:stat_antibiotic (family)')
    # held-out EDS prediction
    ed = D['eds']; ci = ed['ci'].set_index(['grouping', 'model'])
    for t in ed['tests'].itertuples():
        c = ci.loc[(t.grouping, t.model)]
        add(claim=f'Held-out EDS prediction, {t.grouping} grouping, {t.model} model (first split seed)', v28d_line='803-805',
            test='Spearman correlation of held-out predictions with EDS; label-permutation test with model selection repeated on each permuted label set',
            sided='two-sided', unit='patient', n=t.n, distinct_or_repeated='distinct: one observation per patient', covariates='none',
            correction='Benjamini-Hochberg', family='eight grouping-by-model tests (two grouping rules x four models), first split seed',
            family_size=ed['n_tests'], estimate=f'Spearman rho {_g(t.rho)}', uncertainty=f'95% CI {_g(c.rho_ci_low)} to {_g(c.rho_ci_high)}',
            uncertainty_construction='2.5th and 97.5th percentiles over 2,000 patient resamples of fixed held-out predictions (seed 20260918); conditional on the fitted folds',
            statistic=f'Spearman rho = {_g(t.rho)}', df='NA', exact_P=f'{t.permutation_p!r} (q {t.permutation_q_8!r})',
            effect_size=_g(t.rho), effect_size_how='Spearman correlation of held-out predictions with EDS',
            draws='199 label permutations; 2,000 patient resamples', seed='20260915 (outer split and label permutations); 20260918 (bootstrap)',
            software_version=ENV_BIOLOGY, code='text_statistics/supplementary_methods_8_11/stats.py:stat_eds_prediction')
    for t in ed['increment'].itertuples():
        add(claim=f'Held-out EDS prediction, {t.grouping} grouping: baseline-minus-{t.model} mean absolute error', v28d_line='805',
            test='interval only', sided='NA', unit='patient', n=ed['n'], distinct_or_repeated='distinct: one observation per patient',
            covariates='none', correction='none', family='NA', family_size='NA', estimate=f'MAE difference {_g(t.MAE_baseline_minus_model)}',
            uncertainty=f'95% CI {_g(t.ci_low)} to {_g(t.ci_high)}',
            uncertainty_construction='percentile bootstrap of the mean absolute-error difference, 2,000 patient resamples of fixed held-out predictions (seed 20260918)',
            statistic='NA', df='NA', exact_P='NA', effect_size=_g(t.MAE_baseline_minus_model), effect_size_how='mean of |baseline error| minus |model error| over patients (EDS units)',
            draws='2,000', seed='20260918', software_version=ENV_BIOLOGY, code='text_statistics/supplementary_methods_8_11/stats.py:stat_eds_prediction')
    # cross-omics native and functional
    fun, nat = crossomics_tables(D['_sd_root'])
    for t in nat.itertuples():
        cap = t.cohort == 'CAPSCAN'
        add(claim=f'{t.cohort} native profiles ({t.analysis}): {t.layer_a} versus {t.layer_b}', v28d_line='809-815',
            test='Spearman correlation of distance-matrix upper triangles; ' + ('whole-participant five-region block permutation' if cap else 'patient-label permutation (rows and columns together)') + ', plus-one',
            sided='two-sided', unit='participant' if cap else 'patient', n=f'{t.n_participants} {"participants" if cap else "patients"}; {t.n_profiles} profiles',
            distinct_or_repeated='repeated: five sampling-region profiles within participant' if cap else 'distinct: one profile per patient',
            covariates='sampling region (feature-wise residualization)' if cap else 'none', correction='Benjamini-Hochberg',
            family=f'{t.cohort} native profiles, {t.analysis}', family_size=6 if cap else 10, estimate=f'rho {_g(t.rho)}',
            uncertainty='none computed', uncertainty_construction='NA', statistic=f'Spearman rho = {_g(t.rho)}', df='NA',
            exact_P=f'{t.p!r} (q {t.q!r})', effect_size=_g(t.rho), effect_size_how='Spearman correlation between distance upper triangles',
            draws=int(t.permutations), seed=CROSSOMICS_SEED, software_version=ENV_BIOLOGY,
            code='text_statistics/supplementary_methods_8_11/stats.py:stat_crossomics (family and q); ' + ('figures/Fig5/panels.py:stat_g_native (test)' if cap else 'figures/SFig11/panels.py (test)'))
    for t in fun.itertuples():
        cap = t.cohort == 'CAPSCAN'
        add(claim=f'{t.cohort} functional profiles ({t.scope}, prevalence {t.prevalence}): {t.level} {t.layer_a} versus {t.layer_b}',
            v28d_line='833-837' if cap else '827',
            test='Spearman correlation of distance-matrix upper triangles; ' + ('whole-participant five-region block permutation' if cap else 'patient-label permutation (rows and columns together)') + ', plus-one',
            sided='two-sided', unit='participant' if cap else 'patient', n=f'{t.n_donors} {"participants" if cap else "patients"}; {t.n_units} profiles',
            distinct_or_repeated='repeated: five sampling-region profiles within participant' if cap else 'distinct: one profile per patient',
            covariates='sampling region (feature-wise OLS residualization)' if cap else 'none', correction='Benjamini-Hochberg',
            family=f'{t.cohort} functional profiles, {t.scope}, prevalence {t.prevalence}', family_size=16 if cap else 28,
            estimate=f'rho {_g(t.rho)}', uncertainty='none computed', uncertainty_construction='NA', statistic=f'Spearman rho = {_g(t.rho)}',
            df='NA', exact_P=f'{t.p!r} (q {t.q!r})', effect_size=_g(t.rho), effect_size_how='Spearman correlation between distance upper triangles',
            draws=CROSSOMICS_PERMUTATIONS, seed=CROSSOMICS_SEED, software_version=ENV_BIOLOGY,
            code='text_statistics/supplementary_methods_8_11/stats.py:stat_crossomics (family and q); ' + ('figures/Fig5/panels.py:stat_g_functional (test)' if cap else 'figures/SFig11/panels.py (test)'))
    for t in D['controlled']['families'].itertuples():
        cap = t.cohort == 'CAPSCAN'
        add(claim=f'Controlled distance sensitivity, {t.cohort}: {t.features}, {t.distance}, {t.adjustment}, {t.outer_correlation} outer correlation',
            v28d_line='819-823', test=f'{t.outer_correlation.capitalize()} correlation of distance upper triangles; design-preserving permutation, plus-one; positive-tail P retained separately',
            sided='two-sided', unit='participant' if cap else 'patient', n='13 participants; 65 profiles' if cap else '77 patients',
            distinct_or_repeated='repeated: five sampling-region profiles within participant' if cap else 'distinct: one profile per patient',
            covariates='sampling region (rank-space mean removal)' if (cap and t.adjustment == 'region_adjusted') else 'none',
            correction='Benjamini-Hochberg', family='this feature/distance/adjustment/outer-correlation variant', family_size=t.size,
            estimate=f'{t.significant} of {t.size} with q<0.05', uncertainty='none computed', uncertainty_construction='NA',
            statistic='correlation coefficient per comparison (source table)', df='NA', exact_P='per comparison in source_data/SFig11/controlled_176_comparisons.tsv',
            effect_size='per comparison', effect_size_how='correlation coefficient', draws=D['controlled']['permutations'], seed=CROSSOMICS_SEED,
            software_version=ENV_BIOLOGY, code='text_statistics/supplementary_methods_8_11/stats.py:stat_controlled')
    for t in D['individual']['families'].itertuples():
        add(claim=f'Individual-function correlations, {t.level}: {t.layer_a} versus {t.layer_b}', v28d_line='829',
            test='Spearman correlation across patients (scipy.stats.spearmanr, t approximation)', sided='two-sided', unit='patient',
            n='35 to 77 patients per function (positive in both layers)', distinct_or_repeated='distinct: one profile per patient',
            covariates='none', correction='Benjamini-Hochberg', family='this representation-by-layer-pair family', family_size=t.size,
            estimate=f'{t.significant} of {t.size} with q<0.05', uncertainty='none computed', uncertainty_construction='NA',
            statistic='Spearman rho per function (source table)', df='n - 2 per function (t approximation)',
            exact_P='per function in source_data/STab5/crossomics_same_function_correlations.tsv', effect_size='per function',
            effect_size_how='Spearman rho', draws='NA', seed='NA', software_version=ENV_BIOLOGY,
            code='text_statistics/supplementary_methods_8_11/stats.py:stat_individual_functions')
    for t in D['screen']['families'].itertuples():
        adj = t.model != 'unadjusted'
        add(claim=f'EDS screen, {t.layer}, {t.normalization} signal, {t.model}', v28d_line='799',
            test='partial rank correlation with EDS, t test on residual df' if adj else 'Spearman correlation with EDS, t test', sided='two-sided',
            unit='patient', n=f'43 to 85 patients per feature ({t.tested} of {t.size} features eligible)', distinct_or_repeated='distinct: earliest matched observation per patient',
            covariates='ranked human-assigned share, ranked log2 protein-group depth, antibiotic status, clinical group' if adj else 'none',
            correction='Benjamini-Hochberg', family='this layer, normalization and model (all eligible features)', family_size=t.tested,
            estimate=f'{t.significant} of {t.tested} with q<0.05', uncertainty='none computed', uncertainty_construction='NA',
            statistic='rank correlation per feature (source table)', df='n - k - 2 (k = covariate rank)' if adj else 'n - 2',
            exact_P='per feature in source_data/SFig10/all_EDS_associations.tsv', effect_size='per feature', effect_size_how='rank correlation',
            draws='NA', seed='NA', software_version=ENV_BIOLOGY, code='text_statistics/supplementary_methods_8_11/stats.py:stat_eds_screen')
    # CapScan regional
    for t in D['regional']['tests'].itertuples():
        add(claim=f'CapScan human-assigned LFQ share, {t.region} versus stool', v28d_line='843', test='exact paired sign-flip test of the mean difference (all 2^15 sign patterns)',
            sided='two-sided', unit='participant', n=f'{t.n} paired participants', distinct_or_repeated='repeated: region and stool within participant (acquisitions averaged within specimen, specimens within participant and region)',
            covariates='none', correction='Benjamini-Hochberg', family='four capsule-region versus stool contrasts', family_size=4,
            estimate=f'mean difference {_g(t.mean_difference)} percentage points', uncertainty=f'95% CI {_g(t.ci_low)} to {_g(t.ci_high)}',
            uncertainty_construction='percentile bootstrap, 2,000 participant resamples (seed 20261003)',
            statistic=f'sum of paired differences = {_g(t.statistic_sum)}; {t.n_positive} of {t.n} positive', df='NA',
            exact_P=f'{t.p!r} (q {t.q!r})', effect_size=_g(t.mean_difference), effect_size_how='mean paired difference (percentage points)',
            draws='32,768 (exact); 2,000 bootstrap', seed='20261003 (bootstrap only)', software_version=f'{ENV_BIOLOGY}; recomputed: {me}',
            code='text_statistics/supplementary_methods_8_11/stats.py:stat_regional_human_share')
    for t in D['own']['tests'].itertuples():
        add(claim=f'CapScan own-stool advantage of capsule OG profiles, OG prevalence filter {t.prevalence}', v28d_line='847-849',
            test='stool-label permutation of the mean own-versus-other distance difference, plus-one', sided='two-sided', unit='participant',
            n=f'{t.participants} participants; 75 profiles', distinct_or_repeated='repeated: four capsule profiles and one stool profile per participant',
            covariates='none', correction='none', family='NA (prevalence filters are sensitivities)', family_size='NA',
            estimate=f'mean advantage {_g(t.mean_advantage)} (Euclidean distance on Hellinger profiles)', uncertainty=f'95% CI {_g(t.ci_low)} to {_g(t.ci_high)}',
            uncertainty_construction='percentile bootstrap, 2,000 participant resamples (seed 20261003)',
            statistic=f'mean(other minus own) = {_g(t.mean_advantage)}; {t.own_closer} of {t.participants} closer to own stool', df='NA',
            exact_P=f'{t.p!r} (0 of 99,999 permuted values as extreme)', effect_size=_g(t.mean_advantage), effect_size_how='mean over participants of other-stool minus own-stool distance',
            draws='99,999 permutations; 2,000 bootstrap', seed='20260910 (permutation); 20261003 (bootstrap)', software_version=f'{ENV_BIOLOGY}; recomputed: {me}',
            code='text_statistics/supplementary_methods_8_11/stats.py:stat_own_stool')
    pt = D['partition']
    add(claim='CapScan OG profile variation: participant, sampling-group, capsule-stool and among-capsule shares (descriptive)', v28d_line='849',
        test='additive participant + sampling-group sums-of-squares partition (no test)', sided='NA', unit='participant', n='15 participants; 75 profiles (60 for capsules only)',
        distinct_or_repeated='repeated: five sampling-group profiles within participant', covariates='NA', correction='none', family='NA', family_size='NA',
        estimate='; '.join(f'{r.population} {r.component} {100 * r.fraction:.4f}%' for r in pt.itertuples()), uncertainty='none computed',
        uncertainty_construction='NA', statistic='NA', df='NA', exact_P='NA', effect_size='sum-of-squares shares', effect_size_how='component sum of squares over total sum of squares',
        draws='NA', seed='NA', software_version=ENV_BIOLOGY, code='figures/Fig5/panels.py:stat_cd_profiles (recompute from Zenodo profiles)')
    for t in D['evenness'].itertuples():
        add(claim=f'CapScan intestinal fluid versus stool, {t.metric} (named-species protein signal)', v28d_line='849 (Results 159)',
            test='Wilcoxon signed-rank, exact (scipy.stats.wilcoxon method="exact"); exact sign-flip also computed', sided='two-sided',
            unit='participant', n=f'{t.paired_participants} paired participants', distinct_or_repeated='repeated: fluid and stool within participant; per-acquisition metric, median within participant and compartment',
            covariates='none', correction='none', family='OPEN (not stated)', family_size='OPEN',
            estimate=f'participant medians {_g(t.participant_median_fluid)} (fluid) and {_g(t.participant_median_stool)} (stool)',
            uncertainty='none computed', uncertainty_construction='NA', statistic=f'W = {_g(t.wilcoxon_W)}', df='NA',
            exact_P=f'{t.wilcoxon_p_exact!r} (sign-flip {t.signflip_p_exact!r})', effect_size=_g(t.median_paired_difference),
            effect_size_how='median paired difference, fluid minus stool', draws=f'{t.signflip_patterns} sign patterns (exact)', seed='NA',
            software_version=ENV_BIOLOGY, code='figures/Fig5/panels.py:stat_evenness')
    for t in D['his']['families'].itertuples():
        add(claim=f'{t.cohort} {t.table} family {t.family}', v28d_line='869-873, 893',
            test='paired or one-sample sign-flip test of participant-level (or patient-level) differences; exact for n<=16, 200,000 assignments plus-one otherwise',
            sided='two-sided', unit='participant' if t.cohort == 'CAPSCAN' else 'patient', n='per test in source table (CapScan up to 15 participants; PREDICT up to 93 patients)',
            distinct_or_repeated='repeated: regions within participant (CapScan); specimens averaged within patient (PREDICT)', covariates='none',
            correction='Benjamini-Hochberg', family=t.family, family_size=t.size, estimate=f'{t.significant} of {t.size} with q<0.05',
            uncertainty='per test in source table', uncertainty_construction='percentile bootstrap over participants or patients (2,000, 20,000 or 50,000 multinomial resamples by endpoint)',
            statistic='mean difference per test', df='NA', exact_P='per test in source table', effect_size='per test', effect_size_how='mean observed-minus-expected or paired difference',
            draws='exact or 200,000', seed='terminal: sha256-derived per n; histidine: 20260914 + n (sign flips), 20261014 + n (bootstrap)',
            software_version=ENV_HISTIDINE, code='text_statistics/supplementary_methods_8_11/stats.py:stat_histidine_families (family and q); figures/Fig6/panels.py (tests)')
    return rows
