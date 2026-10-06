"""Supplementary Figure 5: error controls at different selection stages.

Panels a-f (one composite drawn by one producer) and the continuation page
(decomposition of the CAMPI fusion diagnostic).

Source data (``source_data/SFig05/``):

* ``fusion_pooled_current.tsv``          panels a-c (drawn values)
* ``fusion_acquisition_components.tsv``  acquisition-level components; recomputes a-c values
* ``donor_entrapment_seed_means.tsv``    panel d (drawn values)
* ``donor_recount.tsv``                  per-seed donor records; recomputes panel d values
* ``simulation_mixture_outcomes.tsv``    panels e,f dots (72 ordinary outcomes)
* ``simulation_uncertainty.tsv``         panels e,f squares and intervals
* ``simulation_all_outcomes.tsv``        all 288 simulation outcomes (72 ordinary are a subset)
* ``simulation_design.tsv``              design constants of the simulation, traced to its code
* ``fusion_decomposition_current375.tsv`` continuation page (drawn values)

The plotting code is moved without change of logic from the v30 house-style build
(``build_sfig05_style.py``; see PROVENANCE.md), which restyled the v29 producers (the
error-control block of the supplement build script for a-f and the release-packaging
script for the continued page) without changing any plotted value. Style comes from
``style/fl_style``: Söhne, ink, bold panel letters, no subpanel titles, placed sizes
160 x 189.3 mm (a-f) and 180 x 90 mm (continued).
"""

import os
from pathlib import Path

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from style import fl_style as fl
from scipy.stats import t as student_t

# ---------------------------------------------------------------- choices
Q_CUTOFF = 0.01                     # nominal peptide q, legend visible_v30f line 227; Methods line 367
FUSION_WINDOWS_AB = (5, 9)          # minimum residues of a,b (max 50), legend line 227
FUSION_MAX_RESIDUES = 50            # Methods line 367
SHORT_WINDOW = (5, 8)               # panel c, legend line 227; Methods line 367
FUSION_SEEDS = (20260907, 20260908) # two positional-permutation seeds; visible_v30f line 45 ("two randomizations"), values from data
DONOR_TIERS = ("T0_self", "T1_within_group")  # self-selected (5) and donor-selected (58), legend line 227
DONOR_SEEDS_PER_ACQUISITION = 3     # "averaged over three seeds", lines 227, 374
SIM_ARMS = ("fixed", "independent", "self")   # Methods line 376
SIM_REGIMES = ("baseline", "stress")          # Methods line 376
BOOTSTRAP_DRAWS = 20000             # Methods line 376
BOOTSTRAP_SEED = 2026091051         # not stated in the text; from the producing collect.py
BOOTSTRAP_LEVEL = 0.95              # Methods line 376
NOMINAL_LINE_PCT = 1                # dotted line, legend lines 227 and 229

# v30 restyle (build_sfig05_style.py): colours by meaning. Paired (resolved) contribution =
# FastaLake blue, unresolved = light grey; every dot in c-f is a FastaLake-pipeline control.
PAIRED_COLOR = fl.FL
UNRESOLVED_COLOR = fl.lighten(fl.GREY, 0.45)
DOT_COLOR = fl.FL
AF_SIZE_MM = (160.0, 189.3)         # placed size of a-f (docx extent)
CONT_SIZE_MM = (180.0, 90.0)        # placed size of the continued page
LETTER_X_MM = (1.0, 82.0)           # panel-letter x per column
PNG_DPI = 300                       # PNGs rasterised from the PDFs at 300 dpi


# ---------------------------------------------------------------- loading
def load(sd):
    """Read every source table of this unit.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig05``.

    Returns
    -------
    dict of pandas.DataFrame
    """
    names = ['fusion_pooled_current', 'fusion_acquisition_components', 'donor_entrapment_seed_means',
             'donor_recount', 'simulation_mixture_outcomes', 'simulation_uncertainty',
             'simulation_all_outcomes', 'simulation_design', 'fusion_decomposition_current375']
    return {n: pd.read_csv(Path(sd) / f'{n}.tsv', sep='\t') for n in names}


# ---------------------------------------------------------------- drawing (v30 build, unchanged)
def _dots(ax, groups, ylabel, percent=False):
    for i, (label, values) in enumerate(groups):
        v = np.asarray(values, dtype=float); v = v[np.isfinite(v)] * (100 if percent else 1)
        x = i + np.linspace(-.16, .16, len(v)); ax.scatter(x, v, s=8, alpha=.5, color=DOT_COLOR, lw=0, rasterized=len(v) > 100)
        if len(v): ax.hlines(np.median(v), i - .25, i + .25, color=fl.INK, lw=1.5)
    ax.set_xticks(range(len(groups)), [g[0] for g in groups]); ax.set_ylabel(ylabel)


def align_letters(fig, items, w_mm):
    """Draw panel letters with ``fl.panel`` and move each to ``x_mm`` from the left edge.

    Moved unchanged from ``_sf456_common.align_letters`` (v30 restyle helper).

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    items : list of (Axes, str, float)
    w_mm : float
        Figure width in mm.
    """
    for ax, letter, x_mm in items:
        fl.panel(ax, letter)
        pos = ax.get_position()
        t = ax._left_title
        t.set_x((x_mm / w_mm - pos.x0) / pos.width)
        t.set_ha('left')


def panel_ab(axes, fusion):
    """Draw panels a and b: fusion-linked diagnostic split into paired and unresolved parts.

    Parameters
    ----------
    axes : sequence of two matplotlib Axes
    fusion : pandas.DataFrame
        ``fusion_pooled_current.tsv``.
    """
    for i, lo in enumerate([5, 9]):
        g = fusion[fusion.q_threshold.eq(.01) & fusion.min_residues.eq(lo) & fusion.max_residues.eq(50)].copy(); g['label'] = g.cohort + '\n' + g.seed.astype(str).str[-2:]
        axes[i].bar(np.arange(len(g)), g.paired_component * 100, label='Paired', color=PAIRED_COLOR)
        axes[i].bar(np.arange(len(g)), g.unpaired_component * 100, bottom=g.paired_component * 100, label='Unresolved / unpaired', color=UNRESOLVED_COLOR)
        axes[i].set_xticks(np.arange(len(g)), g.label, rotation=45, ha='right', rotation_mode='anchor')
        axes[i].tick_params(axis='x', labelsize=6)
        axes[i].set_ylabel(f'Fusion-linked diagnostic,\n{lo}–50 residues (%)')
        axes[i].legend(frameon=False, loc='upper right', bbox_to_anchor=(1.0, 1.07), handlelength=1.2)


def panel_c(ax, fusion):
    """Draw panel c: diagnostic of the nested 5-8-residue population, one dot per seed.

    Parameters
    ----------
    ax : matplotlib Axes
    fusion : pandas.DataFrame
        ``fusion_pooled_current.tsv``.
    """
    g = fusion[fusion.q_threshold.eq(.01) & fusion.min_residues.eq(5) & fusion.max_residues.eq(8)]
    _dots(ax, [(co, g[g.cohort.eq(co)].diagnostic) for co in g.cohort.unique()], 'Fusion-linked diagnostic,\n5–8 residues (%)', True)


def panel_d(ax, don):
    """Draw panel d: explicit-peptide self-selected and donor-selected controls.

    Parameters
    ----------
    ax : matplotlib Axes
    don : pandas.DataFrame
        ``donor_entrapment_seed_means.tsv``.
    """
    g = don[don.q_threshold.eq(.01)]
    _dots(ax, [(('Self' if t == 'T0_self' else 'Donor') + '\n' + str(lo) + '–50', g[g.tier.eq(t) & g.min_length.eq(lo)].paired_estimate) for t in ['T0_self', 'T1_within_group'] for lo in [5, 9]], 'Conditional paired estimate (%)', True)


def panel_ef(axes, sim, u):
    """Draw panels e and f: realised FDP in baseline and stress simulations.

    Parameters
    ----------
    axes : sequence of two matplotlib Axes
    sim : pandas.DataFrame
        ``simulation_mixture_outcomes.tsv``.
    u : pandas.DataFrame
        ``simulation_uncertainty.tsv``.
    """
    for j, reg in enumerate(['baseline', 'stress']):
        ax = axes[j]; g = sim[sim.regime.eq(reg)]
        _dots(ax, [(arm, g[g.arm.eq(arm)].realised_FDP) for arm in ['fixed', 'independent', 'self']], f'Realised FDP, {reg}\nsynthetic mixtures (%)', True)
        for k, arm in enumerate(['fixed', 'independent', 'self']):
            r = u[u.regime.eq(reg) & u.endpoint.eq(arm)].iloc[0]
            ax.errorbar(k + .27, r['mean'] * 100, yerr=[[100 * (r['mean'] - r.bootstrap_low)], [100 * (r.bootstrap_high - r['mean'])]], fmt='s', ms=3, c=fl.INK, capsize=3)
        ax.axhline(1, c=fl.GREY, ls=':', lw=1)


def _rasterise(pdf, png, dpi=PNG_DPI):
    import fitz
    pix = fitz.open(pdf)[0].get_pixmap(dpi=dpi, alpha=False); pix.set_dpi(dpi, dpi); pix.save(png)
    return png


def figure_a_to_f(tables, out):
    """Draw panels a-f as one composite, as the v30 build did, and rasterise it at 300 dpi.

    Parameters
    ----------
    tables : dict of pandas.DataFrame
        Output of :func:`load`.
    out : pathlib.Path
        Output directory. Writes ``SFig05.pdf`` and ``SFig05.png``.

    Returns
    -------
    pathlib.Path
        The 300 dpi PNG.
    """
    out = Path(out)
    fusion, don = tables['fusion_pooled_current'], tables['donor_entrapment_seed_means']
    sim, u = tables['simulation_mixture_outcomes'], tables['simulation_uncertainty']
    fl.init()
    W, H = AF_SIZE_MM
    fig, axes = plt.subplots(3, 2, figsize=fl.mm(W, H)); axes = axes.ravel()
    fig.subplots_adjust(left=.105, right=.985, top=.965, bottom=.06, wspace=.36, hspace=.5)
    for ax in axes: fl.style(ax)
    panel_ab(axes[:2], fusion)
    panel_c(axes[2], fusion)
    panel_d(axes[3], don)
    panel_ef(axes[4:], sim, u)
    XL, XR = LETTER_X_MM
    align_letters(fig, [(axes[0], 'a', XL), (axes[1], 'b', XR), (axes[2], 'c', XL), (axes[3], 'd', XR), (axes[4], 'e', XL), (axes[5], 'f', XR)], W)
    fl.save(fig, out / 'SFig05', dpi=PNG_DPI, formats=('pdf',))
    return _rasterise(out / 'SFig05.pdf', out / 'SFig05.png')


# ---------------------------------------------------------------- continuation
def figure_continued(tables, out):
    """Draw the continuation page (CAMPI fusion decomposition) as the v30 build did.

    Parameters
    ----------
    tables : dict of pandas.DataFrame
        Output of :func:`load`.
    out : pathlib.Path
        Output directory; writes ``SFig05_continued.pdf`` and ``.png`` (300 dpi).

    Returns
    -------
    pathlib.Path
        The 300 dpi PNG.
    """
    out = Path(out)
    f = tables['fusion_decomposition_current375']
    fl.init()
    W2, H2 = CONT_SIZE_MM
    fig, ax = plt.subplots(figsize=fl.mm(W2, H2)); fl.style(ax)
    fig.subplots_adjust(left=.075, right=.99, top=.95, bottom=.2)
    x = list(range(len(f)))
    ax.bar(x, f.paired_pct, label='Weighted paired contribution', color=PAIRED_COLOR); ax.bar(x, f.unresolved_pct, bottom=f.paired_pct, label='Unresolved contribution', color=UNRESOLVED_COLOR)
    ax.axhline(1, color=fl.GREY, linestyle=':', linewidth=.8)
    ax.set_xticks(x, [f"{r['min']}–{r['max']} aa\n{r['seed']}" for _, r in f.iterrows()]); ax.set_ylabel('Fusion diagnostic (%)'); ax.set_xlabel('Peptide length and control seed')
    ax.legend(frameon=False)
    for i, r in f.iterrows(): ax.text(i, r.total_pct + .18, f'{r.total_pct:.2f}', ha='center', fontsize=6.5)
    ax.set_ylim(0, 14)
    fl.save(fig, out / 'SFig05_continued', dpi=PNG_DPI, formats=('pdf',))
    return _rasterise(out / 'SFig05_continued.pdf', out / 'SFig05_continued.png')


# ---------------------------------------------------------------- statistics
def stat_ac_fusion_pooled(components):
    """Recompute the pooled fusion diagnostics of panels a-c from acquisition components.

    Cohort summaries sum acquisition-level numerators and denominators before division
    (Methods line 715). Statistical unit: acquisition, summed within cohort and seed
    (CAMPI 15, healthy stool 58, PREDICT 306, CapScan 331 acquisitions; two seeds).

    Parameters
    ----------
    components : pandas.DataFrame
        ``fusion_acquisition_components.tsv``.

    Returns
    -------
    pandas.DataFrame
        One row per cohort, seed, q threshold and length window with ``acquisitions``,
        summed counts, ``diagnostic``, ``paired_component``, ``unpaired_component`` and
        ``median_acquisition_diagnostic``.
    """
    c = components.copy()
    c['cohort'] = c.cohort.replace({'CAMPI_MOCK': 'CAMPI'})
    fk = ['cohort', 'seed', 'q_threshold', 'min_residues', 'max_residues']
    ff = c.groupby(fk)[['n_canonical', 'paired_numerator', 'n_unpaired_accepted', 'n_entrapment_q_tied']].sum().reset_index()
    ff['acquisitions'] = c.groupby(fk)['sample'].nunique().to_numpy()
    ff['diagnostic'] = (ff.paired_numerator + ff.n_unpaired_accepted) / ff.n_canonical.replace(0, np.nan)
    ff['paired_component'] = ff.paired_numerator / ff.n_canonical.replace(0, np.nan)
    ff['unpaired_component'] = ff.n_unpaired_accepted / ff.n_canonical.replace(0, np.nan)
    acq = (c.paired_numerator + c.n_unpaired_accepted) / c.n_canonical.replace(0, np.nan)
    ff['median_acquisition_diagnostic'] = acq.groupby([c[k] for k in fk]).median().to_numpy()
    return ff


def compare_pooled(recomputed, drawn):
    """Compare recomputed pooled fusion values with the drawn table.

    Returns
    -------
    pandas.DataFrame
        Per column: rows compared and maximum absolute difference.
    """
    fk = ['cohort', 'seed', 'q_threshold', 'min_residues', 'max_residues']
    j = drawn.merge(recomputed, on=fk, suffixes=('_drawn', '_recomputed'), how='outer', indicator=True)
    rows = [dict(column='_merge', rows=len(j), max_abs_difference=float((j._merge != 'both').sum()))]
    for col in ['n_canonical', 'paired_numerator', 'n_unpaired_accepted', 'diagnostic', 'paired_component',
                'unpaired_component', 'acquisitions', 'median_acquisition_diagnostic']:
        d = (j[col + '_drawn'] - j[col + '_recomputed']).abs()
        rows.append(dict(column=col, rows=int(d.notna().sum()), max_abs_difference=float(d.max())))
    return pd.DataFrame(rows)


def stat_cont_fusion_decomposition(components, q=Q_CUTOFF):
    """Recompute the CAMPI fusion decomposition (continuation page; Methods line 715).

    Ported from the independent recount script: sums E, E0, B (strictly better plus tied
    q), U and N over the 15 CAMPI acquisitions of each seed; paired = (E+E0+2B)/N,
    unresolved = U/N. Statistical unit: acquisition, summed within seed (n = 15 per seed,
    two seeds); no test.

    Parameters
    ----------
    components : pandas.DataFrame
        ``fusion_acquisition_components.tsv`` (rows with cohort ``CAMPI_MOCK``).
    q : float
        Peptide q threshold.

    Returns
    -------
    pandas.DataFrame
        Columns as in ``fusion_decomposition_current375.tsv``.
    """
    d = components[components.cohort.eq('CAMPI_MOCK') & components.q_threshold.eq(q)]
    rows = []
    for (seed, lo, hi), g in d.groupby(['seed', 'min_residues', 'max_residues']):
        n = int(g.n_canonical.sum()); e = int(g.n_paired_entrapment.sum()); e0 = int(g.n_entrapment_alone.sum()); b = int((g.n_entrapment_q_strict_better + g.n_entrapment_q_tied).sum()); u = int(g.n_unpaired_accepted.sum())
        num = e + e0 + 2 * b
        assert num == g.paired_numerator.sum()
        rows.append(dict(seed=int(seed), min=int(lo), max=int(hi), n=n, E=e, E0=e0, B=b, q_ties=int(g.n_entrapment_q_tied.sum()), paired_numerator=num, U=u, paired_pct=100 * num / n, unresolved_pct=100 * u / n, total_pct=100 * (num + u) / n, unresolved_share_pct=100 * u / (num + u), acquisitions=g['sample'].nunique()))
    return pd.DataFrame(rows)


def stat_d_donor_seed_means(recount):
    """Average the three entrapment seeds within target acquisition (panel d; Methods line 729).

    Ported from the control-reconciliation script. Statistical unit: acquisition
    (58 donor-selected targets, 5 self-selected targets), three seeds averaged within each.

    Parameters
    ----------
    recount : pandas.DataFrame
        ``donor_recount.tsv`` (one row per target, donor, seed, length and q).

    Returns
    -------
    pandas.DataFrame
        Seed means with the columns of ``donor_entrapment_seed_means.tsv``.
    """
    gkeys = ['target', 'donor', 'tier', 'group', 'min_length', 'q_threshold']
    g = recount.groupby(gkeys, dropna=False)
    assert (g.size() == DONOR_SEEDS_PER_ACQUISITION).all()
    return g[['n_target', 'n_entrapment', 'paired_estimate', 'lower_fraction', 'combined_estimate']].mean().reset_index()


def stat_d_donor_summary(seed_means, q=Q_CUTOFF):
    """Median and quartiles of acquisition seed means per tier and length (line 65, panel d).

    Statistical unit: acquisition (n = 58 donor-selected, 5 self-selected); descriptive,
    no test. Quartiles use pandas linear interpolation, as in the producer.

    Returns
    -------
    pandas.DataFrame
        tier, min_length, acquisitions, median, q25, q75, minimum, maximum.
    """
    rows = []
    for (tier, lo), x in seed_means[seed_means.q_threshold.eq(q)].groupby(['tier', 'min_length']):
        z = x.paired_estimate
        rows.append(dict(tier=tier, min_length=lo, acquisitions=len(x), median=z.median(), q25=z.quantile(.25), q75=z.quantile(.75), minimum=z.min(), maximum=z.max()))
    return pd.DataFrame(rows)


def _summarize(values, rng):
    """Mean, percentile bootstrap and t interval of a vector (moved unchanged from collect.py)."""
    values = np.asarray(values, dtype=float)
    n = len(values)
    mean = float(values.mean())
    if n < 2:
        return dict(n_units=n, mean=mean, bootstrap_low=None, bootstrap_high=None, t_low=None, t_high=None)
    draws = values[rng.integers(n, size=(BOOTSTRAP_DRAWS, n))].mean(axis=1)
    half = float(student_t.ppf(0.975, n - 1) * values.std(ddof=1) / np.sqrt(n))
    return dict(n_units=n, mean=mean, bootstrap_low=float(np.quantile(draws, 0.025)),
                bootstrap_high=float(np.quantile(draws, 0.975)), t_low=mean - half, t_high=mean + half)


def stat_ef_simulation_uncertainty(primary):
    """Recompute mixture means and 95% whole-mixture bootstrap intervals (panels e,f; line 733).

    Moved from the simulation collect step with the same random generator, seed and draw
    order. Statistical unit: synthetic mixture (n = 12 per regime); bootstrap resamples
    whole mixtures; the paired difference self minus independent keeps mixture pairing.

    Parameters
    ----------
    primary : pandas.DataFrame
        ``simulation_all_outcomes.tsv`` in its original row order (288 rows).

    Returns
    -------
    pandas.DataFrame
        Columns as in ``simulation_uncertainty.tsv``.
    """
    rng = np.random.default_rng(BOOTSTRAP_SEED)
    estimates = []
    ordinary = primary[primary.arm.isin(['fixed', 'self', 'independent'])]
    for (regime, arm), frame in ordinary.groupby(['regime', 'arm'], sort=True):
        estimates.append(dict(regime=regime, endpoint=arm, **_summarize(frame.realised_FDP, rng),
                              total_discoveries=int(frame.discoveries.sum()), total_false=int(frame.false_discoveries.sum()),
                              pooled_FDP=float(frame.false_discoveries.sum() / max(frame.discoveries.sum(), 1))))
    paired = ordinary.pivot(index=['regime', 'unit_id'], columns='arm', values='realised_FDP')
    for regime, frame in paired.groupby(level='regime', sort=True):
        estimates.append(dict(regime=regime, endpoint='self_minus_independent', **_summarize(frame.self - frame.independent, rng),
                              total_discoveries=None, total_false=None, pooled_FDP=None))
    return pd.DataFrame(estimates)


def stat_ef_fdp_identity(sim):
    """Check FDP = false discoveries / max(discoveries, 1) for all 72 ordinary outcomes.

    Returns
    -------
    float
        Maximum absolute difference.
    """
    return float((sim.realised_FDP - sim.false_discoveries / sim.discoveries.clip(lower=1)).abs().max())


def stat_stress_self_minus_fixed(sim):
    """Mean paired difference self minus fixed under stress (Methods line 733).

    Statistical unit: synthetic mixture, n = 12 stress mixtures; point estimate only
    (no interval is printed for it).

    Returns
    -------
    float
        Mean difference in FDP (fraction).
    """
    s = sim[sim.regime.eq('stress')].pivot(index='unit_id', columns='arm', values='realised_FDP')
    return float((s.self - s.fixed).mean())


def compare_frames(a, b, keys):
    """Maximum absolute difference of every shared numeric column of two frames joined on ``keys``."""
    j = a.merge(b, on=keys, suffixes=('_a', '_b'), how='outer', indicator=True)
    out = {'unmatched_rows': int((j._merge != 'both').sum())}
    for col in a.columns:
        if col in keys or col + '_b' not in j:
            continue
        x = pd.to_numeric(j[col + '_a'], errors='coerce') - pd.to_numeric(j[col + '_b'], errors='coerce')
        out[col] = float(x.abs().max()) if x.notna().any() else None
    return out


def numbers(tables):
    """Compute every value listed in NUMBERS.tsv, in file order (keys name visible_v30f lines).

    Returns
    -------
    list of (key, value)
    """
    fu = tables['fusion_decomposition_current375']
    sim, u = tables['simulation_mixture_outcomes'], tables['simulation_uncertainty']
    dsg = tables['simulation_design'].set_index(['parameter', 'regime']).value
    sm = stat_d_donor_seed_means(tables['donor_recount'])
    ds = stat_d_donor_summary(sm).set_index(['tier', 'min_length'])
    ur = u.set_index(['regime', 'endpoint'])
    w = lambda lo, hi: fu[fu['min'].eq(lo) & fu['max'].eq(hi)]
    comp = tables['fusion_acquisition_components']
    q01 = sm[sm.q_threshold.eq(Q_CUTOFF)]
    vals = [
        ('L30 donor median 5-50 %', 100 * ds.loc[('T1_within_group', 5), 'median']),
        ('L30 donor median 9-50 %', 100 * ds.loc[('T1_within_group', 9), 'median']),
        ('L45 CAMPI 9-50 mean of seeds %', w(9, 50).total_pct.mean()),
        ('L45 seeds', fu.seed.nunique()),
        ('L46 baseline self mean %', 100 * ur.loc[('baseline', 'self'), 'mean']),
        ('L46 baseline fixed mean %', 100 * ur.loc[('baseline', 'fixed'), 'mean']),
        ('L46 stress self mean %', 100 * ur.loc[('stress', 'self'), 'mean']),
        ('L46 stress self-independent pp', 100 * ur.loc[('stress', 'self_minus_independent'), 'mean']),
        ('L227 q', Q_CUTOFF),
        ('L227 self-selected acquisitions', q01[q01.tier.eq('T0_self')].target.nunique()),
        ('L227 donor-selected acquisitions', q01[q01.tier.eq('T1_within_group')].target.nunique()),
        ('L227 seeds per acquisition', int(tables['donor_recount'].groupby(['target', 'donor', 'tier', 'min_length', 'q_threshold']).seed.nunique().unique()[0])),
        ('L227 mixtures per regime', int(sim.groupby('regime').unit_id.nunique().unique()[0])),
        ('L227 nominal line %', NOMINAL_LINE_PCT),
        ('L229 CAMPI acquisitions per seed', int(comp[comp.cohort.eq('CAMPI_MOCK')].groupby('seed')['sample'].nunique().unique()[0])),
        ('L229 q', Q_CUTOFF),
        ('L229 nominal line %', NOMINAL_LINE_PCT),
        ('L367 seeds', fu.seed.nunique()),
        ('L367 5-50 total min %', w(5, 50).total_pct.min()),
        ('L367 5-50 total max %', w(5, 50).total_pct.max()),
        ('L367 9-50 total min %', w(9, 50).total_pct.min()),
        ('L367 9-50 total max %', w(9, 50).total_pct.max()),
        ('L367 5-8 paired min %', w(5, 8).paired_pct.min()),
        ('L367 5-8 paired max %', w(5, 8).paired_pct.max()),
        ('L367 5-8 unresolved min %', w(5, 8).unresolved_pct.min()),
        ('L367 5-8 unresolved max %', w(5, 8).unresolved_pct.max()),
        ('L367 5-8 total min %', w(5, 8).total_pct.min()),
        ('L367 5-8 total max %', w(5, 8).total_pct.max()),
        ('L367 5-8 unresolved share min %', w(5, 8).unresolved_share_pct.min()),
        ('L367 5-8 unresolved share max %', w(5, 8).unresolved_share_pct.max()),
        ('L367 9-50 paired min %', w(9, 50).paired_pct.min()),
        ('L367 9-50 paired max %', w(9, 50).paired_pct.max()),
        ('L367 9-50 unresolved min %', w(9, 50).unresolved_pct.min()),
        ('L367 9-50 unresolved max %', w(9, 50).unresolved_pct.max()),
        ('L373 donor-selected acquisitions', q01[q01.tier.eq('T1_within_group')].target.nunique()),
        ('L374 seeds per acquisition', int(tables['donor_recount'].groupby(['target', 'donor', 'tier', 'min_length', 'q_threshold']).seed.nunique().unique()[0])),
        ('L374 self-selected acquisitions', q01[q01.tier.eq('T0_self')].target.nunique()),
        ('L376 spectra per acquisition', dsg[('spectra_per_acquisition', 'all')]),
        ('L376 baseline mixtures', sim[sim.regime.eq('baseline')].unit_id.nunique()),
        ('L376 stress mixtures', sim[sim.regime.eq('stress')].unit_id.nunique()),
        ('L376 generating min length', dsg[('generating_min_length', 'all')]),
        ('L376 generating max length', dsg[('generating_max_length', 'all')]),
        ('L376 reference proteins', dsg[('reference_proteins', 'all')]),
        ('L376 ordinary outcomes', len(sim)),
        ('L376 bootstrap draws', BOOTSTRAP_DRAWS),
        ('L376 all outcomes', len(tables['simulation_all_outcomes'])),
        ('L376 baseline fixed %', 100 * ur.loc[('baseline', 'fixed'), 'mean']),
        ('L376 baseline independent %', 100 * ur.loc[('baseline', 'independent'), 'mean']),
        ('L376 baseline self %', 100 * ur.loc[('baseline', 'self'), 'mean']),
        ('L376 stress fixed %', 100 * ur.loc[('stress', 'fixed'), 'mean']),
        ('L376 stress independent %', 100 * ur.loc[('stress', 'independent'), 'mean']),
        ('L376 stress self %', 100 * ur.loc[('stress', 'self'), 'mean']),
        ('L376 stress absent %', 100 * (1 - dsg[('present_spectra', 'stress')] / dsg[('spectra_per_acquisition', 'all')])),
        ('L376 baseline absent %', 100 * (1 - dsg[('present_spectra', 'baseline')] / dsg[('spectra_per_acquisition', 'all')])),
        ('L376 stress self-independent pp', 100 * ur.loc[('stress', 'self_minus_independent'), 'mean']),
        ('L376 stress self-independent CI low pp', 100 * ur.loc[('stress', 'self_minus_independent'), 'bootstrap_low']),
        ('L376 stress self-independent CI high pp', 100 * ur.loc[('stress', 'self_minus_independent'), 'bootstrap_high']),
        ('L376 stress self-fixed pp', 100 * stat_stress_self_minus_fixed(sim)),
        ('L376 reference proteins, nearest thousand', int(round(dsg[('reference_proteins', 'all')], -3))),
        ('bootstrap level %', int(round(100 * BOOTSTRAP_LEVEL))),
    ]
    by = dict(vals)
    words = [
        ('L45 seeds word', as_word(by['L45 seeds'])),
        ('L227 self-selected word', as_word(by['L227 self-selected acquisitions'])),
        ('L227 seeds word', as_word(by['L227 seeds per acquisition'])),
        ('L229 CAMPI acquisitions word', as_word(by['L229 CAMPI acquisitions per seed'], capital=True)),
        ('L367 seeds word', as_word(by['L367 seeds'])),
        ('L374 seeds word', as_word(by['L374 seeds per acquisition'], capital=True)),
        ('L374 self-selected word', as_word(by['L374 self-selected acquisitions'])),
        ('L376 baseline mixtures word', as_word(by['L376 baseline mixtures'])),
        ('L376 stress mixtures word', as_word(by['L376 stress mixtures'])),
    ]
    return vals + words


_WORDS = {2: 'two', 3: 'three', 5: 'five', 12: 'twelve', 15: 'fifteen'}


def as_word(n, capital=False):
    """Spell a small count as printed in the manuscript (``5`` -> ``five``)."""
    w = _WORDS[int(n)]
    return w.capitalize() if capital else w
