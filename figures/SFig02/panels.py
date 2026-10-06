"""Supplementary Figure 2: specificity and sensitivity of peptide-guided database selection.

Panels a-d (one producer figure) and the continuation panels e,f (one producer figure).

Source data (``source_data/SFig02/``):

* ``selection_length_lookup.tsv``  per arm, seed and length: distinct queries and matches (a,b)
* ``lookup_query_outcomes.tsv``    one row per distinct query, arm and seed (recount of a,b)
* ``selection_grid_acquisitions.tsv``  19 arms x 58 acquisitions (c,d)
* ``selection_grid_summary.tsv``   per-arm medians and exclusive quartiles (c,d)
* ``bags_bounded48_per_acquisition.tsv``  7 arms x 48 acquisitions, gains and losses (e)
* ``bags_bounded48_orderings.tsv``  per acquisition: lookup-eligible predictions by number of
  bag orderings (Supplementary Methods 2 coverage figure)
* ``bags_fragment_audit_summary.tsv``  six-acquisition fragment audit (f)

The plotting code is moved from the original producers unchanged (see PROVENANCE.md):
panels a-d from ``build_revision.py`` (selection block) followed by the title-removal and
300-dpi rasterisation of ``prepare_supplement_figures.py``; panels e,f from the plotting
block of ``audit_and_plot.py`` as executed by ``integrate_v18.py`` (panel letters e,f),
followed by the 400-dpi re-export of ``reexport.py``. Only file paths were changed. Both
producers used matplotlib's default typeface (DejaVu Sans), not the house style, so the
functions below set the producers' own rcParams inside a context and do not use fl_style.

All statistics in this unit are descriptive counts, fractions, medians and quartiles. No
inferential test is made.

v30: ``panel_abcd_v30``, ``panel_ef_v30`` and ``build_v30`` redraw both canvases in the house
style (fl_style from the repository's ``style/``) and prove every plotted value identical to the
producers above; moved unchanged from the v30 rebuild ``build_sfig02_style.py`` together with the
restyle helpers of ``restyle_common.py``.
"""

import io
import math
import re
import collections
import statistics as st
from pathlib import Path

import matplotlib as mpl
import numpy as np
import pandas as pd

# ---- choices (each is reported in README.md, Choices) -------------------------------
SCRAMBLE_SEEDS = (20260901, 20260913, 20260914, 20260915, 20260916)  # Supp. Methods 6, v28d line 665
MIN_LOOKUP_LENGTH = 9            # production minimum query length, v28d line 557
MAX_LOOKUP_LENGTH = 50           # production maximum query length, v28d line 557
TOP_SCORE_FRACTION = 0.30        # production top-scoring fraction, v28d line 557
LENGTH_GRID = (6, 7, 8, 9, 10, 11, 12)          # sensitivity grid, v28d line 679
SCORE_GRID = (0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00)  # v28d line 679
QUARTILE_METHOD = "exclusive"    # statistics.quantiles default; legend v28d line 409
PEPTIDE_Q = 0.01                 # nominal peptide q threshold, v28d lines 413 and 579
ORDERING_BOUNDS = (2, 4, 8)      # bag-ordering bounds, v28d lines 413 and 579
FRAGMENT_TOLERANCE_PPM = 20      # upstream fragment audit only (not used here), v28d line 581
FRAGMENT_CHARGES = (1, 2)        # upstream fragment audit only (not used here), v28d line 581
SHORT_LENGTH_CUTOFF = 13         # "no match at >=13 residues", v28d line 59
PDF_RASTER_DPI = 300             # prepare_supplement_figures.py rasterisation of panels a-d
EF_PLACED_WIDTH_IN = 6.1         # placed width of the continuation in the manuscript (154.94 mm)
EF_TARGET_DPI = 400              # reexport.py TARGET_DPI
# Title written by build_revision.py and removed again by prepare_supplement_figures.py.
NUMBER_WORDS = {1: 'one', 2: 'two', 3: 'three', 5: 'five', 6: 'six', 7: 'seven'}  # for numbers printed as words
ABCD_TITLE = "Evidence selection and exact-lookup specificity"
ABCD_TITLE_BBOX = (108.6468734741211, 6.780389785766602, 402.5669860839844, 19.59539031982422)

# build_revision.py, lines 13-14 (rcParams and colours of the a-d producer)
_ABCD_RC = {'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.titlesize': 9, 'axes.labelsize': 8,
            'xtick.labelsize': 7, 'ytick.labelsize': 7, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
            'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 180}
COL = ['#2878a3', '#df8f44', '#619b68', '#b35872', '#8277ac']


# ---- loading -----------------------------------------------------------------------
def load(sd, name, exact=False):
    """Read one source_data table.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig02``.
    name : str
        File name.
    exact : bool
        Parse floats with round-trip precision (used by the statistics functions). The
        plotting functions read as the original producers did (pandas default parser).

    Returns
    -------
    pandas.DataFrame
    """
    kw = {"float_precision": "round_trip"} if exact else {}
    return pd.read_csv(Path(sd) / name, sep="\t", **kw)


# ---- panels a-d ----------------------------------------------------------------------
def panel_abcd(sd, out_pdf):
    """Draw panels a-d exactly as ``build_revision.py`` (selection block) did.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig02``.
    out_pdf : pathlib.Path
        Output PDF (the producer's ``SELECTION_SPECIFICITY.pdf``, with its in-image title).

    Returns
    -------
    pathlib.Path
        ``out_pdf``.
    """
    import matplotlib.pyplot as plt
    q = load(sd, 'selection_length_lookup.tsv')
    a = load(sd, 'selection_grid_summary.tsv')
    matched = a[a.kind.eq('length') | a.arm.eq('tp0.30')].sort_values('min_length').copy()

    def title(ax, l, t):
        ax.set_title(l + '  ' + t, loc='left', pad=9)

    with mpl.rc_context():
        mpl.rcdefaults()
        plt.rcParams.update(_ABCD_RC)
        fig, ax = plt.subplots(2, 2, figsize=(7.1, 6.0)); ax = ax.ravel(); real = q[q.arm.eq('REAL')]
        ax[0].plot(real.length, real.fraction_matched * 100, 'o-', c=COL[0], label='Authentic')
        for _, g in q[~q.arm.eq('REAL')].groupby('arm'): ax[0].plot(g.length, g.fraction_matched * 100, c=COL[1], alpha=.6)
        ax[0].plot([], [], c=COL[1], label='Five scrambled seeds'); ax[0].set(xlabel='Peptide length (residues)', ylabel='Distinct queries matching (%)', xlim=(8.5, 50.5)); fig.legend(*ax[0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(.5, .947), ncol=2, frameon=False, fontsize=7); title(ax[0], 'a', 'Exact reservoir lookup')
        inset = ax[0].inset_axes([.45, .15, .5, .20])
        for _, g in q[~q.arm.eq('REAL')].groupby('arm'): inset.plot(g.length, g.fraction_matched * 100, 'o-', ms=2, c=COL[1], alpha=.6)
        inset.set(xlim=(8.8, 13.2), ylim=(-.3, 10), xticks=[9, 11, 13], yticks=[0, 5, 10]); inset.tick_params(labelsize=7); inset.set_title('Scrambled: 9–13 residues', fontsize=7)
        ax[1].plot(real.length, real.distinct_queries, 'o-', c=COL[0]); ax[1].set(xlabel='Peptide length (residues)', ylabel='Distinct authentic queries', xlim=(8.5, 50.5)); title(ax[1], 'b', 'Length-specific denominator')
        for axis, kind, x in [(ax[2], 'length', 'min_length'), (ax[3], 'score', 'top_percent')]:
            g = (matched if kind == 'length' else a[a.kind.eq(kind)]).sort_values(x); xx = g[x].to_numpy() * (100 if kind == 'score' else 1)
            axis.errorbar(xx, g.peptides_q01_median / 1000, yerr=np.array([g.peptides_q01_median - g.peptides_q01_q1, g.peptides_q01_q3 - g.peptides_q01_median]) / 1000, fmt='o-', c=COL[0], label='Accepted peptides')
            axis.plot(xx, g.razor_selected_median / 1000, 's--', c=COL[1], label='Selected proteins (median)'); axis.set(xlabel='Minimum lookup length (residues)' if kind == 'length' else 'Top-scoring predictions (%)', ylabel='Count (thousands)', ylim=(0, 48)); axis.legend(frameon=False, fontsize=7, loc='upper right')
            axis.axvline(9 if kind == 'length' else 30, c='grey', lw=.8, ls=':')
            title(ax[2], 'c', 'Length sensitivity'); title(ax[3], 'd', 'Score-fraction sensitivity')
        # build_revision.py save(): in-image title, tight layout, vector output
        fig.suptitle(ABCD_TITLE, fontsize=11, fontweight='bold')
        fig.tight_layout(rect=(0, 0, 1, .96), pad=1.6)
        fig.savefig(out_pdf)
        plt.close(fig)
    return Path(out_pdf)


def compose_abcd(pdf_in, pdf_out, png_out):
    """Remove the in-image title and rasterise at 300 dpi (``prepare_supplement_figures.py``).

    Parameters
    ----------
    pdf_in : pathlib.Path
        Output of :func:`panel_abcd`.
    pdf_out, png_out : pathlib.Path
        Title-free PDF and its 300-dpi PNG (the image placed in the manuscript).

    Returns
    -------
    dict
        Title bbox used and whether every numeric label was unchanged.
    """
    import fitz
    doc = fitz.open(pdf_in); p = doc[0]; before = p.get_text()
    spans = [s for b in p.get_text('dict')['blocks'] for l in b.get('lines', []) for s in l['spans']
             if s['text'] == ABCD_TITLE]
    assert len(spans) == 1, 'title span not found once'
    bbox = tuple(spans[0]['bbox'])
    assert all(abs(x - y) < 1e-3 for x, y in zip(bbox, ABCD_TITLE_BBOX)), ('title moved', bbox)
    assert bbox[1] < 20
    p.add_redact_annot(fitz.Rect(bbox), fill=(1, 1, 1)); p.apply_redactions(images=0, graphics=0)
    assert ABCD_TITLE not in p.get_text()
    before_nums = collections.Counter(re.findall(r'\d[\d,.]*', before.replace(ABCD_TITLE, '')))
    after_nums = collections.Counter(re.findall(r'\d[\d,.]*', p.get_text()))
    assert before_nums == after_nums, 'numeric label changed'
    doc.save(pdf_out, garbage=4, deflate=True)
    p = fitz.open(pdf_out)[0]; p.get_pixmap(dpi=PDF_RASTER_DPI, alpha=False).save(png_out)
    return dict(title_bbox=bbox, numeric_labels_unchanged=True)


# ---- panels e,f ----------------------------------------------------------------------
def bags_summary(sd):
    """Sum gained, lost and net identifications per arm over the 48 acquisitions.

    Same aggregation as ``audit_and_plot.py`` line 19. Unit: peptide-acquisition
    identification; n = 48 acquisitions per arm (336 arm-acquisition rows).

    Returns
    -------
    pandas.DataFrame
        Index ``arm``; columns ``gained``, ``lost``, ``net``.
    """
    f = load(sd, 'bags_bounded48_per_acquisition.tsv')
    assert len(f) == 336 and not f.duplicated(['cohort', 'sample', 'arm']).any()
    assert (f.gained - f.lost == f.net).all()
    return f.groupby('arm')[['gained', 'lost', 'net']].sum()


def fragment_subset(sd):
    """Gained composition-preserving rearrangement pairs of the bag-vs-exact contrast.

    Same filter as ``audit_and_plot.py`` line 24. Unit: same-spectrum comparison pair;
    n = 215 pairs in six acquisitions.
    """
    fr = load(sd, 'bags_fragment_audit_summary.tsv')
    return fr[(fr.contrast == 'bags_vs_exact') & (fr.change == 'gained') & (fr.relation == 'composition_preserving_rearrangement')]


def panel_ef(sd, out_pdf):
    """Draw panels e,f exactly as ``audit_and_plot.py`` (via ``integrate_v18.py``) did.

    The bar lengths of panel f (193 and 22) are literals in the original plotting code; the
    original guarded them with an assertion against the audit table, kept here.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig02``.
    out_pdf : pathlib.Path
        Output PDF (the producer's ``Supplementary_Figure2ef_bags.pdf``).

    Returns
    -------
    pathlib.Path
    """
    import matplotlib.pyplot as plt
    summary = bags_summary(sd)
    fr = fragment_subset(sd)
    assert fr.spectrum_events.sum() == 215 and fr.zero_diagnostic.sum() == 193  # audit_and_plot.py line 24
    with mpl.rc_context():
        mpl.rcdefaults()
        fig, axes = plt.subplots(2, 1, figsize=(9, 6.8), gridspec_kw={'height_ratios': [4, 1.3]}); ax = axes[0]
        order = ['bags_le2', 'bags_le4', 'bags_le8', 'bags_all', 'mass_all', 'comb_all']; labels = ['Bags: ≤2 orderings', 'Bags: ≤4 orderings', 'Bags: ≤8 orderings', 'Bags: unrestricted', 'Mass substitutions', 'Bags + mass']; ss = summary.loc[order]; yy = np.arange(len(order)); ax.barh(yy, ss.gained, color='#0072B2', label='Gained'); ax.barh(yy, -ss.lost, color='#D55E00', label='Lost'); ax.axvline(0, color='.2', lw=.8); ax.set_yticks(yy, labels); ax.invert_yaxis(); ax.set_xlim(-42000, 48000)
        for i, (_, r) in enumerate(ss.iterrows()): ax.text(-r.lost - 600, i, f'{int(r.lost):,}', ha='right', va='center', fontsize=8); ax.text(r.gained + 600, i, f'{int(r.gained):,}', va='center', fontsize=8); ax.text(47500, i, f'{int(r.net):+,}', ha='right', va='center', fontsize=8, fontweight='bold')
        ax.set_xticks([-30000, 0, 30000], ['30,000 lost', '0', '30,000 gained']); ax.set_xlabel('Peptide–acquisition identifications relative to exact lookup'); ax.set_title('e  Completed matched comparison: 48 acquisitions', loc='left'); ax.legend(loc='upper left', frameon=False, fontsize=8); ax.text(.96, 1.035, 'Net', transform=ax.transAxes, ha='right', fontsize=9)
        bx = axes[1]; bx.barh([0], [193], color='#999999', label='No diagnostic ion'); bx.barh([0], [22], left=[193], color='#009E73', label='≥1 diagnostic ion'); bx.text(96.5, 0, '193', ha='center', va='center'); bx.text(204, 0, '22', ha='center', va='center'); bx.set_xlim(0, 215); bx.set_yticks([]); bx.set_xticks([0, 50, 100, 150, 200, 215]); bx.set_xlabel('Order-related same-spectrum comparison pairs'); bx.set_title('f  Separate fragment audit: six acquisitions', loc='left'); bx.legend(loc='upper center', bbox_to_anchor=(.5, -.55), ncol=2, frameon=False, fontsize=8)
        fig.tight_layout(h_pad=2.0, rect=[0, .07, 1, 1]); fig.savefig(out_pdf); plt.close(fig)
    return Path(out_pdf)


def compose_ef(pdf_in, png_out):
    """Render the e,f PDF at 400 dpi of its placed width (``reexport.py`` ``render``).

    Returns
    -------
    tuple
        Output width and height in pixels.
    """
    import fitz
    from PIL import Image
    width_px = math.ceil(EF_PLACED_WIDTH_IN * EF_TARGET_DPI)
    pg = fitz.open(pdf_in)[0]
    zoom = width_px / pg.rect.width
    pix = pg.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
    im = Image.open(io.BytesIO(pix.tobytes('png')))
    im.save(png_out, dpi=(EF_TARGET_DPI, EF_TARGET_DPI), optimize=True)
    return im.size


# ---- statistics (all descriptive) ----------------------------------------------------
def stat_ab_lookup(sd):
    """Length-specific exact-lookup fractions, recounted from per-query outcomes.

    Unit: distinct query (I/L-normalised) within arm and seed, one healthy-stool
    acquisition; authentic arm plus five scrambled seeds. The per-query table is first
    aggregated and checked against ``selection_length_lookup.tsv`` (every arm, seed and
    length, denominators and numerators).

    Returns
    -------
    dict
        Recomputed values keyed by name.
    """
    q = load(sd, 'lookup_query_outcomes.tsv', exact=True)
    t = load(sd, 'selection_length_lookup.tsv', exact=True)
    q['seed'] = q['seed'].astype('Int64'); t['seed'] = t['seed'].astype('Int64')
    assert not q.duplicated(['arm', 'seed', 'query']).any()
    assert (q['query'].str.len() == q['length']).all()
    g = q.groupby(['arm', 'seed', 'length'], dropna=False).agg(n=('query', 'size'), m=('matched', 'sum'), rows=('retained_rows', 'sum')).reset_index()
    j = t.merge(g, on=['arm', 'seed', 'length'], how='left')
    j[['n', 'm', 'rows']] = j[['n', 'm', 'rows']].fillna(0)
    assert (j.n == j.distinct_queries).all() and (j.m == j.matching_queries).all() and (j.rows == j.retained_rows).all()
    sim = t[t.arm.ne('REAL')].copy()
    assert sorted(int(x) for x in sim.seed.unique()) == sorted(SCRAMBLE_SEEDS)
    assert set(t.length) == set(range(MIN_LOOKUP_LENGTH, MAX_LOOKUP_LENGTH + 1))
    frac = lambda L: (sim[sim.length.eq(L)].matching_queries / sim[sim.length.eq(L)].distinct_queries)
    long_ = sim[sim.length.ge(SHORT_LENGTH_CUTOFF)]
    real = t[t.arm.eq('REAL')].set_index('length')
    den = lambda L: ';'.join(str(int(x)) for x in sim[sim.length.eq(L)].sort_values('seed').distinct_queries)
    return {
        'n_scrambled_seeds': int(sim.seed.nunique()),
        'n_scrambled_seeds_word': NUMBER_WORDS[int(sim.seed.nunique())],
        'authentic_len9_denominator': int(t[t.arm.eq('REAL') & t.length.eq(9)].distinct_queries.iloc[0]),
        'scrambled_len9_denominators': den(9), 'scrambled_len10_denominators': den(10),
        'scrambled_len9_min_pct': float(frac(9).min() * 100), 'scrambled_len9_max_pct': float(frac(9).max() * 100),
        'scrambled_len10_min_pct': float(frac(10).min() * 100), 'scrambled_len10_max_pct': float(frac(10).max() * 100),
        'scrambled_ge13_query_seed_evaluations': int(long_.distinct_queries.sum()),
        'scrambled_ge13_matches': int(long_.matching_queries.sum()),
        'authentic_len9_pct': float(real.loc[9, 'matching_queries'] / real.loc[9, 'distinct_queries'] * 100),
        'authentic_distinct_queries_total': int(real.distinct_queries.sum()),
    }


def stat_cd_grid(sd):
    """Selection-sensitivity medians and exclusive quartiles from all 1,102 observations.

    Unit: acquisition; n = 58 acquisitions per arm, the same acquisitions in every arm
    (repeated measures). 19 arms: 6 minimum lengths at top 30%, 12 score fractions at
    minimum length nine, and the older 97-pool reference at nine/30%. Quartiles use
    ``statistics.quantiles(n=4)`` (method 'exclusive'), as the upstream sweep did.

    Returns
    -------
    dict
    """
    obs = load(sd, 'selection_grid_acquisitions.tsv', exact=True)
    a = load(sd, 'selection_grid_summary.tsv', exact=True).set_index('arm')
    worst = 0.0
    for arm, g in obs.groupby('arm'):
        for c in ['peptides_q01', 'razor_selected', 'proteins_q01', 'decoy_frac_q01']:
            v = sorted(g[c].tolist())
            qq = st.quantiles(v, n=4, method=QUARTILE_METHOD)
            for stat_, val in [('median', st.median(v)), ('q1', qq[0]), ('q3', qq[2])]:
                ref = a.loc[arm, f'{c}_{stat_}']
                worst = max(worst, abs(val - ref) / max(1.0, abs(ref)))
    assert worst < 1e-12, worst
    score = sorted(a[a.kind.eq('score')].top_percent.round(2))
    length = sorted(a[a.kind.eq('length')].min_length.tolist() + [9])
    assert tuple(score) == SCORE_GRID and tuple(length) == LENGTH_GRID
    assert set(obs.groupby('arm').size()) == {58} and obs['sample'].nunique() == 58
    tp = a.loc['tp0.30']
    return {
        'observations': int(len(obs)), 'arms': int(obs.arm.nunique()), 'acquisitions': int(obs['sample'].nunique()),
        'distinct_settings': int(len({(r.min_length, round(r.top_percent, 2)) for r in a.itertuples()})),
        'max_relative_difference_vs_summary': worst,
        'production_peptides_median': float(tp.peptides_q01_median),
        'production_peptides_q1': float(tp.peptides_q01_q1), 'production_peptides_q3': float(tp.peptides_q01_q3),
        'production_selected_proteins_median': float(tp.razor_selected_median),
        'production_selected_proteins_q1': float(tp.razor_selected_q1), 'production_selected_proteins_q3': float(tp.razor_selected_q3),
    }


def stat_e_bags(sd):
    """Gains, losses, nets and per-acquisition positives of the 48-acquisition comparison.

    Unit: peptide-acquisition identification (canonical peptide, deduplicated within an
    acquisition, nominal peptide q <= 0.01); n = 48 acquisitions (15 CAMPI, 32 healthy
    stool, 1 CapScan) per arm. Coverage of the eight-ordering bound is counted over
    lookup-eligible prediction rows of the same 48 acquisitions.

    Returns
    -------
    dict
    """
    f = load(sd, 'bags_bounded48_per_acquisition.tsv')
    s = bags_summary(sd)
    acq = f[['cohort', 'sample']].drop_duplicates()
    pos = f[f.net > 0].groupby('arm').size()
    hs = f[f.cohort.eq('HSTOOL')]
    hs_mass = hs[hs.arm.eq('mass_all')]
    o = load(sd, 'bags_bounded48_orderings.tsv')
    elig = o[['cohort', 'sample', 'lookup_eligible_predictions']].drop_duplicates()
    assert len(elig) == 48
    assert (o.groupby(['cohort', 'sample']).predictions.sum().values == elig.set_index(['cohort', 'sample']).sort_index().lookup_eligible_predictions.values).all()
    le8 = int(o[o.orderings.le(8)].predictions.sum()); tot = int(elig.lookup_eligible_predictions.sum())
    out = {'acquisitions': int(len(acq)), 'arm_rows': int(len(f)), 'arms': int(f.arm.nunique()),
           'arms_word': NUMBER_WORDS[int(f.arm.nunique())]}
    out.update({f'cohort_{c}': int(n) for c, n in acq.cohort.value_counts().items()})
    out['cohort_CAPSCAN_word'] = NUMBER_WORDS[out['cohort_CAPSCAN']]
    for arm in s.index:
        out.update({f'{arm}_gained': int(s.loc[arm, 'gained']), f'{arm}_lost': int(s.loc[arm, 'lost']),
                    f'{arm}_net': int(s.loc[arm, 'net']), f'{arm}_positive_acquisitions': int(pos.get(arm, 0))})
    out.update({'hstool_acquisitions': int(hs_mass['sample'].nunique()), 'hstool_mass_net': int(hs_mass.net.sum()),
                'hstool_mass_positive': int((hs_mass.net > 0).sum()),
                'mass_all_net_loss': -int(s.loc['mass_all', 'net']),
                'le8_eligible_rows': le8, 'eligible_rows': tot, 'le8_eligible_pct': le8 / tot * 100})
    return out


def stat_f_fragment(sd):
    """Six-acquisition fragment audit totals.

    Unit: same-spectrum comparison pair (gained bag-arm sequence against the exact-arm
    rank-one sequence); n = 215 pairs in six acquisitions.

    Returns
    -------
    dict
    """
    fr = fragment_subset(sd)
    assert (fr.zero_diagnostic + fr.positive_diagnostic == fr.spectrum_events).all()
    return {'acquisitions': int(len(fr)), 'acquisitions_word': NUMBER_WORDS[int(len(fr))], 'pairs': int(fr.spectrum_events.sum()),
            'zero_diagnostic': int(fr.zero_diagnostic.sum()), 'positive_diagnostic': int(fr.positive_diagnostic.sum()),
            **{f'cohort_{c}': int(n) for c, n in fr.cohort.value_counts().items()}}


# =====================================================================================
# v30 house-style rebuild: shared helpers
# Moved unchanged from the v30 rebuild helper module (restyle_common.py, used for
# Supplementary Figures 1-3; see PROVENANCE.md). Only the lab paths and environment
# setup at the top of that module were dropped; fl_style is imported from the
# repository's style package instead. Functions:
#   fingerprint()   plotted arrays + axis limits/scales of a list of axes (values-unchanged proof)
#   compare()       exact comparison of two fingerprints
#   texts()         every visible string (texts, tick labels, legend entries) of a figure
#   capture_saved() run a producer and capture the Figure objects it saves (nothing written)
#   align_letters() put the bold panel letters of one figure column on one x position
#   draw_rel_label() draw a label containing <= or >= in Söhne (the face has no such glyphs)
#   render_png()    rasterise page 1 of a PDF at 300 dpi
# =====================================================================================
import hashlib  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import Collection, LineCollection, PathCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402


def sha(p):
    """Return the hex sha256 of a file."""
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ---------------------------------------------------------------- values-unchanged proof
def _arrays(ax):
    out = []
    for ln in ax.get_lines():
        xy = np.asarray(ln.get_xydata(), dtype=float)
        if xy.size:
            out.append(('line', xy))
    for c in ax.collections:
        if isinstance(c, PathCollection):
            off = np.asarray(c.get_offsets(), dtype=float)
            if off.size:
                out.append(('points', off))
        elif isinstance(c, LineCollection):
            segs = c.get_segments()
            if len(segs):
                out.append(('segments', np.concatenate([np.asarray(s, dtype=float) for s in segs])))
    for p in ax.patches:
        if isinstance(p, Rectangle) and not getattr(p, '_restyle_decor', False):
            out.append(('bar', np.array([p.get_x(), p.get_y(), p.get_width(), p.get_height()], dtype=float)))
    return out


def fingerprint(axes):
    """Per axes: sorted list of (kind, array) plus x/y limits and scales."""
    fp = []
    for ax in axes:
        arr = sorted(_arrays(ax), key=lambda t: (t[0], t[1].shape, t[1].tobytes()))
        fp.append(dict(arrays=arr, xlim=tuple(map(float, ax.get_xlim())), ylim=tuple(map(float, ax.get_ylim())),
                       xscale=ax.get_xscale(), yscale=ax.get_yscale()))
    return fp


def compare(old, new, names):
    """Exact equality of every plotted array, limit and scale; returns a JSON-able report."""
    rep = []
    assert len(old) == len(new), (len(old), len(new))
    for nm, o, n in zip(names, old, new):
        r = dict(panel=nm, n_arrays_old=len(o['arrays']), n_arrays_new=len(n['arrays']),
                 n_values=int(sum(a.size for _, a in o['arrays'])))
        r['arrays_identical'] = (len(o['arrays']) == len(n['arrays']) and all(
            ko == kn and ao.shape == an.shape and np.array_equal(ao, an, equal_nan=True) for (ko, ao), (kn, an) in zip(o['arrays'], n['arrays'])))
        r['xlim_identical'] = o['xlim'] == n['xlim']
        r['ylim_identical'] = o['ylim'] == n['ylim']
        r['scales_identical'] = (o['xscale'], o['yscale']) == (n['xscale'], n['yscale'])
        r['xlim'], r['ylim'] = o['xlim'], o['ylim']
        rep.append(r)
    ok = all(r['arrays_identical'] and r['xlim_identical'] and r['ylim_identical'] and r['scales_identical'] for r in rep)
    return ok, rep


def texts(fig):
    """Visible strings: free texts, legend entries, and tick labels whose tick lies inside
    the axis view interval (matplotlib keeps unseen ticks beyond the limits)."""
    fig.canvas.draw()
    tick_ids, keep = set(), []
    axes = list(fig.axes)
    for ax in list(axes):
        axes.extend(a for a in ax.child_axes if a not in axes)
    for ax in axes:
        for axis, i in ((ax.xaxis, 0), (ax.yaxis, 1)):
            lo, hi = sorted(axis.get_view_interval())
            span = hi - lo
            for t in axis.get_ticklabels(which='both'):
                tick_ids.add(id(t))
                v = t.get_position()[i]
                if t.get_visible() and t.get_text().strip() and lo - 1e-9 * abs(span) <= v <= hi + 1e-9 * abs(span):
                    keep.append(t.get_text().strip())
    for t in fig.findobj(match=lambda o: hasattr(o, 'get_text') and hasattr(o, 'get_fontsize')):
        if id(t) in tick_ids:
            continue
        try:
            if t.get_visible() and t.get_text().strip():
                keep.append(t.get_text().strip())
        except Exception:  # noqa: BLE001
            pass
    return keep


class capture_saved:
    """Context manager: Figure.savefig records the figure instead of writing."""

    def __enter__(self):
        from matplotlib.figure import Figure
        self.figs = []
        self._orig = Figure.savefig
        figs = self.figs

        self.texts = []
        texts_at_save = self.texts

        def fake(fig, *a, **k):
            figs.append(fig)
            texts_at_save.append(texts(fig))   # strings as drawn under the producer's own rcParams
        Figure.savefig = fake
        return self

    def __exit__(self, *exc):
        from matplotlib.figure import Figure
        Figure.savefig = self._orig
        return False


# ---------------------------------------------------------------- layout helpers
def align_letters(fig, columns):
    """columns: list of lists of axes. Letters (left titles) of a column share the x of
    the left-most drawn element of that column (tick labels included)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    inv = fig.transFigure.inverted()
    for col in columns:
        x0 = min(inv.transform((ax.get_tightbbox(r, for_layout_only=False).x0, 0))[0] for ax in col)
        for ax in col:
            # remove the letter itself from the bbox calculation: use x of figure -> axes
            ax_x0 = ax.get_position().x0
            ax_w = ax.get_position().width
            ax._left_title.set_x((x0 - ax_x0) / ax_w)
    fig.canvas.draw()


def draw_rel_label(ax, x, y, prefix, rel, suffix, transform, fontsize, color, ha='right', va='center'):
    """Draw prefix + (<= or >=) + suffix in Söhne. Söhne Buch has no U+2264/U+2265, so the
    relation is the Söhne '<' or '>' glyph (as text) raised slightly, with a rule beneath it
    of the type's stem weight. Returns the artists. ha in {'right', 'left'}."""
    from matplotlib.transforms import offset_copy
    fig = ax.figure
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    k = 72.0 / fig.dpi

    def width(s):
        t = ax.text(0, 0, s, fontsize=fontsize)
        w = t.get_window_extent(r).width * k
        t.remove()
        return w
    sym = '<' if rel == 'le' else '>'
    w_pre, w_sym, w_suf = width(prefix) if prefix else 0.0, width(sym), width(suffix)
    total = w_pre + w_sym + w_suf
    start = -total if ha == 'right' else 0.0
    arts = []
    if prefix:
        arts.append(ax.text(x, y, prefix, transform=offset_copy(transform, fig=fig, x=start, y=0, units='points'),
                            ha='left', va=va, fontsize=fontsize, color=color, clip_on=False))
    raise_pt = 0.10 * fontsize
    arts.append(ax.text(x, y, sym, transform=offset_copy(transform, fig=fig, x=start + w_pre, y=raise_pt, units='points'),
                        ha='left', va=va, fontsize=fontsize, color=color, clip_on=False))
    arts.append(ax.text(x, y, suffix, transform=offset_copy(transform, fig=fig, x=start + w_pre + w_sym, y=0, units='points'),
                        ha='left', va=va, fontsize=fontsize, color=color, clip_on=False))
    # rule under the raised chevron, inset like the glyph's side bearings; positioned in points
    yb = -0.30 * fontsize + raise_pt
    x_a, x_b = start + w_pre + 0.14 * w_sym, start + w_pre + 0.86 * w_sym
    from matplotlib.transforms import Affine2D
    # a filled bar, not a stroked line: fl.save raises strokes below 1 pt to 1 pt, which would make
    # the rule twice the weight of the Söhne Buch stems; a fill with zero edge is left alone.
    th = 0.075 * fontsize
    rule = Rectangle((x_a, yb - th / 2), x_b - x_a, th, facecolor=color, edgecolor='none', linewidth=0, clip_on=False,
                     transform=Affine2D().scale(1 / 72.0) + fig.dpi_scale_trans + _point_anchor(x, y, transform))
    rule._restyle_decor = True   # typographic, not data: excluded from the fingerprint
    ax.add_artist(rule)
    arts.append(rule)
    return arts


def _point_anchor(x, y, transform):
    """Transform that translates display units so (0, 0) lands at (x, y) in `transform`."""
    from matplotlib.transforms import ScaledTranslation
    return ScaledTranslation(x, y, transform)


# ---------------------------------------------------------------- raster
def render_png(pdf, png, dpi=300):
    import fitz
    pix = fitz.open(pdf)[0].get_pixmap(dpi=dpi, alpha=False)
    pix.set_dpi(dpi, dpi)
    pix.save(png)
    return pix.width, pix.height




# =====================================================================================
# v30: house-style rebuild of both canvases (moved unchanged from build_sfig02_style.py)
# =====================================================================================
V30_W_ABCD, V30_H_ABCD = 160.0, 135.2   # v30 build: placed size of a-d (v29 docx extent 160.0 x 135.2 mm)
V30_W_EF, V30_H_EF = 155.0, 117.1       # v30 build: placed size of e,f (v29 docx extent 154.9 x 117.1 mm)
V30_MARKER_PT = 3.2                     # v30 build: marker size; the producer used 6 at a larger canvas
V30_PNG_DPI = 300                       # v30 build: PNG rasterised at 300 dpi from the PDF
V30_TITLES_ABCD = {'a  Exact reservoir lookup', 'b  Length-specific denominator', 'c  Length sensitivity',
                   'd  Score-fraction sensitivity', 'Scrambled: 9–13 residues',
                   'Evidence selection and exact-lookup specificity'}   # removed in v30; content is in the legend
V30_TITLES_EF = {'e  Completed matched comparison: 48 acquisitions', 'f  Separate fragment audit: six acquisitions'}
V30_COMPOSED = {'Bags: ≤2 orderings', 'Bags: ≤4 orderings', 'Bags: ≤8 orderings', '≥1 diagnostic ion'}


def panel_abcd_v30(sd):
    """Draw panels a-d in the house style (v30): same data and encodings as ``panel_abcd``.

    Changes from the producer (layout only): Söhne type, ink, Okabe-Ito colours by meaning,
    placed width 160 mm, subpanel titles removed, one c/d key in d's lower right, inset raised.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig02``.

    Returns
    -------
    tuple
        ``(fig, axes)`` with ``axes`` = [a, b, c, d, a inset].
    """
    from style import fl_style as fl
    mpl.rcdefaults()
    fl.init()
    MS = V30_MARKER_PT
    q = load(sd, 'selection_length_lookup.tsv')
    a = load(sd, 'selection_grid_summary.tsv')
    matched = a[a.kind.eq('length') | a.arm.eq('tp0.30')].sort_values('min_length').copy()
    fig, axs = plt.subplots(2, 2, figsize=fl.mm(V30_W_ABCD, V30_H_ABCD)); ax = axs.ravel()
    fig.subplots_adjust(left=0.085, right=0.985, top=0.885, bottom=0.075, wspace=0.30, hspace=0.42)
    real = q[q.arm.eq('REAL')]
    ax[0].plot(real.length, real.fraction_matched * 100, 'o-', c=fl.FL, ms=MS, label='Authentic')
    for _, g in q[~q.arm.eq('REAL')].groupby('arm'):
        ax[0].plot(g.length, g.fraction_matched * 100, c=fl.ORANGE, alpha=.6)
    ax[0].plot([], [], c=fl.ORANGE, label='Five scrambled seeds')
    ax[0].set(xlabel='Peptide length (residues)', ylabel='Distinct queries matching (%)', xlim=(8.5, 50.5))
    fig.legend(*ax[0].get_legend_handles_labels(), loc='upper center', bbox_to_anchor=(.5, .995), ncol=2, frameon=False, fontsize=7)
    inset = ax[0].inset_axes([.45, .21, .5, .20])   # raised 0.06 so its tick labels clear the scrambled trace at 0 %
    for _, g in q[~q.arm.eq('REAL')].groupby('arm'):
        inset.plot(g.length, g.fraction_matched * 100, 'o-', ms=1.6, c=fl.ORANGE, alpha=.6)
    inset.set(xlim=(8.8, 13.2), ylim=(-.3, 10), xticks=[9, 11, 13], yticks=[0, 5, 10])
    ax[1].plot(real.length, real.distinct_queries, 'o-', c=fl.FL, ms=MS)
    ax[1].set(xlabel='Peptide length (residues)', ylabel='Distinct authentic queries', xlim=(8.5, 50.5))
    for axis, kind, x in [(ax[2], 'length', 'min_length'), (ax[3], 'score', 'top_percent')]:
        g = (matched if kind == 'length' else a[a.kind.eq(kind)]).sort_values(x)
        xx = g[x].to_numpy() * (100 if kind == 'score' else 1)
        axis.errorbar(xx, g.peptides_q01_median / 1000,
                      yerr=np.array([g.peptides_q01_median - g.peptides_q01_q1, g.peptides_q01_q3 - g.peptides_q01_median]) / 1000,
                      fmt='o-', c=fl.FL, ms=MS, elinewidth=1, label='Accepted peptides')
        axis.plot(xx, g.razor_selected_median / 1000, 's--', c=fl.ORANGE, ms=MS, label='Selected proteins (median)')
        axis.set(xlabel='Minimum lookup length (residues)' if kind == 'length' else 'Top-scoring predictions (%)',
                 ylabel='Count (thousands)', ylim=(0, 48))
        if kind == 'score':
            # one key for c and d (identical encoding), in d's empty lower right
            axis.legend(frameon=False, fontsize=7, loc='lower right', borderaxespad=0.3)
        axis.axvline(9 if kind == 'length' else 30, c=fl.GREY, lw=1, ls=':')
    ax[0].set_xticks(range(10, 51, 5)); ax[1].set_xticks(range(10, 51, 5))   # v29 tick spacing
    for axx in list(ax) + [inset]:
        fl.style(axx)
    inset.tick_params(labelsize=6.5)
    for axx, l in zip(ax, 'abcd'):
        fl.panel(axx, l)
    fl.normalise_type(fig, 'pre-align', verbose=False)   # letters at their saved size before aligning
    align_letters(fig, [[ax[0], ax[2]], [ax[1], ax[3]]])
    return fig, list(ax) + [inset]


def panel_ef_v30(sd):
    """Draw panels e,f (continuation) in the house style (v30): same data as ``panel_ef``.

    The net values become right-hand tick labels of a twin axis sharing e's y range; the
    ≤/≥ labels are composed from Söhne glyphs (``draw_rel_label``) because Söhne Buch has
    no U+2264/U+2265; f's key is drawn by hand below the axis.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig02``.

    Returns
    -------
    tuple
        ``(fig, [e, f], net_axis, net_strings)``.
    """
    from matplotlib.transforms import blended_transform_factory
    from style import fl_style as fl
    summary = bags_summary(sd)
    fr = fragment_subset(sd)
    assert fr.spectrum_events.sum() == 215 and fr.zero_diagnostic.sum() == 193
    fig2, axes2 = plt.subplots(2, 1, figsize=fl.mm(V30_W_EF, V30_H_EF), gridspec_kw={'height_ratios': [4, 1.3]})
    fig2.subplots_adjust(left=0.205, right=0.915, top=0.955, bottom=0.135, hspace=0.42)
    e = axes2[0]
    order = ['bags_le2', 'bags_le4', 'bags_le8', 'bags_all', 'mass_all', 'comb_all']
    labels = ['Bags: ≤2 orderings', 'Bags: ≤4 orderings', 'Bags: ≤8 orderings', 'Bags: unrestricted', 'Mass substitutions', 'Bags + mass']
    ss = summary.loc[order]; yy = np.arange(len(order))
    e.barh(yy, ss.gained, color=fl.FL, label='Gained')
    e.barh(yy, -ss.lost, color=fl.VERM, label='Lost')
    e.axvline(0, color=fl.INK, lw=1)
    e.set_yticks(yy, [lab if '≤' not in lab else '' for lab in labels]); e.invert_yaxis(); e.set_xlim(-42000, 48000)
    for i, (_, r) in enumerate(ss.iterrows()):
        e.text(-r.lost - 600, i, f'{int(r.lost):,}', ha='right', va='center', fontsize=6.5, color=fl.INK)
        e.text(r.gained + 600, i, f'{int(r.gained):,}', va='center', fontsize=6.5, color=fl.INK)
    e.set_xticks([-30000, 0, 30000], ['30,000 lost', '0', '30,000 gained'])
    e.set_xlabel('Peptide–acquisition identifications relative to exact lookup')
    e.legend(loc='upper left', frameon=False, fontsize=7)
    fl.style(e)
    # net difference: right-hand tick labels on a twin axis that shares e's y range exactly
    net = e.twinx(); net.set_ylim(e.get_ylim()); net.set_yticks(yy, [f'{int(r.net):+,}' for _, r in ss.iterrows()])
    for sp in ('top', 'right', 'left', 'bottom'):
        net.spines[sp].set_visible(False)
    net.tick_params(axis='y', length=0, labelsize=6.5, pad=2); net.set_ylabel('Net', rotation=270, labelpad=9)
    f = axes2[1]
    f.barh([0], [193], color=fl.GREY, label='No diagnostic ion')
    f.barh([0], [22], left=[193], color=fl.GREEN, label='≥1 diagnostic ion')
    f.text(96.5, 0, '193', ha='center', va='center', fontsize=7, color='white')
    f.text(204, 0, '22', ha='center', va='center', fontsize=7, color='white')
    f.set_xlim(0, 215); f.set_yticks([]); f.set_xticks([0, 50, 100, 150, 200, 215])
    f.set_xlabel('Order-related same-spectrum comparison pairs')
    fl.style(f); f.spines['left'].set_visible(False)
    fl.panel(e, 'e'); fl.panel(f, 'f')
    fl.normalise_type(fig2, 'pre-align', verbose=False)
    # <= labels of e (Söhne has no U+2264): composed from Söhne glyphs at the tick-label position
    fig2.canvas.draw()
    tick_tr = blended_transform_factory(e.transAxes, e.transData)
    pad_ax = (e.yaxis.get_major_ticks()[0].get_pad() + 2.6) / 72 / (e.get_position().width * fig2.get_figwidth())
    for i, lab in enumerate(labels):
        if '≤' in lab:
            pre, suf = lab.split('≤')
            draw_rel_label(e, -pad_ax, i, pre, 'le', suf, tick_tr, fontsize=7, color=fl.INK, ha='right')
    # legend for f, hand-drawn below the axis so the >= entry can be composed the same way
    lt = f.transAxes
    for x0, col, txt in [(0.27, fl.GREY, 'No diagnostic ion'), (0.56, fl.GREEN, None)]:
        key = Rectangle((x0, -0.73), 0.035, 0.14, transform=lt, facecolor=col, edgecolor='none', clip_on=False)
        key._restyle_decor = True   # legend key, not data: excluded from the fingerprint
        f.add_artist(key)
        if txt:
            f.text(x0 + 0.05, -0.66, txt, transform=lt, ha='left', va='center', fontsize=7, color=fl.INK)
        else:
            draw_rel_label(f, x0 + 0.05, -0.66, '', 'ge', '1 diagnostic ion', lt, fontsize=7, color=fl.INK, ha='left')
    align_letters(fig2, [[e, f]])
    assert net.get_ylim() == e.get_ylim()
    return fig2, [e, f], net, [f'{int(r.net):+,}' for _, r in ss.iterrows()]


def build_v30(sd, out_dir):
    """Build both v30 canvases and prove the plotted values unchanged.

    Steps of ``build_sfig02_style.py``: (1) run ``panel_abcd`` and ``panel_ef`` with
    ``savefig`` intercepted and fingerprint their axes; (2) draw ``panel_abcd_v30`` and
    ``panel_ef_v30``; (3) assert every plotted array, limit and scale identical and that only
    the removed titles and the composed ≤/≥ labels are lost as strings; (4) write
    ``SFig02_v30.pdf/.png`` and ``SFig02_continued_v30.pdf/.png`` (300 dpi).

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig02``.
    out_dir : pathlib.Path
        Output directory.

    Returns
    -------
    dict
        ``values_identical_abcd``, ``values_identical_ef``, per-panel reports, lost and
        added label strings, and the PNG pixel sizes.
    """
    from style import fl_style as fl
    out_dir = Path(out_dir)
    mpl.rcdefaults()
    with capture_saved() as cap:
        panel_abcd(sd, out_dir / '_unused_abcd.pdf')
        panel_ef(sd, out_dir / '_unused_ef.pdf')
    old_abcd, old_ef = cap.figs
    oa = list(old_abcd.axes[:4]); o_inset = oa[0].child_axes[0]
    OLD_ABCD = fingerprint(oa + [o_inset]); OLD_EF = fingerprint(list(old_ef.axes[:2]))
    old_text_abcd, old_text_ef = cap.texts
    fig, axes = panel_abcd_v30(sd)
    NEW_ABCD = fingerprint(axes)
    new_text_abcd = texts(fig)
    fig2, ef, net, net_strings = panel_ef_v30(sd)
    NEW_EF = fingerprint(ef)
    new_text_ef = texts(fig2)
    ok1, rep1 = compare(OLD_ABCD, NEW_ABCD, ['a', 'b', 'c', 'd', 'a inset'])
    ok2, rep2 = compare(OLD_EF, NEW_EF, ['e', 'f'])
    lost1 = sorted(set(old_text_abcd) - set(new_text_abcd)); added1 = sorted(set(new_text_abcd) - set(old_text_abcd))
    lost2 = sorted(set(old_text_ef) - set(new_text_ef)); added2 = sorted(set(new_text_ef) - set(old_text_ef))
    assert ok1 and ok2, (rep1, rep2)
    assert set(lost1) <= V30_TITLES_ABCD, lost1
    assert set(lost2) <= V30_TITLES_EF | V30_COMPOSED, lost2
    for lab in V30_COMPOSED:
        rel = '≤' if '≤' in lab else '≥'
        pre, suf = lab.split(rel)
        assert (not pre.strip() or pre.strip() in new_text_ef) and suf.strip() in new_text_ef, lab
    assert all(s in new_text_ef for s in net_strings)
    px = {}
    for fg, stem in [(fig, 'SFig02_v30'), (fig2, 'SFig02_continued_v30')]:
        pdf, _ = fl.save(fg, out_dir / stem, formats=('pdf',))
        px[stem] = render_png(pdf, out_dir / f'{stem}.png', V30_PNG_DPI)
        plt.close(fg)
    return dict(values_identical_abcd=ok1, values_identical_ef=ok2, panels=rep1 + rep2,
                label_strings_lost_abcd=lost1, label_strings_added_abcd=added1,
                label_strings_lost_ef=lost2, label_strings_added_ef=added2, png_px=px)
