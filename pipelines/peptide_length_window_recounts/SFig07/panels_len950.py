"""Supplementary Figure 7: common search databases (depth, identity retention, precision).

Panels a-d are drawn exactly as by the original producer (``build_supp_fig7.py``, see
PROVENANCE.md). The plotting code below is that producer's code with file paths made
repository-relative and the code split into one function per panel; the drawing calls,
their order, the palette, the in-figure titles and the layout are unchanged.

Source data (``source_data/SFig07/``):

- ``identity_comparison.tsv``: per-acquisition identity retention, study-level arms (panel a)
- ``identity_per_acquisition.tsv``: per-acquisition identity retention, all arms including
  per-acquisition confidence (panel a, text line 115)
- ``identity_summary.tsv``: per-arm summary written by the upstream comparison (line 115)
- ``completeness_all.tsv``: features complete in all replicates (line 115)
- ``common_search_participants.tsv``: PREDICT and CapScan participant summaries, pseudonymised
  identifiers (panel b)
- ``microbial_group_pairs.tsv``: per acquisition-pair microbial protein-group CV and
  correlation (panels c, d, text line 769)
- ``microbial_group_pair_summary.tsv``: upstream median-over-pairs summary (line 769)
- ``common_search_pair_CV.tsv``: earlier pair-CV table used by the producer as a consistency check
- ``technical_agreement_summary.tsv``: leading-species agreement (line 769)
- ``arm_definitions.tsv``: meaning of each arm label
"""

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import font_manager  # noqa: E402

# ---------------------------------------------------------------------------------------
# Choices (named constants; README.md Choices table)
# ---------------------------------------------------------------------------------------
MIN_LENGTH_PANEL_A = 9          # len950: panel a uses 9-50-residue peptides (accepted 5-50 lists filtered by length); repo value 5
PEPTIDE_Q_CUTOFF = 0.01         # rank-one target peptides at q <= 0.01; v28d line 765 (applied upstream)
N_ACQUISITIONS = 58             # healthy-stool acquisitions; v28d lines 113, 115
TIERS = ("TechnicalReplicates", "WorkflowReplicates")   # panels c, d; v28d line 437
TIER_PAIRS = {"TechnicalReplicates": 66, "WorkflowReplicates": 105}  # v28d line 437
PAIR_CV_DDOF = 0                # population SD / mean for pair CV; v28d lines 437, 769 (applied upstream)
PAIR_SUMMARY = "median"         # median across rows within a pair, then median across pairs; v28d line 769
MICROBIAL_CLASS = "Microbial"   # protein_class used for the producer's consistency check
FIG_SIZE_IN = (7.1, 6.4)        # producer figure size (inches)
PNG_DPI = 300                   # producer PNG resolution
FONT_FILE = "soehne-buch.ttf"   # producer font file (Söhne Buch, TrueType)
FONT_FILE_FALLBACK = "soehne-buch.otf"
COL = ["#2878a3", "#df8f44", "#619b68", "#b35872", "#8277ac"]  # producer palette (legacy, kept unchanged)
RC = {"font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8, "xtick.labelsize": 6.5,
      "ytick.labelsize": 7, "pdf.fonttype": 42, "svg.fonttype": "none",
      "axes.spines.top": False, "axes.spines.right": False}   # producer rcParams


# ---------------------------------------------------------------------------------------
# Style
# ---------------------------------------------------------------------------------------
def producer_rc():
    """Return the producer's rcParams, with the Söhne family when the fonts are available.

    The producer started from matplotlib defaults and applied ``RC`` plus the family of
    ``soehne-buch.ttf``. The font directory is read from ``FASTALAKE_FONT_DIR``. Without it,
    the matplotlib default family is kept and text pixels differ from the published figure.

    Returns
    -------
    tuple of (dict, str or None)
        rcParams to apply on top of ``matplotlib.rcdefaults()``, and the font file name used.
    """
    rc = dict(RC)
    d = os.environ.get("FASTALAKE_FONT_DIR") or os.environ.get("FL_STYLE_DIR") or ""
    used = None
    for name in (FONT_FILE, FONT_FILE_FALLBACK):
        f = Path(d) / name
        if d and f.is_file():
            font_manager.fontManager.addfont(str(f))
            rc["font.family"] = font_manager.FontProperties(fname=str(f)).get_name()
            used = name
            break
    return rc, used


# ---------------------------------------------------------------------------------------
# Producer helpers (unchanged)
# ---------------------------------------------------------------------------------------
def dots(ax, groups, ylabel, percent=False):
    """Draw jittered dot columns with a black median bar (producer function, unchanged).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    groups : list of (label, values, colour)
    ylabel : str
    percent : bool
        Multiply values by 100.

    Returns
    -------
    list of tuple
        ``(label, n, median)`` for each column.
    """
    rec = []
    for i, (label, values, c) in enumerate(groups):
        v = np.asarray(values, dtype=float)
        v = v[np.isfinite(v)] * (100 if percent else 1)
        x = i + np.linspace(-.16, .16, len(v))
        ax.scatter(x, v, s=8, alpha=.48, color=c, rasterized=len(v) > 100)
        if len(v):
            ax.hlines(np.median(v), i - .25, i + .25, color='black', lw=1.5)
        rec.append((label.replace('\n', ' '), len(v), float(np.median(v))))
    ax.set_xticks(range(len(groups)), [g[0] for g in groups])
    ax.set_ylabel(ylabel)
    return rec


def title(ax, l, t):
    """Write the producer's in-figure panel title (kept unchanged, see OPEN_ITEMS.md)."""
    ax.set_title(l + '  ' + t, loc='left', pad=9)


def load(sd):
    """Read every table the figure needs.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig07``.

    Returns
    -------
    dict of pandas.DataFrame
    """
    r = lambda n: pd.read_csv(Path(sd) / n, sep='\t')  # noqa: E731
    return dict(ident=r('identity_comparison.tsv'), mine=r('identity_per_acquisition.tsv'),
                pat=r('common_search_participants.tsv'), pairs=r('microbial_group_pairs.tsv'),
                src=r('common_search_pair_CV.tsv'), summary=r('identity_summary.tsv'),
                completeness=r('completeness_all.tsv'), pair_summary=r('microbial_group_pair_summary.tsv'),
                agreement=r('technical_agreement_summary.tsv'))


# ---------------------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------------------
def panel_a(ax, ident, mine, plotted, checks):
    """Panel a: original identities recovered per acquisition, five arms (58 acquisitions each).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    ident : pandas.DataFrame
        ``identity_comparison.tsv`` (study-level arms).
    mine : pandas.DataFrame
        ``identity_per_acquisition.tsv`` (adds the per-acquisition-confidence arm).
    plotted, checks : list
        Appended with plotted medians and the producer's consistency check.

    Returns
    -------
    list of tuple
        ``(label, n, median)`` per column; unit acquisition, n = 58 per column.
    """
    ident = ident[ident.min_length.eq(MIN_LENGTH_PANEL_A)]
    fq = mine[mine.arm.eq('joint_fileq')]
    a = ident[ident.test.eq('joint_corrected')].set_index('sample').original_identity_recovery
    b = mine[mine.arm.eq('joint_corrected')].set_index('sample').original_identity_recovery
    checks.append(dict(check='panel a joint_corrected recovery equals A21',
                       passed=bool(np.allclose(a.sort_index(), b.sort_index()))))
    A = [('joint corrected\nstudy q', ident[ident.test.eq('joint_corrected')].original_identity_recovery, COL[0]),
         ('joint original\nstudy q', ident[ident.test.eq('joint_original')].original_identity_recovery, COL[1]),
         ('joint corrected\nacquisition q', fq.original_identity_recovery, COL[4]),
         ('repeat\nindividual', ident[ident.test.eq('repeat_individual')].original_identity_recovery, COL[2]),
         ('union\nindividual', ident[ident.test.eq('union_individual')].original_identity_recovery, COL[3])]
    rec = dots(ax, A, 'Original identities recovered (%)', True)
    for r in rec:
        plotted.append(dict(panel='a', column=r[0], n=r[1], median=r[2], unit='percent of original identities'))
    title(ax, 'a', 'Healthy stool, 58 acquisitions')
    return rec


def panel_b(ax, pat, plotted):
    """Panel b: participant-level percentage change, common FASTA versus own FASTA.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    pat : pandas.DataFrame
        ``common_search_participants.tsv``; unit participant, n = 93 PREDICT and 15 CapScan.
    plotted : list

    Returns
    -------
    dict
        Participants per cohort.
    """
    n = {}
    for i, (co, g) in enumerate(pat.groupby('cohort')):
        ax.scatter(g.peptides_percent_difference, g.directlfq_groups_percent_difference, s=12, alpha=.6, color=COL[i], label=co)
        plotted.append(dict(panel='b', column=co, n=len(g), median=float(g.peptides_percent_difference.median()),
                            unit='median peptides percent difference (unchanged panel)'))
        n[co] = len(g)
    ax.axhline(0, c='grey', lw=.6)
    ax.axvline(0, c='grey', lw=.6)
    ax.set(xlabel='Accepted peptides: common vs own (%)', ylabel='DirectLFQ groups: common vs own (%)')
    ax.legend(frameon=False)
    title(ax, 'b', 'Clinical common-FASTA tradeoffs')
    return n


def check_cd(pairs, src, checks):
    """Producer check: study-level pair CVs equal the earlier source table (171 pairs each)."""
    src = src[src.protein_class.eq(MICROBIAL_CLASS)]
    for de_src, de in [('A_arm_alone', 'own_rows'), ('B_matched_both_arms', 'matched')]:
        for arm in ['original', 'joint_corrected']:
            s = src[src.definition.eq(de_src) & src.arm.eq(arm)].set_index(['sample_a', 'sample_b']).median_pair_cv_pct.sort_index()
            q = pairs[pairs.definition.eq(de) & pairs.method.eq(arm) & (pairs.comparison.eq(arm) if de == 'own_rows' else pairs.comparison.eq('original_vs_joint_corrected'))].set_index(['sample_a', 'sample_b']).median_pair_cv_pct.sort_index()
            checks.append(dict(check=f'panel c,d {de} {arm} pair CV equals existing source table',
                               passed=bool(len(s) == len(q) == 171 and np.allclose(s.values, q.values))))


def _sel(pairs, de, comp, method, t):
    return pairs[pairs.definition.eq(de) & pairs.comparison.eq(comp) & pairs.method.eq(method) & pairs.tier.eq(t)].median_pair_cv_pct


def panel_cd(ax, pairs, i, plotted):
    """Panels c (i = 0, technical) and d (i = 1, preparation): per-pair median microbial CV.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    pairs : pandas.DataFrame
        ``microbial_group_pairs.tsv``; unit acquisition pair (dependent), n = 66 or 105.
    i : int
    plotted : list

    Returns
    -------
    list of tuple
        ``(label, n, median)`` per column.
    """
    t = TIERS[i]
    G = [('Own\noriginal', _sel(pairs, 'own_rows', 'original', 'original', t), COL[0]),
         ('Own\njoint\nstudy q', _sel(pairs, 'own_rows', 'joint_corrected', 'joint_corrected', t), COL[1]),
         ('Own\njoint\nacq. q', _sel(pairs, 'own_rows', 'joint_fileq', 'joint_fileq', t), COL[4]),
         ('Matched\noriginal', _sel(pairs, 'matched', 'original_vs_joint_corrected', 'original', t), COL[0]),
         ('Matched\njoint\nstudy q', _sel(pairs, 'matched', 'original_vs_joint_corrected', 'joint_corrected', t), COL[1]),
         ('Matched\noriginal', _sel(pairs, 'matched', 'original_vs_joint_fileq', 'original', t), COL[0]),
         ('Matched\njoint\nacq. q', _sel(pairs, 'matched', 'original_vs_joint_fileq', 'joint_fileq', t), COL[4])]
    rec = dots(ax, G, 'Median pair CV (%)')
    for k, r in enumerate(rec):
        plotted.append(dict(panel='cd'[i], column=r[0] + (' (matched to study-q joint)' if k in (3, 4) else ' (matched to acquisition-q joint)' if k in (5, 6) else ''),
                            n=r[1], median=r[2], unit='median pair CV, percent'))
    ax.axvline(2.5, c='grey', lw=.5)
    ax.axvline(4.5, c='grey', lw=.5, ls=':')
    title(ax, 'cd'[i], ['Technical replicates (66 pairs)', 'Preparation replicates (105 pairs)'][i])
    return rec


def render(sd, out):
    """Draw the whole figure as the producer did and save PNG (300 dpi), PDF and SVG.

    The figure is drawn in one canvas (2 x 2 axes) because the producer drew it that way;
    the per-panel functions above draw into its axes in the producer's order.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig07``.
    out : pathlib.Path
        Output directory; writes ``SFig07.{png,pdf,svg}`` and ``plotted_values.tsv``.

    Returns
    -------
    dict
        ``plotted`` (DataFrame of plotted medians), ``checks`` (producer checks),
        ``font`` (font file used or None), ``axes_bbox_px`` (panel boxes in the 300-dpi PNG).
    """
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    T = load(sd)
    rc, font = producer_rc()
    plotted, checks = [], []
    with plt.rc_context():
        matplotlib.rcdefaults()
        plt.rcParams.update(rc)
        fig, axes = plt.subplots(2, 2, figsize=FIG_SIZE_IN)
        axes = axes.ravel()
        panel_a(axes[0], T['ident'], T['mine'], plotted, checks)
        panel_b(axes[1], T['pat'], plotted)
        check_cd(T['pairs'], T['src'], checks)
        for i in range(2):
            panel_cd(axes[i + 2], T['pairs'], i, plotted)
        fig.tight_layout(pad=1.6)
        for ext in ['pdf', 'svg', 'png']:
            fig.savefig(out / ('SFig07.' + ext), dpi=PNG_DPI if ext == 'png' else None)
        boxes = panel_boxes(fig, axes)
        plt.close(fig)
    pd.DataFrame(plotted).to_csv(out / 'plotted_values.tsv', sep='\t', index=False)
    return dict(plotted=pd.DataFrame(plotted), checks=checks, font=font, axes_bbox_px=boxes)


def panel_boxes(fig, axes):
    """Return each panel's pixel box (left, top, right, bottom) in the 300-dpi PNG.

    The box is the axes' tight bounding box (axes, labels, ticks, title) at ``PNG_DPI``.
    Used only to crop the new and the expected composite at identical pixels for
    per-panel receipts.
    """
    fig.set_dpi(PNG_DPI)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    H = fig.bbox.height
    boxes = {}
    for letter, ax in zip('abcd', axes):
        bb = ax.get_tightbbox(r)
        boxes[letter] = [max(0, int(np.floor(bb.x0))), max(0, int(np.floor(H - bb.y1))),
                         int(np.ceil(bb.x1)), int(np.ceil(H - bb.y0))]
    return boxes


def crop(src_png, box, dst_png):
    """Crop ``src_png`` to ``box`` (left, top, right, bottom) and write ``dst_png``."""
    from PIL import Image
    with Image.open(src_png) as im:
        w, h = im.size
        l, t, r, b = box
        im.crop((l, t, min(r, w), min(b, h))).save(dst_png)


# ---------------------------------------------------------------------------------------
# Statistics (all descriptive; no hypothesis test is made for this figure)
# ---------------------------------------------------------------------------------------
def stat_identity(mine, arm):
    """Median per-acquisition peptide-count gain and identity retention for one arm.

    Statistical unit: acquisition; n = 58 healthy-stool acquisitions. The gain is
    ``100 * (median(test / original) - 1)``; retention is ``100 * median(overlap / original)``,
    medians over acquisitions (numpy median, linear midpoint for even n).

    Parameters
    ----------
    mine : pandas.DataFrame
        ``identity_per_acquisition.tsv``.
    arm : str
        ``union_individual``, ``joint_fileq`` or ``joint_corrected``.

    Returns
    -------
    dict
        ``n``, ``gain_pct``, ``n_higher``, ``retention_pct``.
    """
    x = mine[mine.arm.eq(arm)]
    ratio = x.test / x.original
    rec = x.overlap / x.original
    return dict(n=int(len(x)), gain_pct=100 * (float(np.median(ratio)) - 1), n_higher=int((ratio > 1).sum()),
                retention_pct=100 * float(np.median(rec)))


def stat_study_distinct(summary, arm):
    """Distinct peptides across the 58-acquisition study for one arm (upstream count).

    Statistical unit: peptide (distinct canonical sequence, I/L equivalent); computed
    upstream from the accepted identity tables (``compare_fileq.py``) and read here.
    """
    return int(summary.loc[summary.arm.eq(arm), 'study_distinct'].iloc[0])


def stat_complete_features(completeness, arm, tier='TechnicalReplicates', space='features'):
    """LFQ features with an accepted same-acquisition identification in all replicates.

    Statistical unit: peptide feature; n acquisitions = 12 technical replicates. Computed
    upstream (``compare_fileq.py``) from the gated feature matrices and read here.
    """
    r = completeness[completeness.space.eq(space) & completeness.arm.eq(arm) & completeness.tier.eq(tier)]
    return int(r.rows_complete_all.iloc[0])


def stat_pair_precision(pairs, definition, comparison, method, tier):
    """Median over acquisition pairs of per-pair median microbial CV and log2 Pearson r.

    Statistical unit: acquisition pair (pairs share acquisitions and are dependent);
    n = 66 technical or 105 preparation pairs. Per-pair values (population SD / mean over
    the two raw intensities, median across protein-group rows; Pearson r of log2
    intensities on the same rows) were computed upstream (``compare_fileq.py``).

    Returns
    -------
    dict
        ``pairs``, ``median_cv_pct``, ``median_r``.
    """
    x = pairs[pairs.definition.eq(definition) & pairs.comparison.eq(comparison) & pairs.method.eq(method) & pairs.tier.eq(tier)]
    return dict(pairs=int(len(x)), median_cv_pct=float(np.median(x.median_pair_cv_pct)),
                median_r=float(np.median(x.log2_pearson)))


def stat_leading_species(agreement, tier, column, method='fixed_joint_corrected', threshold='both_resolved'):
    """Function-averaged leading-species agreement (percent), study-level confidence.

    Statistical unit: acquisition pair, averaged over matched resolved functions; computed
    upstream (A21 ``leaders_revision2.py``) and read here.
    """
    r = agreement[agreement.method.eq(method) & agreement.tier.eq(tier) & agreement.threshold.eq(threshold)]
    return 100 * float(r[column].iloc[0])


def stat_participants(pat):
    """Participants per cohort in panel b (unit participant)."""
    return pat.groupby('cohort').size().to_dict()


def stat_pairs_per_tier(pairs):
    """Acquisition pairs per tier in panels c, d (own_rows, original arm)."""
    x = pairs[pairs.definition.eq('own_rows') & pairs.method.eq('original')]
    return x.groupby('tier').size().to_dict()


def all_numbers(T):
    """Recompute every number in NUMBERS.tsv, in file order.

    Parameters
    ----------
    T : dict
        Output of ``load``.

    Returns
    -------
    list of (str, float)
        ``(key, full-precision value)``.
    """
    mine, pairs = T['mine'], T['pairs']
    u = stat_identity(mine, 'union_individual')
    j = stat_identity(mine, 'joint_fileq')
    s = stat_identity(mine, 'joint_corrected')
    pp = lambda d, c, m, t: stat_pair_precision(pairs, d, c, m, t)  # noqa: E731
    tf_o = pp('matched', 'original_vs_joint_fileq', 'original', 'TechnicalReplicates')
    tf_j = pp('matched', 'original_vs_joint_fileq', 'joint_fileq', 'TechnicalReplicates')
    wf_o = pp('matched', 'original_vs_joint_fileq', 'original', 'WorkflowReplicates')
    wf_j = pp('matched', 'original_vs_joint_fileq', 'joint_fileq', 'WorkflowReplicates')
    ts_o = pp('matched', 'original_vs_joint_corrected', 'original', 'TechnicalReplicates')
    ts_j = pp('matched', 'original_vs_joint_corrected', 'joint_corrected', 'TechnicalReplicates')
    ws_o = pp('matched', 'original_vs_joint_corrected', 'original', 'WorkflowReplicates')
    ws_j = pp('matched', 'original_vs_joint_corrected', 'joint_corrected', 'WorkflowReplicates')
    part = stat_participants(T['pat'])
    npairs = stat_pairs_per_tier(pairs)
    ag = T['agreement']
    return [
        ('L115 union_individual gain %', u['gain_pct']),
        ('L115 union_individual retention %', u['retention_pct']),
        ('L115 joint_fileq gain %', j['gain_pct']),
        ('L115 joint_fileq acquisitions higher', j['n_higher']),
        ('L115 joint_fileq retention %', j['retention_pct']),
        ('L115 complete features original', stat_complete_features(T['completeness'], 'original')),
        ('L115 complete features union_individual', stat_complete_features(T['completeness'], 'union_individual')),
        ('L115 complete features joint_fileq', stat_complete_features(T['completeness'], 'joint_fileq')),
        ('L115 distinct peptides original', stat_study_distinct(T['summary'], 'original')),
        ('L115 distinct peptides union_individual', stat_study_distinct(T['summary'], 'union_individual')),
        ('L115 distinct peptides joint_fileq', stat_study_distinct(T['summary'], 'joint_fileq')),
        ('L115 joint_corrected gain %', s['gain_pct']),
        ('L115 joint_corrected retention %', s['retention_pct']),
        ('L115 complete features joint_corrected', stat_complete_features(T['completeness'], 'joint_corrected')),
        ('L115 distinct peptides joint_corrected', stat_study_distinct(T['summary'], 'joint_corrected')),
        ('L437 PREDICT participants', part['PREDICT']),
        ('L437 CapScan participants', part['CAPSCAN']),
        ('L437 technical pairs', npairs['TechnicalReplicates']),
        ('L437 preparation pairs', npairs['WorkflowReplicates']),
        ('L769 acq-q technical CV original', tf_o['median_cv_pct']),
        ('L769 acq-q technical CV joint', tf_j['median_cv_pct']),
        ('L769 acq-q technical r original', tf_o['median_r']),
        ('L769 acq-q technical r joint', tf_j['median_r']),
        ('L769 acq-q preparation CV original', wf_o['median_cv_pct']),
        ('L769 acq-q preparation CV joint', wf_j['median_cv_pct']),
        ('L769 acq-q preparation r original', wf_o['median_r']),
        ('L769 acq-q preparation r joint', wf_j['median_r']),
        ('L769 study-q technical CV original', ts_o['median_cv_pct']),
        ('L769 study-q technical CV joint', ts_j['median_cv_pct']),
        ('L769 study-q technical r original', ts_o['median_r']),
        ('L769 study-q technical r joint', ts_j['median_r']),
        ('L769 study-q preparation CV original', ws_o['median_cv_pct']),
        ('L769 study-q preparation CV joint', ws_j['median_cv_pct']),
        ('L769 study-q preparation r original', ws_o['median_r']),
        ('L769 study-q preparation r joint', ws_j['median_r']),
        ('L769 leading species technical original %', stat_leading_species(ag, 'TechnicalReplicates', 'original_agreement')),
        ('L769 leading species technical joint %', stat_leading_species(ag, 'TechnicalReplicates', 'test_agreement')),
        ('L769 leading species preparation original %', stat_leading_species(ag, 'WorkflowReplicates', 'original_agreement')),
        ('L769 leading species preparation joint %', stat_leading_species(ag, 'WorkflowReplicates', 'test_agreement')),
    ]


def consistency(T):
    """Check that recomputed summaries equal the upstream summary tables.

    Returns
    -------
    list of dict
        One ``{check, passed}`` per comparison.
    """
    out = []
    for r in T['summary'].itertuples(index=False):
        x = T['mine'][T['mine'].arm.eq(r.arm)]
        ratio = x.test / x.original
        out.append(dict(check=f'identity_summary {r.arm} equals recomputation from identity_per_acquisition',
                        passed=bool(np.isclose(np.median(ratio), r.median_count_ratio, rtol=0, atol=1e-12)
                                    and int((ratio > 1).sum()) == r.n_higher
                                    and np.isclose(np.median(x.overlap / x.original), r.median_identity_recovery, rtol=0, atol=1e-12))))
    for r in T['pair_summary'].itertuples(index=False):
        p = stat_pair_precision(T['pairs'], r.definition, r.comparison, r.method, r.tier)
        out.append(dict(check=f'pair summary {r.definition} {r.comparison} {r.method} {r.tier}',
                        passed=bool(p['pairs'] == r.pairs and np.isclose(p['median_cv_pct'], r.median_pair_cv_pct, rtol=0, atol=1e-9)
                                    and np.isclose(p['median_r'], r.log2_pearson, rtol=0, atol=1e-12))))
    ident = T['ident'][T['ident'].min_length.eq(MIN_LENGTH_PANEL_A)]
    for arm in ['union_individual', 'joint_corrected']:
        a = ident[ident.test.eq(arm)].set_index('sample').original_identity_recovery.sort_index()
        b = T['mine'][T['mine'].arm.eq(arm)].set_index('sample').original_identity_recovery.sort_index()
        out.append(dict(check=f'identity_comparison (A21) and identity_per_acquisition agree for {arm}',
                        passed=bool(len(a) == len(b) == N_ACQUISITIONS and (a.index == b.index).all() and np.allclose(a.values, b.values, rtol=0, atol=1e-12))))
    return out
