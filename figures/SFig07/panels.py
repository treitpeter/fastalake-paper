"""Supplementary Figure 7: common search databases (depth, identity retention, precision).

v30 figure. Two steps, both moved from the producers as they ran (see PROVENANCE.md):

1. The data producer (``build_supp_fig7.py``, as used by the 9-50-residue rebuild
   ``panels_len950.py``): functions ``panel_a``, ``panel_b``, ``panel_cd`` and ``render``.
   Panel a and the peptide axis of panel b use peptides of 9-50 residues
   (``MIN_LENGTH_PANEL_A = 9`` and the 9-50 tables below).
2. The v30 house-style restyle (``build_sfig07_style.py``): functions ``v30_*`` and
   ``render_v30``. It runs the producer unchanged into a hidden figure, redraws the same
   values with ``style/fl_style`` (colours by meaning, Söhne, panel letters, no titles),
   and compares every plotted array of the two figures exactly (``values_check``).

Source data (``source_data/SFig07/``):

- ``identity_comparison.tsv``: per-acquisition identity retention, study-level arms (panel a)
- ``identity_per_acquisition.tsv``: per-acquisition identity retention at 9-50 residues, all
  arms including per-acquisition confidence (panel a, Results)
- ``identity_summary_by_length.tsv``: per-arm summary at 5-50 and 9-50 residues (Results, Methods)
- ``completeness_by_length.tsv``: features complete in all replicates, 5-50 and 9-50 (Results)
- ``common_search_participants.tsv``: PREDICT and CapScan participant summaries at 9-50
  residues, pseudonymised identifiers (panel b, Methods)
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
MIN_LENGTH_PANEL_A = 9          # v30: panels a,b use 9-50-residue peptides (legend "Peptide counts in a,b are 9-50 residues"); v28d value 5
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
                src=r('common_search_pair_CV.tsv'), summary_by_length=r('identity_summary_by_length.tsv'),
                completeness_by_length=r('completeness_by_length.tsv'), pair_summary=r('microbial_group_pair_summary.tsv'),
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


def stat_study_distinct(summary_by_length, arm, min_length=MIN_LENGTH_PANEL_A):
    """Distinct peptides across the 58-acquisition study for one arm (upstream count).

    Statistical unit: peptide (distinct canonical sequence, I/L equivalent); n = 58
    acquisitions pooled. Computed upstream (``hstool_recount.py``) and read here.
    """
    r = summary_by_length[summary_by_length.min_length.eq(min_length) & summary_by_length.arm.eq(arm)]
    return int(r.study_distinct.iloc[0])


def stat_distinct_gain(summary_by_length, arm, min_length=MIN_LENGTH_PANEL_A):
    """Percent more distinct study peptides in ``arm`` than in the acquisition-specific searches.

    Statistical unit: peptide; n = 58 acquisitions pooled. ``100 * (arm / original - 1)``.
    """
    return 100 * (stat_study_distinct(summary_by_length, arm, min_length)
                  / stat_study_distinct(summary_by_length, 'original', min_length) - 1)


def stat_complete_features(completeness_by_length, arm, tier='TechnicalReplicates', space='features',
                           min_length=MIN_LENGTH_PANEL_A):
    """LFQ features with an accepted same-acquisition identification in all replicates.

    Statistical unit: peptide feature; n acquisitions = 12 technical replicates. Computed
    upstream (``hstool_recount.py``) from the gated feature matrices and read here.
    """
    c = completeness_by_length
    r = c[c.min_length.eq(min_length) & c.space.eq(space) & c.arm.eq(arm) & c.tier.eq(tier)]
    return int(r.rows_complete_all.iloc[0])


def stat_participant_median(pat, cohort):
    """Median participant-level percent change of accepted 9-50-residue peptides, common vs own FASTA.

    Statistical unit: participant; n = 93 PREDICT patients or 15 CapScan participants.
    Per-participant values are ``100 * (common / original - 1)`` after averaging
    acquisitions within specimen and specimens within participant (upstream).
    """
    return float(pat.loc[pat.cohort.eq(cohort), 'peptides_percent_difference'].median())


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
    cb, sb, pat = T['completeness_by_length'], T['summary_by_length'], T['pat']
    return [
        ('R56 union_individual gain %', u['gain_pct']),
        ('R56 joint_fileq gain %', j['gain_pct']),
        ('R56 union_individual retention %', u['retention_pct']),
        ('R56 joint_fileq retention %', j['retention_pct']),
        ('R56 complete features original', stat_complete_features(cb, 'original')),
        ('R56 complete features joint_fileq', stat_complete_features(cb, 'joint_fileq')),
        ('R56 joint_corrected gain %', s['gain_pct']),
        ('M141 distinct peptides union_individual gain %', stat_distinct_gain(sb, 'union_individual')),
        ('M141 PREDICT median participant peptide change %', stat_participant_median(pat, 'PREDICT')),
        ('M141 PREDICT participants', part['PREDICT']),
        ('M141 CapScan median participant peptide change %', stat_participant_median(pat, 'CAPSCAN')),
        ('M141 CapScan participants', part['CAPSCAN']),
        ('SF7 legend PREDICT participants', part['PREDICT']),
        ('SF7 legend CapScan participants', part['CAPSCAN']),
        ('SF7 legend technical pairs', npairs['TechnicalReplicates']),
        ('SF7 legend preparation pairs', npairs['WorkflowReplicates']),
        ('M394 acq-q technical CV original', tf_o['median_cv_pct']),
        ('M394 acq-q technical CV joint', tf_j['median_cv_pct']),
        ('M394 acq-q technical r original', tf_o['median_r']),
        ('M394 acq-q technical r joint', tf_j['median_r']),
        ('M394 acq-q preparation CV original', wf_o['median_cv_pct']),
        ('M394 acq-q preparation CV joint', wf_j['median_cv_pct']),
        ('M394 acq-q preparation r original', wf_o['median_r']),
        ('M394 acq-q preparation r joint', wf_j['median_r']),
        ('M394 study-q technical CV original', ts_o['median_cv_pct']),
        ('M394 study-q technical CV joint', ts_j['median_cv_pct']),
        ('M394 study-q technical r original', ts_o['median_r']),
        ('M394 study-q technical r joint', ts_j['median_r']),
        ('M394 study-q preparation CV original', ws_o['median_cv_pct']),
        ('M394 study-q preparation CV joint', ws_j['median_cv_pct']),
        ('M394 study-q preparation r original', ws_o['median_r']),
        ('M394 study-q preparation r joint', ws_j['median_r']),
        ('M394 leading species technical original %', stat_leading_species(ag, 'TechnicalReplicates', 'original_agreement')),
        ('M394 leading species technical joint %', stat_leading_species(ag, 'TechnicalReplicates', 'test_agreement')),
        ('M394 leading species preparation original %', stat_leading_species(ag, 'WorkflowReplicates', 'original_agreement')),
        ('M394 leading species preparation joint %', stat_leading_species(ag, 'WorkflowReplicates', 'test_agreement')),
    ]


def consistency(T):
    """Check that recomputed summaries equal the upstream summary tables.

    Returns
    -------
    list of dict
        One ``{check, passed}`` per comparison.
    """
    out = []
    for r in T['summary_by_length'][T['summary_by_length'].min_length.eq(MIN_LENGTH_PANEL_A)].itertuples(index=False):
        x = T['mine'][T['mine'].arm.eq(r.arm)]
        ratio = x.test / x.original
        out.append(dict(check=f'identity_summary_by_length (9-50) {r.arm} equals recomputation from identity_per_acquisition',
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


# =======================================================================================
# v30 house-style restyle (moved from build_sfig07_style.py; logic unchanged)
# =======================================================================================
import json  # noqa: E402

from matplotlib.collections import LineCollection, PathCollection  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from style import fl_style as fl  # noqa: E402

V30_W_MM, V30_H_MM = 160.0, 144.0   # placed size in the docx, 160.0 x 144.2 mm
V30_TICK_PT, V30_CAT_PT, V30_LAB_PT, V30_LEG_PT = 6.5, 5.6, 7.0, 6.5   # type sizes (pt) at placed size
V30_DOT_S, V30_DOT_A = 9, 0.6       # dot size (pt^2) and alpha

C_ORIG = fl.REF        # original / repeat individual search (the reference)
C_STUDY = fl.FL        # joint search, study-level confidence (study q)
C_ACQ = fl.SKY         # joint search, per-acquisition confidence (acquisition q)
C_JOINT_ORIG = fl.PURPLE   # joint search, original (uncorrected) build
C_UNION = fl.GREEN     # common-FASTA (union) individual search
#: v35: drawn names as in the manuscript; the producer/restyle check maps the producer's strings through these.
V35_COHORT = {"CAPSCAN": "CapScan"}
V35_TICK = {"union individual": "common-FASTA individual"}

V30_XLABEL_A = "Healthy stool, 58 acquisitions"
#: v38: d "Preparation replicates" -> "Workflow replicates", the name of this series in Fig. 3a and Supplementary Fig. 12.
V30_XLABEL_CD = ["Technical replicates, 66 pairs", "Workflow replicates, 105 pairs"]


def v30_collect(axes):
    """Every plotted quantity of each axes, colour-free, for the before/after check.

    Parameters
    ----------
    axes : sequence of matplotlib.axes.Axes

    Returns
    -------
    list of dict
        Limits, scales, ticks, scatter offsets, segments, lines, bars, texts, labels, legend.
    """
    axes[0].figure.canvas.draw()
    out = []
    for ax in axes:
        d = dict(
            xlim=list(ax.get_xlim()), ylim=list(ax.get_ylim()),
            xscale=ax.get_xscale(), yscale=ax.get_yscale(),
            xticks=[float(v) for v in ax.get_xticks()], yticks=[float(v) for v in ax.get_yticks()],
            xticklabels=[t.get_text() for t in ax.get_xticklabels()],
            scatter=[np.asarray(c.get_offsets(), float) for c in ax.collections if isinstance(c, PathCollection)],
            segments=[np.asarray(s, float) for c in ax.collections if isinstance(c, LineCollection) for s in c.get_segments()],
            lines=[np.asarray(l.get_xydata(), float) for l in ax.lines],
            bars=[(p.get_x(), p.get_y(), p.get_width(), p.get_height()) for p in ax.patches if isinstance(p, Rectangle)],
            texts=[(tuple(map(float, t.get_position())), t.get_text()) for t in ax.texts],
            ylabel=ax.get_ylabel(), xlabel=ax.get_xlabel(), title=ax.get_title(loc="left"),
            legend=[t.get_text() for t in ax.get_legend().get_texts()] if ax.get_legend() else [])
        out.append(d)
    return out


def v30_compare(before, after, label_changes):
    """Exact comparison of plotted arrays of the producer and the restyled figure.

    Parameters
    ----------
    before, after : list of dict
        ``v30_collect`` of the producer and the restyled axes.
    label_changes : dict
        Panel letter -> expected new x label (documented title moves).

    Returns
    -------
    list of dict
        One ``{panel, check, passed, detail}`` record per check.
    """
    rec = []
    for i, (b, a) in enumerate(zip(before, after)):
        L = "abcd"[i]

        def ok(name, cond, detail=""):
            rec.append(dict(panel=L, check=name, passed=bool(cond), detail=detail))
        ok("scatter arrays identical (count, order, x, y)",
           len(b["scatter"]) == len(a["scatter"]) and all(x.shape == y.shape and np.array_equal(x, y) for x, y in zip(b["scatter"], a["scatter"])),
           f"{len(a['scatter'])} collections, {sum(len(x) for x in a['scatter'])} points")
        ok("median bars identical", len(b["segments"]) == len(a["segments"]) and all(np.array_equal(x, y) for x, y in zip(b["segments"], a["segments"])),
           f"{len(a['segments'])} bars")
        ok("reference lines identical", len(b["lines"]) == len(a["lines"]) and all(np.array_equal(x, y) for x, y in zip(b["lines"], a["lines"])),
           f"{len(a['lines'])} lines")
        ok("bars identical", b["bars"] == a["bars"])
        ok("data-label texts identical", b["texts"] == a["texts"])
        ok("x and y limits identical", np.allclose(b["xlim"], a["xlim"], rtol=0, atol=1e-12) and np.allclose(b["ylim"], a["ylim"], rtol=0, atol=1e-12),
           f"x {a['xlim']} y {a['ylim']}")
        ok("scales identical", (b["xscale"], b["yscale"]) == (a["xscale"], a["yscale"]))
        ok("x tick positions and labels identical (line breaks may differ)", b["xticks"] == a["xticks"] and [V35_TICK.get(" ".join(t.split()), " ".join(t.split())) for t in b["xticklabels"]] == [" ".join(t.split()) for t in a["xticklabels"]])
        ok("y tick positions identical", np.allclose(b["yticks"], a["yticks"], rtol=0, atol=1e-12))
        ok("y label identical", b["ylabel"] == a["ylabel"], a["ylabel"])
        exp_x = label_changes.get(L, b["xlabel"])
        ok("x label identical or the documented title move", a["xlabel"] == exp_x, f"{b['xlabel']!r} -> {a['xlabel']!r}")
        ok("legend entries identical", [V35_COHORT.get(t, t) for t in b["legend"]] == a["legend"], str(a["legend"]))
    return rec


def v30_reference(sd):
    """Run the producer unchanged into a hidden figure and collect its plotted values.

    Returns
    -------
    tuple
        ``(tables, collected values, plotted DataFrame, producer checks)``.
    """
    T = load(sd)
    rc, _ = producer_rc()
    with plt.rc_context():
        matplotlib.rcdefaults()
        plt.rcParams.update(rc)
        fig, axes = plt.subplots(2, 2, figsize=FIG_SIZE_IN)
        axes = axes.ravel()
        plotted, checks = [], []
        panel_a(axes[0], T["ident"], T["mine"], plotted, checks)
        panel_b(axes[1], T["pat"], plotted)
        for i in range(2):
            panel_cd(axes[i + 2], T["pairs"], i, plotted)
        fig.tight_layout(pad=1.6)
        vals = v30_collect(axes)
        plt.close(fig)
    return T, vals, pd.DataFrame(plotted), checks


def v30_dots(ax, groups, ylabel, percent=False):
    """Producer ``dots()`` in house style: same x (linspace jitter), values and median bars."""
    for i, (label, values, c) in enumerate(groups):
        v = np.asarray(values, dtype=float)
        v = v[np.isfinite(v)] * (100 if percent else 1)
        x = i + np.linspace(-.16, .16, len(v))
        ax.scatter(x, v, s=V30_DOT_S, alpha=V30_DOT_A, color=c, linewidths=0, rasterized=len(v) > 100)
        if len(v):
            ax.hlines(np.median(v), i - .25, i + .25, color=fl.INK, lw=1.5)
    ax.set_xticks(range(len(groups)), [g[0] for g in groups])
    ax.set_ylabel(ylabel)


def v30_panel_a(ax, ident, mine):
    """Panel a (v30): original identities recovered, 58 acquisitions per column, 9-50 residues."""
    ident = ident[ident.min_length.eq(MIN_LENGTH_PANEL_A)]
    fq = mine[mine.arm.eq("joint_fileq")]
    A = [("joint\ncorrected\nstudy q", ident[ident.test.eq("joint_corrected")].original_identity_recovery, C_STUDY),
         ("joint\noriginal\nstudy q", ident[ident.test.eq("joint_original")].original_identity_recovery, C_JOINT_ORIG),
         ("joint\ncorrected\nacquisition q", fq.original_identity_recovery, C_ACQ),
         ("repeat\nindividual", ident[ident.test.eq("repeat_individual")].original_identity_recovery, C_ORIG),
         ("common-FASTA\nindividual", ident[ident.test.eq("union_individual")].original_identity_recovery, C_UNION)]   # v35: legend name
    v30_dots(ax, A, "Original identities recovered (%)", True)
    ax.set_xlabel(V30_XLABEL_A)


def v30_panel_b(ax, pat):
    """Panel b (v30): participant changes, common vs own FASTA (93 PREDICT, 15 CapScan participants)."""
    for co, g in pat.groupby("cohort"):
        if co == "CAPSCAN":
            ax.scatter(g.peptides_percent_difference, g.directlfq_groups_percent_difference, s=12,
                       facecolors="none", edgecolors=C_UNION, linewidths=1.0, label=V35_COHORT.get(co, co), zorder=3)
        else:
            ax.scatter(g.peptides_percent_difference, g.directlfq_groups_percent_difference, s=12,
                       alpha=V30_DOT_A, color=C_UNION, linewidths=0, label=V35_COHORT.get(co, co), zorder=2)
    ax.axhline(0, c=fl.GREY, lw=1.0, zorder=1)
    ax.axvline(0, c=fl.GREY, lw=1.0, zorder=1)
    ax.set(xlabel="Accepted peptides: common vs own (%)", ylabel="DirectLFQ groups: common vs own (%)")
    ax.legend(frameon=False, fontsize=V30_LEG_PT, loc="upper left", handletextpad=0.3, borderaxespad=0.2)


def v30_panel_cd(ax, pairs, i):
    """Panels c (i = 0, 66 technical pairs) and d (i = 1, 105 preparation pairs), v30 style."""
    t = TIERS[i]
    s = _sel
    G = [("Own\noriginal", s(pairs, "own_rows", "original", "original", t), C_ORIG),
         ("Own\njoint\nstudy q", s(pairs, "own_rows", "joint_corrected", "joint_corrected", t), C_STUDY),
         ("Own\njoint\nacq. q", s(pairs, "own_rows", "joint_fileq", "joint_fileq", t), C_ACQ),
         ("Matched\noriginal", s(pairs, "matched", "original_vs_joint_corrected", "original", t), C_ORIG),
         ("Matched\njoint\nstudy q", s(pairs, "matched", "original_vs_joint_corrected", "joint_corrected", t), C_STUDY),
         ("Matched\noriginal", s(pairs, "matched", "original_vs_joint_fileq", "original", t), C_ORIG),
         ("Matched\njoint\nacq. q", s(pairs, "matched", "original_vs_joint_fileq", "joint_fileq", t), C_ACQ)]
    v30_dots(ax, G, "Median pair CV (%)")
    ax.axvline(2.5, c=fl.GREY, lw=1.0)
    ax.axvline(4.5, c=fl.GREY, lw=1.0, ls=":")
    ax.set_xlabel(V30_XLABEL_CD[i])


def render_v30(sd, out):
    """Draw the v30 figure and check its plotted values against the producer.

    Call ``fl.init()`` first. Writes ``SupplementaryFigure7.pdf/.png`` (300 dpi),
    ``values_check.json`` and ``plotted_values.tsv`` to ``out``.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig07``.
    out : pathlib.Path

    Returns
    -------
    dict
        ``checks`` (all records), ``png``, ``pdf``, ``plotted``.
    """
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    T, before, plotted_ref, prod_checks = v30_reference(sd)
    fl.init()
    fig, axes = plt.subplots(2, 2, figsize=fl.mm(V30_W_MM, V30_H_MM))
    axes = axes.ravel()
    fig.subplots_adjust(left=0.075, right=0.985, bottom=0.095, top=0.965, wspace=0.22, hspace=0.36)
    v30_panel_a(axes[0], T["ident"], T["mine"])
    v30_panel_b(axes[1], T["pat"])
    for i in range(2):
        v30_panel_cd(axes[i + 2], T["pairs"], i)
    for ax, L in zip(axes, "abcd"):
        fl.style(ax)
        ax.tick_params(labelsize=V30_TICK_PT)
        ax.xaxis.label.set_size(V30_LAB_PT)
        ax.yaxis.label.set_size(V30_LAB_PT)
        fl.panel(ax, L)
    for ax in (axes[0], axes[2], axes[3]):
        ax.tick_params(axis="x", labelsize=V30_CAT_PT)
    # keep the producer tick values (the auto locator depends on axes size in pixels)
    for ax, b in zip(axes, before):
        ax.set_yticks(b["yticks"])
        if ax is axes[1]:
            ax.set_xticks(b["xticks"])
        ax.set_xlim(b["xlim"])
        ax.set_ylim(b["ylim"])
    fig.canvas.draw()
    after = v30_collect(axes)
    rec = v30_compare(before, after, {"a": V30_XLABEL_A, "c": V30_XLABEL_CD[0], "d": V30_XLABEL_CD[1]})
    for L, b, a in zip("abcd", before, after):
        rec.append(dict(panel=L, check="producer title removed; content placed (x label or legend)", passed=a["title"] in ("", L),
                        detail=f"{b['title']!r} -> x label {a['xlabel']!r}" if L != "b" else f"{b['title']!r} -> legend b 'Clinical common-FASTA comparison'"))
    rec += [dict(panel="producer", check=c["check"], passed=c["passed"], detail="") for c in prod_checks]
    (out / "values_check.json").write_text(json.dumps(rec, indent=1) + "\n")
    plotted_ref.to_csv(out / "plotted_values.tsv", sep="\t", index=False)
    bad = [r for r in rec if not r["passed"]]
    if bad:
        raise SystemExit("VALUE CHECK FAILED: " + json.dumps(bad, indent=1))
    pdf, png = fl.save(fig, out / "SupplementaryFigure7", formats=("pdf", "png"))
    boxes = panel_boxes(fig, axes)
    plt.close(fig)
    return dict(checks=rec, pdf=Path(pdf), png=Path(png), plotted=plotted_ref, axes_bbox_px=boxes)
