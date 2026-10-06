"""Supplementary Figure 6 (v30): plotting and statistics functions.

Panels
------
a
    Schematic of the spectrum-split experiment: databases selected from all predictions
    (all-prediction database), odd-scan predictions (same-half database; ``self`` in the
    tables) or even-scan predictions (disjoint database), each searched against the
    identical odd-scan MS2 spectra.
b, c
    Disjoint / same-half (b) and disjoint / all-prediction (c) peptide-count ratio per
    acquisition at 9-50 residues, 726 acquisitions (15 CAMPI, 58 healthy stool, 306
    PREDICT, 347 CapScan).
d, e
    Same-spectrum benchmark against the matched metagenome search at 9-50 residues:
    peptides per acquisition (d) and the fraction of metagenome-search peptide identities
    also found by FastaLake (e); 296 PREDICT and 282 CapScan acquisitions.

Source data (``source_data/SFig06/``)
-------------------------------------
split_half_acquisitions_726.tsv (``SPLIT_TABLE``)
    panels b, c and the spectrum-split numbers.
same_spectrum_benchmark_0_predict.tsv, same_spectrum_benchmark_1_capscan.tsv (``BENCH_TABLES``)
    panels d, e and the same-spectrum numbers. These inputs are parameters: the d,e
    benchmark will be replaced by the CapScan metagenome benchmark extended to 295
    acquisitions (OPEN_ITEMS.md); pass ``bench_tables`` and ``bench_n`` to ``figure``.
split_half_structure.tsv
    acquisitions nested in specimens and participants (STATS.tsv).

How the figure is drawn
-----------------------
The drawing code is moved without change of logic from the v30 house-style build
(``build_sfig06_style.py``; PROVENANCE.md), which restyled the v29 producer
(``panels_len950.py``, ``_render_abc`` and ``_render_de``, minimum length 9, roster 726)
without changing any plotted value, and drew all five panels on one 180 x 186 mm canvas.
Style comes from ``style/fl_style``.
"""

import json
import os
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

from style import fl_style as fl

UNIT = Path(__file__).resolve().parent
REPO = UNIT.parents[1]
SD = REPO / "source_data" / UNIT.name

# ---------------------------------------------------------------- choices
MIN_LENGTH = 9                 # peptides of 9-50 residues; v30 legend (visible_v30f line 231)
MAX_LENGTH = 50                # applied upstream (record only, not re-applied here)
Q_CUTOFF = 0.01                # Sage peptide q <= 0.01; applied upstream (record only)
COHORT_ORDER = ["CAMPI_MOCK", "HSTOOL", "PREDICT", "CAPSCAN"]   # panel b,c x order (producer)
SPLIT_ROSTER_N = 726           # acquisitions in panels b,c (v30 legend line 231)
N_CAPSCAN = 347                # CapScan acquisitions in the b,c roster (tick label)
SPLIT_TABLE = "split_half_acquisitions_726.tsv"           # panels b,c input (parameter)
BENCH_TABLES = ("same_spectrum_benchmark_0_predict.tsv",  # panels d,e inputs (parameter; to be
                "same_spectrum_benchmark_1_capscan.tsv")  # replaced by the 295-acquisition MG benchmark)
JITTER_SEED_BC = 20260927      # producer: default_rng(20260927 + panel index), b = +0, c = +1
JITTER_SEED_DE = 4             # producer: default_rng(4) for d and e
BENCH_TEST_ARM = "FL_public"   # public-sequence FastaLake
BENCH_REF_ARM = "MG"           # matched metagenome search
BENCH_N = {"PREDICT": 296, "CAPSCAN": 282}   # matched acquisitions (v30 legend line 231)
FIG_W_MM, FIG_H_MM = 180.0, 186.0             # one canvas at the placed size (v30 build)
PNG_DPI = 300                  # PNG rasterised from the PDF at 300 dpi
PAIR_LINE_COLOR = "#cfcfcf"    # d pairing lines (v30 build LG)


# ---------------------------------------------------------------- data
def load_split(table=SPLIT_TABLE):
    """Read the spectrum-split acquisition table.

    Parameters
    ----------
    table : str
        File name in ``source_data/SFig06`` (default ``SPLIT_TABLE``, the 726 roster).

    Returns
    -------
    pandas.DataFrame
        One row per acquisition and length view (``min_length`` 5 or 9).
    """
    return pd.read_csv(SD / table, sep="\t")


def load_structure():
    """Read acquisitions, participants and specimens per cohort of the spectrum-split run.

    Returns
    -------
    pandas.DataFrame
        One row per cohort, indexed by cohort; counted by ``prepare_source_data.py`` on the
        710-acquisition run (CapScan 331) before the identifier columns were removed.
    """
    return pd.read_csv(SD / "split_half_structure.tsv", sep="\t").set_index("cohort")


def load_benchmark(tables=BENCH_TABLES, n=None):
    """Read the same-spectrum benchmark rows used by panels d and e.

    Parameters
    ----------
    tables : sequence of str
        File names in ``source_data/SFig06`` (default ``BENCH_TABLES``).
    n : dict or None
        Expected acquisitions per cohort (default ``BENCH_N``).

    Returns
    -------
    pandas.DataFrame
        Rows with ``min_length`` 9, ``test_arm`` FL_public and ``reference_arm`` MG, one
        per acquisition, exactly the producer's filter.
    """
    n = BENCH_N if n is None else n
    pa = pd.concat([pd.read_csv(SD / f, sep="\t") for f in tables])
    pa = pa[pa.min_length.eq(MIN_LENGTH) & pa.test_arm.eq(BENCH_TEST_ARM) & pa.reference_arm.eq(BENCH_REF_ARM)]
    assert pa.groupby("cohort").size().to_dict() == {"CAPSCAN": n["CAPSCAN"], "PREDICT": n["PREDICT"]}
    return pa


def font_dir():
    """Return the Söhne font directory from ``FASTALAKE_FONT_DIR`` (or ``FL_STYLE_DIR``), or None.

    Returns
    -------
    pathlib.Path or None
        The directory if it holds ``soehne-buch.ttf``, else None.
    """
    d = os.environ.get("FASTALAKE_FONT_DIR") or os.environ.get("FL_STYLE_DIR")
    return Path(d) if d and (Path(d) / "soehne-buch.ttf").is_file() else None


# ---------------------------------------------------------------- drawing (v30 build, unchanged logic)
# Colours by meaning: all-prediction database = neutral grey, same-half = FastaLake blue,
# disjoint = green; d,e: MG orange, FL blue.
COL = {'all': fl.GREY, 'same': fl.FL, 'disjoint': fl.GREEN}
LABELS_BC = ['CAMPI\n15 runs', 'Healthy stool\n58 runs', 'PREDICT\n306 runs', f'CapScan\n{N_CAPSCAN} runs']


def _axmm(fig, x0, y0, x1, y1, W=FIG_W_MM, H=FIG_H_MM):
    """Axes from a box in mm measured from the top-left corner."""
    return fig.add_axes([x0 / W, 1 - y1 / H, (x1 - x0) / W, (y1 - y0) / H])


def panel_a(fig):
    """Draw panel a, the spectrum-split schematic (no data).

    Parameters
    ----------
    fig : matplotlib.figure.Figure

    Returns
    -------
    matplotlib.axes.Axes
    """
    a = _axmm(fig, 20, 8, 176, 36); a.set(xlim=(0, 10), ylim=(0, 2.55)); a.axis('off')
    colors = [COL['all'], COL['same'], COL['disjoint']]
    for y, title, db, col in zip([2.1, 1.2, .3], ['All predictions', 'Odd-scan predictions', 'Even-scan predictions'], ['All-prediction database', 'Same-half database', 'Disjoint database'], colors):
        for x, width, text in [(0, 2.45, title), (3.15, 2.35, db)]:
            a.add_patch(FancyBboxPatch((x, y - .27), width, .59, boxstyle='round,pad=0.025,rounding_size=0.06', facecolor=col, edgecolor='none', alpha=.13))
            a.text(x + width / 2, y + .025, text, ha='center', va='center', fontsize=7, color=col)
        a.add_patch(FancyArrowPatch((2.5, y + .025), (3.1, y + .025), arrowstyle='-|>', mutation_scale=7, color=col, linewidth=1))
        a.add_patch(FancyArrowPatch((5.56, y + .025), (7.04, 1.225), arrowstyle='-|>', mutation_scale=7, color=col, linewidth=1))
    a.add_patch(FancyBboxPatch((7.1, .76), 2.55, .94, boxstyle='round,pad=0.035,rounding_size=.06', facecolor='#EEF1F3', edgecolor=fl.lighten(fl.GREY, .4), linewidth=1))
    a.text(8.375, 1.23, 'Search identical\nodd-scan spectra', ha='center', va='center', fontsize=7)
    return a


def panel_bc(fig, df):
    """Draw panels b (disjoint / same-half) and c (disjoint / all-prediction) peptide-count ratios.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    df : pandas.DataFrame
        ``load_split()`` rows at ``MIN_LENGTH`` (726 acquisitions).

    Returns
    -------
    list of matplotlib.axes.Axes
        Axes of b and c.
    """
    nzero = int((df[df.cohort == 'CAPSCAN'].disjoint == 0).sum())
    axes_bc = []
    for panel, col, ylab, letter, color in [(0, 'D_over_S', 'Disjoint / same-half\npeptide-count ratio (%)', 'b', COL['same']),
                                             (1, 'D_over_F', 'Disjoint / all-prediction\npeptide-count ratio (%)', 'c', COL['disjoint'])]:
        ax = _axmm(fig, 22, 48, 86, 101) if panel == 0 else _axmm(fig, 113, 48, 177, 101); fl.style(ax); axes_bc.append(ax)
        rng = np.random.default_rng(JITTER_SEED_BC + panel)
        for i, co in enumerate(COHORT_ORDER):
            v = df[df.cohort == co][col].to_numpy() * 100
            ax.scatter(i + rng.uniform(-.16, .16, len(v)), v, s=4, c=color, alpha=.3, edgecolors='none', rasterized=False, zorder=2)
            q1, med, q3 = np.quantile(v, [.25, .5, .75]); ax.plot([i, i], [q1, q3], color=fl.INK, lw=3, zorder=3); ax.plot([i - .12, i + .12], [med, med], color='white', lw=1.2, zorder=4)
            ax.text(i, 105, f'{med:.1f}%', ha='center', fontsize=6.5, color=color)
        ax.axhline(100, color=fl.GREY, lw=1, ls=(0, (3, 3)), zorder=1)
        ax.set(xlim=(-.5, 3.5), ylim=(-3, 112), xticks=range(4), xticklabels=LABELS_BC, yticks=[0, 25, 50, 75, 100], ylabel=ylab)
        ax.tick_params(axis='x', labelsize=6.5)
        if panel == 0:
            ax.annotate({1: 'One zero-yield run\nretained', 2: 'Two zero-yield runs\nretained'}[nzero], xy=(3, 0), xytext=(2, 19), fontsize=6, ha='center', arrowprops=dict(arrowstyle='-', lw=1, color=fl.GREY))
    return axes_bc


def panel_de(fig, pa):
    """Draw panels d (peptides per acquisition, MG vs FL) and e (MG identities recovered by FL).

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    pa : pandas.DataFrame
        ``load_benchmark()`` rows.

    Returns
    -------
    (Axes, Axes, list of dict)
        Axes of d and e, and the plotted values of d,e (one row per acquisition).
    """
    axD = _axmm(fig, 26, 122, 108, 171); axE = _axmm(fig, 131, 122, 177, 171); fl.style(axD); fl.style(axE)
    rng = np.random.default_rng(JITTER_SEED_DE); xs_ = {'PREDICT': (0, 1), 'CAPSCAN': (2.6, 3.6)}; pv6 = []
    for coh, (x0, x1) in xs_.items():
        d = pa[pa.cohort.eq(coh)]; j0 = rng.uniform(-.14, .14, len(d)); j1 = rng.uniform(-.14, .14, len(d))
        for a_, b_, ja, jb in zip(d.n_reference / 1e3, d.n_test / 1e3, j0, j1):
            axD.plot([x0 + ja, x1 + jb], [a_, b_], color=PAIR_LINE_COLOR, lw=.35, alpha=.35, zorder=1)
        axD.scatter(x0 + j0, d.n_reference / 1e3, s=3, color=fl.METAG, alpha=.7, lw=0, zorder=2)
        axD.scatter(x1 + j1, d.n_test / 1e3, s=3, color=fl.FL, alpha=.7, lw=0, zorder=2)
        for xx, v in ((x0, d.n_reference.median()), (x1, d.n_test.median())):
            axD.plot([xx - .3, xx + .3], [v / 1e3] * 2, color=fl.INK, lw=1.4, zorder=3)
            axD.text(xx + .34, v / 1e3, f'{v / 1e3:.1f}k', ha='left', va='center', fontsize=6.5, zorder=4)
        for _, r in d.iterrows():
            pv6.append(dict(panel='SuppFig6d,e', cohort=coh, sample=r['sample'], participant=r.subject_id, **{f'mg_peptides_{MIN_LENGTH}_50': int(r.n_reference)},
                            **{f'fl_peptides_{MIN_LENGTH}_50': int(r.n_test), f'recovery_pct_{MIN_LENGTH}_50': 100 * r.reference_recovery},
                            source='Figure4/benchmark_0_per_acquisition.tsv + benchmark_1_per_acquisition.tsv (min_length ' + str(MIN_LENGTH) + ', FL_public vs MG)'))
    axD.set_xticks([0, 1, 2.6, 3.6], ['MG', 'FL', 'MG', 'FL']); axD.set_ylabel(f'Peptides per acquisition\n(thousands; same spectra,\n{MIN_LENGTH}–50 residues)', labelpad=3)
    axD.set_ylim(0, 55); axD.set_xlim(-.6, 4.2)
    nP, nC = int(pa.cohort.eq('PREDICT').sum()), int(pa.cohort.eq('CAPSCAN').sum())
    axD.text(.5, -.11, f'PREDICT (n = {nP})', ha='center', va='top', transform=axD.get_xaxis_transform())
    axD.text(3.1, -.11, f'CapScan (n = {nC})', ha='center', va='top', transform=axD.get_xaxis_transform())
    for i, coh in enumerate(['PREDICT', 'CAPSCAN']):
        d = pa[pa.cohort.eq(coh)]; j = rng.uniform(-.18, .18, len(d))
        axE.scatter(i + j, 100 * d.reference_recovery, s=3, color=fl.FL, alpha=.7, lw=0)
        m = 100 * d.reference_recovery.median(); axE.plot([i - .3, i + .3], [m, m], color=fl.INK, lw=1.4, zorder=3)
        axE.text(i + .34, m, f'{m:.1f}%', ha='left', va='center', fontsize=6.5, zorder=4)
    axE.set_xticks([0, 1], ['PREDICT', 'CapScan']); axE.set_ylabel(f'MG peptides found by FL (%)\nsame spectra, {MIN_LENGTH}–50 residues', labelpad=3)
    axE.set_ylim(0, 105); axE.set_xlim(-.6, 1.9)
    return axD, axE, pv6


def align_letters(fig, items, w_mm=FIG_W_MM):
    """Draw panel letters with ``fl.panel`` and move each to ``x_mm`` from the left edge.

    Moved unchanged from ``_sf456_common.align_letters`` (v30 restyle helper).

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    items : list of (Axes, str, float)
    w_mm : float
    """
    for ax, letter, x_mm in items:
        fl.panel(ax, letter)
        pos = ax.get_position()
        t = ax._left_title
        t.set_x((x_mm / w_mm - pos.x0) / pos.width)
        t.set_ha('left')


def figure(out_dir, split_table=SPLIT_TABLE, bench_tables=BENCH_TABLES, bench_n=None):
    """Draw the v30 figure (a-e on one canvas), save PDF and 300 dpi PNG and the d,e plotted values.

    Parameters
    ----------
    out_dir : pathlib.Path
        Output directory; writes ``SFig06.pdf``, ``SFig06.png`` and ``panel_de_plotted_values.tsv``.
    split_table : str
        Panels b,c input (default ``SPLIT_TABLE``).
    bench_tables : sequence of str
        Panels d,e inputs (default ``BENCH_TABLES``).
    bench_n : dict or None
        Expected d,e acquisitions per cohort (default ``BENCH_N``).

    Returns
    -------
    dict
        Paths of the PDF, PNG and plotted-values table, and the PNG size.
    """
    import fitz
    out_dir = Path(out_dir); out_dir.mkdir(parents=True, exist_ok=True)
    fl.init()
    df = load_split(split_table); df = df[df.min_length == MIN_LENGTH]; assert len(df) == SPLIT_ROSTER_N
    pa = load_benchmark(bench_tables, bench_n)
    fig = plt.figure(figsize=fl.mm(FIG_W_MM, FIG_H_MM))
    a = panel_a(fig)
    axes_bc = panel_bc(fig, df)
    axD, axE, pv6 = panel_de(fig, pa)

    def x0mm(ax): return ax.get_position().x0 * FIG_W_MM
    align_letters(fig, [(a, 'a', x0mm(axes_bc[0]) - 19), (axes_bc[0], 'b', x0mm(axes_bc[0]) - 19), (axes_bc[1], 'c', x0mm(axes_bc[1]) - 19),
                        (axD, 'd', x0mm(axD) - 23), (axE, 'e', x0mm(axE) - 19)])
    pv = out_dir / 'panel_de_plotted_values.tsv'; pd.DataFrame(pv6).to_csv(pv, sep='\t', index=False)
    fl.save(fig, out_dir / 'SFig06', dpi=PNG_DPI, formats=('pdf',))
    pdf = out_dir / 'SFig06.pdf'
    pix = fitz.open(pdf)[0].get_pixmap(dpi=PNG_DPI, alpha=False); pix.set_dpi(PNG_DPI, PNG_DPI); pix.save(out_dir / 'SFig06.png')
    return {"pdf": pdf, "png": out_dir / 'SFig06.png', "plotted_values": pv, "png_px": [pix.width, pix.height]}


# ---------------------------------------------------------------- statistics
def stat_split_summary(df):
    """Cohort summaries of the spectrum-split experiment at ``MIN_LENGTH``-50 residues (9-50).

    Statistical unit: acquisition (n = 15 CAMPI, 58 healthy stool, 306 PREDICT and 347
    CapScan; 726 roster). Descriptive only: no
    test, no interval beyond the interquartile range. Medians and quartiles are
    ``pandas.Series.median`` and ``quantile`` (linear interpolation), as in the
    producer's ``review.py``; the plotted bars use ``numpy.quantile`` (same definition).

    Parameters
    ----------
    df : pandas.DataFrame
        Output of ``load_split``.

    Returns
    -------
    pandas.DataFrame
        One row per cohort with n, zero-yield counts, number of acquisitions with fewer
        disjoint than joint peptides, and median/q25/q75 of D/S, D/F and of the fraction
        of joint-search peptide identities also accepted with disjoint selection, all as
        fractions.
    """
    g5 = df[df.min_length == MIN_LENGTH]
    rows = []
    for co in COHORT_ORDER:
        g = g5[g5.cohort == co]
        r = dict(cohort=co, n=len(g), zero_disjoint=int((g.disjoint == 0).sum()),
                 D_lower_S=int((g.disjoint < g.self).sum()))
        for col in ("D_over_S", "D_over_F", "fraction_self_identities_retained"):
            v = g[col]
            r[f"median_{col}"] = float(v.median()); r[f"q25_{col}"] = float(v.quantile(.25)); r[f"q75_{col}"] = float(v.quantile(.75))
        rows.append(r)
    return pd.DataFrame(rows).set_index("cohort")


def stat_split_text(df):
    """Numbers printed for the spectrum-split experiment (visible_v30f lines 49, 64, 231).

    Statistical unit: acquisition; n = 726 (15 + 58 + 306 + 347). Descriptive medians.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of ``load_split``.

    Returns
    -------
    dict
        ``n_total``; ``median_DS_pct`` per cohort (100 x median D/S); ``loss_pct`` per
        cohort (100 x (1 - median D/S)); ``retention_pct`` per cohort (100 x median
        fraction of joint identities retained) and its minimum and maximum over cohorts;
        ``n_capscan_zero_disjoint``; ``all_fewer`` (True if every acquisition has fewer
        disjoint than joint peptides).
    """
    s = stat_split_summary(df)
    ret = 100 * s["median_fraction_self_identities_retained"]
    return dict(n_total=int(s.n.sum()),
                n_by_cohort=s.n.to_dict(),
                median_DS_pct=(100 * s.median_D_over_S).to_dict(),
                median_DF_pct=(100 * s.median_D_over_F).to_dict(),
                loss_pct=(100 * (1 - s.median_D_over_S)).to_dict(),
                retention_pct=ret.to_dict(), retention_min_pct=float(ret.min()), retention_max_pct=float(ret.max()),
                n_capscan_zero_disjoint=int(s.loc["CAPSCAN", "zero_disjoint"]),
                all_fewer=bool((s.D_lower_S == s.n).all()))


def stat_bench(pa):
    """Same-spectrum benchmark summaries for panels d,e (visible_v30f lines 64, 65).

    Statistical unit: acquisition (PREDICT n = 296, CapScan n = 282). Descriptive
    medians; no test and no interval. Gain is ``100 * (median(count_ratio) - 1)`` with
    ``count_ratio`` = FastaLake / metagenome peptides in the same acquisition, as in the
    producer. Recovery is the median over acquisitions of |FL and MG| / |MG|.

    Parameters
    ----------
    pa : pandas.DataFrame
        Output of ``load_benchmark`` (one row per acquisition, 9-50 residues,
        FL_public vs MG).

    Returns
    -------
    dict
        Per cohort: n, median MG and FL peptides, median recovery (%), number of
        acquisitions with more FL than MG peptides, median gain (%).
    """
    out = {}
    for coh in ("PREDICT", "CAPSCAN"):
        d = pa[pa.cohort.eq(coh)]
        out[coh] = dict(n=len(d), median_MG=float(d.n_reference.median()), median_FL=float(d.n_test.median()),
                        median_recovery_pct=float(100 * d.reference_recovery.median()),
                        n_FL_gt_MG=int((d.n_test > d.n_reference).sum()),
                        median_gain_pct=float(100 * (d.count_ratio.median() - 1)))
    return out


NUMBER_WORDS = {0: "zero", 1: "one", 2: "two", 3: "three"}


def as_word(n):
    """Return the English word the manuscript uses for a small count (``1`` -> ``'one'``).

    Parameters
    ----------
    n : int
        Count from 0 to 3.

    Returns
    -------
    str
        The word.
    """
    return NUMBER_WORDS[int(n)]
