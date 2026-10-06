"""Supplementary Figure 4: peptide-length sensitivity of FastaLake results (panels a-g), v30.

The drawing code is moved from the v30 house-style build (``build_sfig04_style.py``, which
restyled the v29 producer ``sf4.py`` without changing any plotted value; both listed in
PROVENANCE.md) with paths made repository-relative and the style taken from
``style/fl_style``. The figure is drawn as one composite: the panels share one figure, one
grid and one random-number generator (the jitter of panels a and b), so ``build_figure``
draws all seven panels in the producer's order. Per-panel images are crops of that
composite at the boxes computed by ``panel_boxes``.

Source data (``source_data/SFig04/``):

- ``benchmark_matched_per_acquisition.tsv``: a, b and the Supplementary Methods numbers
- ``per_acquisition_length_counts.tsv`` and ``per_acquisition_length_counts_participant15.tsv``: c
- ``campi_member_counts.tsv``: d
- ``campi_fusion_pooled_q01.tsv`` (drawn) and ``campi_fusion_per_acquisition_q01.tsv`` (pooling check): e
- ``healthy_strategies_per_acquisition.tsv``: f
- ``grouping_study_groups.tsv``: g
- ``search_inventory.tsv``: the legend search count
"""

import collections
import hashlib
import json
import os
import warnings
from functools import lru_cache
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import fitz  # PyMuPDF
import matplotlib.pyplot as plt
import matplotlib.text as mtext
import numpy as np
import pandas as pd
from matplotlib import font_manager as fm

from style import fl_style as fl

# ---------------------------------------------------------------- choices (named constants)
PLACED_W_MM = 160.0          # placed width in the manuscript (producer constant)
PLACED_H_MM = 208.0          # placed height (producer figsize)
INK = "#4e4f51"              # text and axis colour (producer common.py)
JITTER_SEED = 1              # np.random.default_rng(1), horizontal jitter of a and b only
JITTER_HALF_WIDTH = 0.25     # uniform jitter in x +- 0.25
GAIN_CLIP = 250              # legend v28d line 421: "Values above 250% are drawn at 250% and counted"
LOG_OFFSET = 0.5             # legend v28d line 421: "with 0.5 added for the logarithmic axis"
FUSION_REF_LINE = 1          # legend v28d line 421: "Dashed line, 1%"
FUSION_Q = 0.01              # legend v28d line 421: "at peptide q<=0.01" (file POOLED_q01.tsv)
SLANT_SHEAR = 0.20           # synthetic slant of the three taxon labels in d (no Soehne italic)
PNG_DPI = 300                # v30 build: 300 dpi PNG from the PDF (1890 x 2457 px)
MIN_TEXT_PT = 5              # smallest text at placed size, asserted by the producer
CAPSCAN_ROSTER = 347         # v28f CapScan roster for c (FL_public) and g (FL_cohortMG)
MATCHED_N = {"PREDICT": 296, "CAPSCAN": 282}   # legend v28d line 421, matched FL_public vs MG pairs
POPS_AB = ["baseline5", "fresh7", "nested9"]   # populations drawn in a and b
D_MEMBERS = ["Lactobacillus plantarum", "Clostridium butyricum", "Bifidobacterium longum"]  # three least represented
TAXON_LABELS = {"L. plantarum", "C. butyricum", "B. longum"}
# Colours by meaning (v30 restyle): the five length settings are not database arms, so FL blue
# and MG orange are not used; the production baseline stays neutral grey.
C = {"baseline5": fl.GREY, "nested7": fl.lighten(fl.GREY, 0.55), "fresh7": fl.SKY, "nested9": fl.VERM, "fresh7_nested9": fl.PURPLE}
# '7-50' / '9-50' in place of '>=7' / '>=9' (Soehne has no >= glyph; same meaning, maximum 50 residues)
LAB = {"baseline5": "5–50 (published)", "nested7": "7–50, filtered", "fresh7": "7–50, new search",
       "nested9": "9–50, filtered", "fresh7_nested9": "9–50 subset of new search"}
LETTER_X_MM = (1.0, 82.0)    # panel-letter x per column (v30 build XL, XR)
C_ARMS = [("PREDICT", "FL_public", "PREDICT FL"), ("PREDICT", "MG", "PREDICT MG"), ("CAPSCAN", "FL_public", "CapScan FL"),
          ("CAPSCAN", "MG", "CapScan MG"), ("HSTOOL", "FL", "Healthy FL"), ("CAMPI", "FL", "CAMPI FL")]
MM = 72 / 25.4

UNIT = Path(__file__).resolve().parent
SD_DEFAULT = UNIT.parents[1] / "source_data" / UNIT.name

# ---------------------------------------------------------------- fonts and style (common.py)
_FONT_DIR = os.environ.get("FASTALAKE_FONT_DIR") or os.environ.get("FL_STYLE_DIR")
BUCH_TTF = Path(_FONT_DIR) / "soehne-buch.ttf" if _FONT_DIR else None
FONTS_PRESENT = bool(BUCH_TTF and BUCH_TTF.is_file())
if FONTS_PRESENT:
    fm.fontManager.addfont(str(BUCH_TTF))
    BUCH = fm.FontProperties(fname=str(BUCH_TTF)).get_name()
else:
    warnings.warn("soehne-buch.ttf not found in FASTALAKE_FONT_DIR; text is drawn in DejaVu Sans and pixel checks are skipped")
    BUCH = "DejaVu Sans"


def sha(p):
    """Return the hex sha256 of a file."""
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def init():
    """Apply the house style (``style/fl_style.init``), as the v30 build did."""
    fl.init()


def spans(pdf):
    """Text spans of the first PDF page (common.py)."""
    pg = fitz.open(pdf)[0]
    return [dict(s, direction=l['dir']) for b in pg.get_text('dict')['blocks'] for l in b.get('lines', []) for s in l['spans']]


def slant(pdf_in, pdf_out, labels, shear=SLANT_SHEAR):
    """Re-set the given labels in Soehne Buch with a synthetic slant (common.py, unchanged).

    Removes only the label text and inserts it again with a shear matrix.
    """
    doc = fitz.open(pdf_in); page = doc[0]
    hits = [s for s in spans(pdf_in) if s['text'] in labels]
    for s in hits:
        page.add_redact_annot(fitz.Rect(s['bbox']), fill=None, cross_out=False)
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE, text=fitz.PDF_REDACT_TEXT_REMOVE)
    if FONTS_PRESENT:
        page.insert_font(fontname='SB', fontfile=str(BUCH_TTF)); fname = 'SB'
    else:
        fname = 'helv'
    for s in hits:
        col = tuple(((s['color'] >> k) & 255) / 255 for k in (16, 8, 0))
        pt = fitz.Point(s['origin'])
        page.insert_text(pt, s['text'], fontname=fname, fontsize=s['size'], color=col, morph=(pt, fitz.Matrix(1, 0, shear, 1, 0, 0)))
    doc.save(pdf_out, garbage=4, deflate=True)
    before = collections.Counter(s['text'] for s in spans(pdf_in) if s['text'].strip())
    after = collections.Counter(s['text'] for s in spans(pdf_out) if s['text'].strip())
    assert before == after, (before - after, after - before)
    return len(hits)


def render_png(pdf, png, dpi=PNG_DPI):
    """Rasterise the first PDF page (common.py)."""
    pix = fitz.open(pdf)[0].get_pixmap(dpi=dpi, alpha=False)
    pix.set_dpi(dpi, dpi); pix.save(png)
    return pix.width, pix.height


def min_text_pt(pdf, placed_width_mm):
    """Smallest text size at the placed width (common.py)."""
    pg = fitz.open(pdf)[0]; k = placed_width_mm * MM / pg.rect.width
    sz = [s['size'] * k for s in spans(pdf) if s['text'].strip()]
    return round(min(sz), 3), k


# ---------------------------------------------------------------- data
@lru_cache(maxsize=None)
def load(sd=SD_DEFAULT):
    """Read every source-data table of the unit, as the producer did.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig04``.

    Returns
    -------
    dict of pandas.DataFrame
        ``mc`` matched comparisons, ``bm`` the a/b subset, ``w`` per-acquisition counts
        (1,990 rows plus the 16 participant-15 CapScan FL_public rows), ``mem``, ``fu``,
        ``fa``, ``hs``, ``gr`` and ``inv``.
    """
    sd = Path(sd)
    mc = pd.read_csv(sd / "benchmark_matched_per_acquisition.tsv", sep="\t")
    bm = mc[mc.arm.eq('FL_public') & mc.comparator.eq('MG') & mc.population.isin(POPS_AB)]
    w = pd.read_csv(sd / "per_acquisition_length_counts.tsv", sep="\t")
    n16 = pd.read_csv(sd / "per_acquisition_length_counts_participant15.tsv", sep="\t")
    assert len(n16) == 16 and set(n16.cohort) == {'CAPSCAN'} and set(n16.arm) == {'FL_public'}
    assert not set(n16['sample']) & set(w['sample'])
    w = pd.concat([w, n16[[c for c in n16.columns if c in w.columns]]], ignore_index=True)
    assert len(w[(w.cohort == 'CAPSCAN') & (w.arm == 'FL_public')]) == CAPSCAN_ROSTER
    mem = pd.read_csv(sd / "campi_member_counts.tsv", sep="\t")
    fu = pd.read_csv(sd / "campi_fusion_pooled_q01.tsv", sep="\t")
    fa = pd.read_csv(sd / "campi_fusion_per_acquisition_q01.tsv", sep="\t")
    hs = pd.read_csv(sd / "healthy_strategies_per_acquisition.tsv", sep="\t")
    gr = pd.read_csv(sd / "grouping_study_groups.tsv", sep="\t")
    inv = pd.read_csv(sd / "search_inventory.tsv", sep="\t")
    return dict(mc=mc, bm=bm, w=w, mem=mem, fu=fu, fa=fa, hs=hs, gr=gr, inv=inv)


# ---------------------------------------------------------------- panels (drawing moved from build_sfig04_style.py, v30)
def _panel_ab(fig, gs, k, d, rng, rec):
    bm = d['bm']
    col, yl, letter = [('gain_pct', 'Peptide gain over MG (%)', 'a'),
                       ('comparator_identity_recovery_pct', 'MG identities recovered (%)', 'b')][k]
    ax = fig.add_subplot(gs[0, k]); fl.style(ax)
    for ci, co in enumerate(['PREDICT', 'CAPSCAN']):
        for pi, pop in enumerate(POPS_AB):
            sel = bm[(bm.cohort == co) & (bm.population == pop)]; v = sel[col].values; x = ci * 4 + pi
            vv = np.clip(v, None, GAIN_CLIP) if col == 'gain_pct' else v
            jx = x + rng.uniform(-JITTER_HALF_WIDTH, JITTER_HALF_WIDTH, len(v))
            ax.scatter(jx, vv, s=2, color=C[pop], alpha=.55, lw=0)
            for smp, xx, y0, y1 in zip(sel['sample'] if 'sample' in sel else [''] * len(v), jx, v, vv):
                rec(letter, cohort=co, population=pop, item='acquisition', sample=smp, x=xx, value=y0, plotted=y1, label='')
            if col == 'gain_pct' and (v > GAIN_CLIP).sum():
                ax.text(x, 262, f'{int((v > GAIN_CLIP).sum())}\n>{GAIN_CLIP}', ha='center', fontsize=5.5)
                rec(letter, cohort=co, population=pop, item='n_above_250', sample='', x=x, value=int((v > GAIN_CLIP).sum()), plotted='', label=f'{int((v > GAIN_CLIP).sum())} >{GAIN_CLIP}')
            ax.plot([x - .35, x + .35], [np.median(v)] * 2, color=fl.INK, lw=1.2)
            ax.annotate(f'{np.median(v):.1f}', (x, np.median(v)), xytext=(0, 1.6), textcoords='offset points', fontsize=5.5, va='bottom', ha='center',
                        bbox=dict(boxstyle='square,pad=0.08', fc='white', ec='none', alpha=.85), zorder=5)
            rec(letter, cohort=co, population=pop, item='median', sample='', x=x, value=float(np.median(v)), plotted=float(np.median(v)), label=f'{np.median(v):.1f}')
    ax.set_xticks([1, 5], ['PREDICT (n = 296)', 'CapScan (n = 282)']); ax.set_ylabel(yl)
    if k == 0: ax.axhline(0, color=fl.GREY, lw=.5); ax.set_ylim(-60, 280)
    return ax


def panel_a(fig, gs, d, rng, rec):
    """Panel a: per-acquisition peptide gain of public FastaLake over the matched metagenome database.

    Parameters
    ----------
    fig, gs : matplotlib Figure and GridSpec of the composite.
    d : dict from ``load``.
    rng : numpy.random.Generator shared with panel b (seed ``JITTER_SEED``).
    rec : callable that records each plotted value.

    Returns
    -------
    matplotlib.axes.Axes
    """
    return _panel_ab(fig, gs, 0, d, rng, rec)


def panel_b(fig, gs, d, rng, rec):
    """Panel b: recovery of metagenome-search peptide identities; same layout and rng as a."""
    return _panel_ab(fig, gs, 1, d, rng, rec)


def panel_c(fig, gs, d, rec):
    """Panel c: median change in accepted peptides vs 5-50 after filtering to 7-50 or re-searching at 7-50.

    Returns
    -------
    matplotlib.axes.Axes
    """
    w = d['w']
    ax = fig.add_subplot(gs[1, 0]); fl.style(ax)
    for i, (co, arm, lab) in enumerate(C_ARMS):
        g = w[(w.cohort == co) & (w.arm == arm)]
        a, b = g.nested7_vs_base5.median(), g.fresh7_vs_base5.median()
        ax.barh(i + .18, a, height=.34, color=C['nested7']); ax.barh(i - .18, b, height=.34, color=C['fresh7'])
        rec('c', cohort=co, population='nested7', item=lab, sample='', x=i + .18, value=a, plotted=a, label='')
        rec('c', cohort=co, population='fresh7', item=lab, sample='', x=i - .18, value=b, plotted=b, label='')
    ax.set_yticks(range(len(C_ARMS)), [a[2] for a in C_ARMS]); ax.invert_yaxis(); ax.axvline(0, color=fl.GREY, lw=.5)
    ax.set_xlabel('Median change in peptides vs 5–50 (%)')
    return ax


def panel_d(fig, gs, d, rec):
    """Panel d: pooled member-specific peptides of the three least represented CAMPI members (log axis, +0.5).

    Returns
    -------
    matplotlib.axes.Axes
    """
    mem = d['mem']
    ax = fig.add_subplot(gs[1, 1]); fl.style(ax)
    sel = [('FL', 'FastaLake'), ('MG', 'Metagenome'), ('REF', 'Reference')]
    for mi, m in enumerate(D_MEMBERS):
        for si, (arm, lab) in enumerate(sel):
            for pi, pop in enumerate(['baseline5', 'fresh7', 'nested9']):
                v = int(mem[(mem.arm == arm) & (mem.member == m) & (mem.population == pop)].pooled_specific.iloc[0])
                x = mi * 4 + si + (pi - 1) * .27; ax.bar(x, v + LOG_OFFSET, width=.26, color=C[pop], bottom=0)
                rec('d', cohort='CAMPI', population=pop, item=f'{m} | {lab}', sample='', x=x, value=v, plotted=v + LOG_OFFSET, label='')
    ax.set_yscale('log'); ax.set_ylabel('Member-specific peptides\n(pooled, +0.5)'); ax.set_ylim(.3, 1000)
    ax.tick_params(axis='y', labelsize=fl.LOG_TICK_PT)        # mathtext exponents at 0.7 x 7.2 = 5.04 pt
    ax.set_xticks([m * 4 + 1 for m in range(3)], ['L. plantarum', 'C. butyricum', 'B. longum'])
    sec = ax.secondary_xaxis('bottom')
    sec.set_xticks([m * 4 + s for m in range(3) for s in range(3)], ['FL', 'MG', 'Reference'] * 3, fontsize=5.5, rotation=270, va='top', ha='center')
    sec.tick_params(length=0, pad=2)
    ax.tick_params(axis='x', length=0, pad=31, labelsize=7)
    return ax


def panel_e(fig, gs, d, rec):
    """Panel e: CAMPI fusion-control diagnostic at q<=0.01, pooled over 15 acquisitions per seed.

    Returns
    -------
    matplotlib.axes.Axes
    """
    fu = d['fu']
    ax = fig.add_subplot(gs[2, 0]); fl.style(ax)
    for si, seed in enumerate(sorted(fu.seed.unique())):
        order = [('baseline', 5, 'baseline5'), ('baseline', 7, 'nested7'), ('fresh7', 7, 'fresh7'), ('baseline', 9, 'nested9'), ('fresh7', 9, 'fresh7_nested9')]
        for pi, (s, lo, pop) in enumerate(order):
            r = fu[(fu.seed == seed) & (fu.search == s) & (fu.min_residues == lo)].iloc[0]; x = si * 6 + pi
            ax.bar(x, r.pooled_diag_pct, color=C[pop]); ax.text(x, r.pooled_diag_pct, f'{r.pooled_diag_pct:.2f}', ha='center', va='bottom', fontsize=5.5)
            rec('e', cohort='CAMPI', population=pop, item=f'seed {seed}', sample='', x=x, value=r.pooled_diag_pct, plotted=r.pooled_diag_pct, label=f'{r.pooled_diag_pct:.2f}')
    ax.set_yticks(np.arange(0, 2.6, .5)); ax.set_xticks([2, 8], ['seed 20260907', 'seed 20260908']); ax.set_ylabel('Pooled fusion diagnostic (%)'); ax.axhline(FUSION_REF_LINE, color=fl.GREY, lw=.5, ls='--')
    return ax


def panel_f(fig, gs, d, rec):
    """Panel f: healthy-stool common-FASTA and joint searches relative to each acquisition's own database (58).

    Returns
    -------
    matplotlib.axes.Axes
    """
    hs = d['hs']
    ax = fig.add_subplot(gs[2, 1]); fl.style(ax)
    for si, strat in enumerate(['common', 'joint']):
        for pi, pop in enumerate(['baseline5', 'fresh7']):
            g = hs[(hs.strategy == strat) & (hs.population == pop)]
            x = si * 3 + pi; m1, m2 = g.count_change_pct.median(), g.own_identity_recovery_pct.median()
            ax.bar(x, m1, color=C[pop]); ax.text(x, m1, f'{m1:.1f}%\n{m2:.1f}% kept', ha='center', va='bottom', fontsize=5.5)
            rec('f', cohort='HSTOOL', population=pop, item=strat, sample='', x=x, value=m1, plotted=m1, label=f'{m1:.1f}% / {m2:.1f}% kept')
    ax.set_xticks([.5, 3.5], ['Common FASTA', 'Joint search']); ax.set_ylabel('Median peptide change\nvs own FASTA (%)'); ax.set_ylim(0, 45)
    return ax


def panel_g(fig, gs, d, rec):
    """Panel g: change in study-level protein groups, new 7-50 search vs 5-50 baseline (CapScan on 347).

    Returns
    -------
    matplotlib.axes.Axes
    """
    gr = d['gr'][['task', 'baseline', 'fresh7']].copy()
    ax = fig.add_subplot(gs[3, :]); fl.style(ax)
    gr['change_pct'] = 100 * (gr.fresh7 / gr.baseline - 1)
    names = {'HSTOOL_noclust': 'Healthy (58)', 'PREDICT_noclust': 'PREDICT (306)', 'CAPSCAN_noclust': f'CapScan ({CAPSCAN_ROSTER})', 'CAMPI_FL_UNIPROT': 'CAMPI\nFL',
             'CAMPI_FL_REF': 'CAMPI\nFL + Reference', 'CAMPI_FL_REF_MG': 'CAMPI\nFL + Reference\n+ MG', 'CAMPI_REF': 'CAMPI\nReference', 'CAMPI_MG': 'CAMPI\nMG'}
    gr['label'] = gr.task.map(names); gr = gr.set_index('task').loc[list(names)].reset_index()
    ax.bar(range(len(gr)), gr.change_pct, color=C['fresh7']); ax.set_xticks(range(len(gr)), gr.label, rotation=0)
    ax.tick_params(axis='x', labelsize=6)   # 6 pt as in v29; at 7 pt neighbouring labels touch
    ax.axhline(0, color=fl.GREY, lw=.5); ax.set_ylim(-3, 3); ax.set_yticks(range(-3, 4)); ax.set_ylabel('Study groups, 7–50 vs 5–50 (%)')
    for i, v in enumerate(gr.change_pct):
        ax.text(i, v + (.15 if v >= 0 else -.15), f'{v:+.1f}', ha='center', va='bottom' if v >= 0 else 'top', fontsize=5.5)
        rec('g', cohort=gr.task[i], population='fresh7 vs baseline5', item=gr.label[i].replace('\n', ' '), sample='', x=i, value=v, plotted=v, label=f'{v:+.1f}')
    return ax


def align_letters(fig, items, w_mm):
    """Draw each panel letter with ``fl.panel`` and move it to ``x_mm`` from the figure's left edge.

    Moved unchanged from ``_sf456_common.align_letters`` (the v30 restyle helper).

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    items : list of (Axes, str, float)
        Axes, letter and x position in mm; letters of one column share one x.
    w_mm : float
        Figure width in mm.
    """
    for ax, letter, x_mm in items:
        fl.panel(ax, letter)
        pos = ax.get_position()
        t = ax._left_title
        t.set_x((x_mm / w_mm - pos.x0) / pos.width)
        t.set_ha('left')


def _draw(sd):
    """Draw the v30 composite exactly as ``build_sfig04_style.py`` did; return (fig, axes dict, plotted rows)."""
    init()
    d = load(Path(sd))
    a = d['gr'].set_index('task')
    assert a.loc['CAPSCAN_noclust', 'acquisitions_fresh7'] == CAPSCAN_ROSTER
    rows = []

    def rec(panel, **kw):
        rows.append(dict(panel=panel, **kw))

    fig = plt.figure(figsize=fl.mm(PLACED_W_MM, PLACED_H_MM))
    gs = fig.add_gridspec(4, 2, hspace=.58, wspace=.42, top=.93, bottom=.075, left=.115, right=.985)
    rng = np.random.default_rng(JITTER_SEED)
    axes = {}
    axes['a'] = panel_a(fig, gs, d, rng, rec)
    axes['b'] = panel_b(fig, gs, d, rng, rec)
    axes['c'] = panel_c(fig, gs, d, rec)
    axes['d'] = panel_d(fig, gs, d, rec)
    axes['e'] = panel_e(fig, gs, d, rec)
    axes['f'] = panel_f(fig, gs, d, rec)
    axes['g'] = panel_g(fig, gs, d, rec)
    handles = [plt.Rectangle((0, 0), 1, 1, color=C[p]) for p in C]
    fig.legend(handles, [LAB[p] for p in C], loc='upper center', ncol=5, frameon=False, fontsize=6.5, bbox_to_anchor=(.53, .998), handlelength=1.4, columnspacing=1.2)
    align_letters(fig, [(axes['a'], 'a', LETTER_X_MM[0]), (axes['b'], 'b', LETTER_X_MM[1]), (axes['c'], 'c', LETTER_X_MM[0]), (axes['d'], 'd', LETTER_X_MM[1]),
                        (axes['e'], 'e', LETTER_X_MM[0]), (axes['f'], 'f', LETTER_X_MM[1]), (axes['g'], 'g', LETTER_X_MM[0])], PLACED_W_MM)
    n_geq = sum(any(c in tx.get_text() for c in '≥≤') for tx in fig.findobj(mtext.Text))
    assert n_geq == 0, n_geq    # Soehne has no >= or <= glyph
    return fig, axes, rows


def panel_boxes(sd=SD_DEFAULT, pad_px=6):
    """Pixel boxes (at ``PNG_DPI``) of each panel in the composite: axes, ticks, labels and letter.

    Returns
    -------
    pandas.DataFrame
        Columns ``panel, x0, y0, x1, y1`` (PIL crop convention, origin top left).
    """
    fig, axes, _ = _draw(sd)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    H = fig.get_figheight() * PNG_DPI
    W = fig.get_figwidth() * PNG_DPI
    out = []
    for p, ax in axes.items():
        bbs = [ax.get_tightbbox(r)] + [c.get_tightbbox(r) for c in ax.child_axes]
        x0 = min(b.x0 for b in bbs); x1 = max(b.x1 for b in bbs); y0 = min(b.y0 for b in bbs); y1 = max(b.y1 for b in bbs)
        k = PNG_DPI / fig.dpi
        out.append(dict(panel=p, x0=max(0, int(np.floor(x0 * k)) - pad_px), y0=max(0, int(np.floor(H - y1 * k)) - pad_px),
                        x1=min(int(round(W)), int(np.ceil(x1 * k)) + pad_px), y1=min(int(round(H)), int(np.ceil(H - y0 * k)) + pad_px)))
    plt.close(fig)
    return pd.DataFrame(out)


@lru_cache(maxsize=None)
def build_figure(sd, out):
    """Draw, save (``fl.save``), slant and rasterise the composite, as the v30 build did.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig04``.
    out : pathlib.Path
        ``build/SFig04``.

    Returns
    -------
    dict
        Paths of the final PDF, the 300 dpi PNG and the plotted-values table, and the
        build checks (smallest text size, number of slanted labels).
    """
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    fig, axes, rows = _draw(sd)
    raw = out / 'SFig04_raw'
    fl.save(fig, raw, dpi=PNG_DPI, formats=('pdf',))
    final = out / 'SFig04.pdf'
    n_sl = slant(raw.with_suffix('.pdf'), final, TAXON_LABELS)
    assert n_sl == 3, n_sl
    pv = out / 'SFig04_plotted_values.tsv'
    pd.DataFrame(rows).to_csv(pv, sep='\t', index=False)
    mt, k = min_text_pt(final, PLACED_W_MM)
    assert abs(k - 1) < 1e-6 and mt >= MIN_TEXT_PT, (mt, k)
    png = out / 'SFig04.png'; px = render_png(final, png, PNG_DPI)
    return dict(pdf=final, png=png, plotted_values=pv, png_px=px, smallest_text_pt=mt, slanted_labels=n_sl)


def save_panel(panel, sd, out):
    """Crop one panel from the composite render and save it as ``out/panel_<x>.png``.

    The crop box is computed by ``panel_boxes`` from the drawn axes.

    Returns
    -------
    pathlib.Path
    """
    from PIL import Image
    b = build_figure(Path(sd), Path(out))
    bx = _boxes(Path(sd)).set_index('panel').loc[panel]
    p = Path(out) / f'panel_{panel}.png'
    Image.open(b['png']).crop((int(bx.x0), int(bx.y0), int(bx.x1), int(bx.y1))).save(p)
    return p


@lru_cache(maxsize=None)
def _boxes(sd):
    return panel_boxes(sd)


def plotted(panel, out):
    """Return the plotted values of one panel (rows of the plotted-values table)."""
    pv = pd.read_csv(Path(out) / 'SFig04_plotted_values.tsv', sep='\t')
    return pv[pv.panel == panel]


# ---------------------------------------------------------------- statistics (descriptive)
_WORDS = {0: 'zero', 1: 'one', 2: 'two', 3: 'three', 4: 'four', 5: 'five', 6: 'six', 7: 'seven', 8: 'eight',
          9: 'nine', 10: 'ten', 11: 'eleven', 12: 'twelve'}


def stat_search_count(d):
    """Number of repeated searches behind the legend count (unit: search, one per acquisition and database).

    Returns
    -------
    dict
        ``total`` (all searches with a passing receipt, v28f roster), ``lfq_repaired`` (the
        five repaired PREDICT searches with LFQ), ``standard`` (total minus those five),
        ``total_v28d`` (the M09 searches alone, 331 roster).
    """
    inv = d['inv']
    assert inv.receipt_pass.all()
    lfq5 = int((inv.arm == 'FL_cohortMG_lfq').sum())
    return dict(total=len(inv), lfq_repaired=lfq5, lfq_repaired_word=_WORDS[lfq5], standard=len(inv) - lfq5,
                total_v28d=int((inv.run == 'M09_331').sum()), standard_v28d=int((inv.run == 'M09_331').sum()) - lfq5)


def stat_ab_matched(d, cohort, population):
    """Median paired gain and identity recovery of public FastaLake over the matched metagenome database.

    Statistical unit: acquisition; n = 296 PREDICT or 282 CapScan matched pairs.

    Returns
    -------
    dict
        ``n``, ``median_gain_pct`` (median over acquisitions of 100 (FL/MG - 1)),
        ``median_recovery_pct``, ``n_more_fl`` (acquisitions with more FL than MG
        peptides), ``n_above_clip``.
    """
    s = d['mc'][(d['mc'].arm == 'FL_public') & (d['mc'].comparator == 'MG') & (d['mc'].cohort == cohort) & (d['mc'].population == population)]
    return dict(n=len(s), median_gain_pct=float(np.median(s.gain_pct)), median_recovery_pct=float(np.median(s.comparator_identity_recovery_pct)),
                n_more_fl=int((s.arm_peptides > s.comparator_peptides).sum()), n_above_clip=int((s.gain_pct > GAIN_CLIP).sum()))


def stat_gain_5_8(d, cohort):
    """Median per-acquisition gain restricted to peptides of 5-8 residues.

    The 5-8 count of each arm is the 5-50 count minus the nested 9-50 count of the same
    search (production q values). Statistical unit: acquisition; n = 296 or 282.

    Returns
    -------
    dict
        ``n``, ``median_gain_pct``.
    """
    b = d['mc'][(d['mc'].arm == 'FL_public') & (d['mc'].comparator == 'MG') & (d['mc'].cohort == cohort)]
    p = b.pivot(index='sample', columns='population', values=['arm_peptides', 'comparator_peptides'])
    a58 = p['arm_peptides']['baseline5'] - p['arm_peptides']['nested9']
    c58 = p['comparator_peptides']['baseline5'] - p['comparator_peptides']['nested9']
    g = 100 * (a58 / c58 - 1)
    return dict(n=len(g), median_gain_pct=float(g.median()))


def stat_c_median_change(d, cohort, arm):
    """Median per-acquisition change in accepted peptides vs 5-50 (panel c).

    Statistical unit: acquisition (n per database as in ``n``).

    Returns
    -------
    dict
        ``n``, ``filtered_to_7`` and ``new_search_7`` medians (%).
    """
    g = d['w'][(d['w'].cohort == cohort) & (d['w'].arm == arm)]
    return dict(n=len(g), filtered_to_7=float(g.nested7_vs_base5.median()), new_search_7=float(g.fresh7_vs_base5.median()))


def stat_d_members(d, member, arm, population):
    """Pooled member-specific peptides of one CAMPI member (unit: peptide, pooled over 15 acquisitions)."""
    m = d['mem']
    v = int(m[(m.arm == arm) & (m.member == member) & (m.population == population)].pooled_specific.iloc[0])
    return dict(count=v, word=_WORDS.get(v, str(v)))


def stat_d_fold_ref_over_fl(d, population='baseline5'):
    """Fold difference reference / FastaLake in pooled member-specific peptides, three least represented members.

    Statistical unit: peptide, pooled over 15 CAMPI acquisitions.

    Returns
    -------
    dict
        ``by_member`` ratios, ``min`` and ``max``.
    """
    r = {m: stat_d_members(d, m, 'REF', population)['count'] / stat_d_members(d, m, 'FL', population)['count'] for m in D_MEMBERS}
    return dict(by_member=r, min=min(r.values()), max=max(r.values()))


def stat_e_pooled_diagnostic(d):
    """Pooled fusion diagnostic per seed and population, recomputed from the per-acquisition rows.

    Pooled diagnostic = 100 (sum paired_numerator + sum n_unpaired_accepted) / sum n_canonical
    over the 15 acquisitions of a seed, at q <= 0.01. Statistical unit: acquisition pooled
    within control seed; n = 15 per seed, 2 seeds.

    Returns
    -------
    pandas.DataFrame
        One row per seed, search and length floor, with the recomputed value, the value in
        the drawn table and their absolute difference.
    """
    fa = d['fa'][np.isclose(d['fa'].q_threshold, FUSION_Q)]
    g = fa.groupby(['seed', 'search', 'min_residues']).agg(acq=('sample', 'nunique'), n=('n_canonical', 'sum'),
                                                          num=('paired_numerator', 'sum'), unp=('n_unpaired_accepted', 'sum')).reset_index()
    g['recomputed_pct'] = 100 * (g.num + g.unp) / g.n
    m = g.merge(d['fu'], on=['seed', 'search', 'min_residues'], suffixes=('', '_drawn'))
    m['abs_diff'] = (m.recomputed_pct - m.pooled_diag_pct).abs()
    return m[['seed', 'search', 'min_residues', 'acq', 'n', 'num', 'unp', 'recomputed_pct', 'pooled_diag_pct', 'abs_diff']]


def stat_f_healthy(d, strategy, population):
    """Median peptide change and identities kept vs own database (unit: acquisition, n = 58)."""
    g = d['hs'][(d['hs'].strategy == strategy) & (d['hs'].population == population)]
    return dict(n=len(g), median_change_pct=float(g.count_change_pct.median()), median_kept_pct=float(g.own_identity_recovery_pct.median()))


def stat_g_change(d, task):
    """Change in study-level protein groups, new >=7 search vs 5-50 (unit: protein group, one study-level grouping)."""
    r = d['gr'].set_index('task').loc[task]
    return dict(baseline=int(r.baseline), fresh7=int(r.fresh7), change_pct=float(100 * (r.fresh7 / r.baseline - 1)),
                acquisitions=int(r.acquisitions_fresh7))


def numbers(sd=SD_DEFAULT):
    """Every number of NUMBERS.tsv, recomputed, in NUMBERS.tsv order.

    Keys name the v28d line; the v30 paragraph of each row is in the ``v30_para`` column.

    Returns
    -------
    list of (str, object)
        ``(key, value)``; the key names the NUMBERS.tsv row.
    """
    d = load(Path(sd))
    sc = stat_search_count(d)
    out = [('L83 L. plantarum 5-50', stat_d_members(d, 'Lactobacillus plantarum', 'FL', 'baseline5')['word']),
           ('L83 L. plantarum 9-50', stat_d_members(d, 'Lactobacillus plantarum', 'FL', 'nested9')['word'])]
    f = stat_d_fold_ref_over_fl(d)
    out += [('L83 fold min', f['min']), ('L83 fold max', f['max'])]
    out += [('L239 searches', sc['total']),
            ('L421 searches', sc['total']), ('L421 standard', sc['standard']), ('L421 LFQ repaired', sc['lfq_repaired_word'])]
    out += [('L421 PREDICT n', stat_ab_matched(d, 'PREDICT', 'baseline5')['n']), ('L421 CapScan n', stat_ab_matched(d, 'CAPSCAN', 'baseline5')['n']),
            ('L421 clip', GAIN_CLIP), ('L421 log offset', LOG_OFFSET),
            ('L421 fusion acquisitions per seed', int(stat_e_pooled_diagnostic(d).acq.unique()[0])), ('L421 fusion q', FUSION_Q),
            ('L421 dashed line', FUSION_REF_LINE), ('L421 healthy n', stat_f_healthy(d, 'common', 'baseline5')['n'])]
    p7, c7 = stat_ab_matched(d, 'PREDICT', 'fresh7'), stat_ab_matched(d, 'CAPSCAN', 'fresh7')
    p9, c9 = stat_ab_matched(d, 'PREDICT', 'nested9'), stat_ab_matched(d, 'CAPSCAN', 'nested9')
    out += [('L777 PREDICT gain fresh7', p7['median_gain_pct']), ('L777 CapScan gain fresh7', c7['median_gain_pct']),
            ('L777 PREDICT recovery fresh7', p7['median_recovery_pct']), ('L777 CapScan recovery fresh7', c7['median_recovery_pct']),
            ('L777 PREDICT gain nested9', p9['median_gain_pct']), ('L777 CapScan gain nested9', c9['median_gain_pct']),
            ('L777 PREDICT more FL nested9', p9['n_more_fl']), ('L777 PREDICT n', p9['n']),
            ('L777 CapScan more FL nested9', c9['n_more_fl']), ('L777 CapScan n', c9['n']),
            ('L777 PREDICT gain 5-8', stat_gain_5_8(d, 'PREDICT')['median_gain_pct']), ('L777 CapScan gain 5-8', stat_gain_5_8(d, 'CAPSCAN')['median_gain_pct'])]
    return out
