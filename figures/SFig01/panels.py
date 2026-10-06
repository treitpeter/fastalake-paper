"""Supplementary Figure 1: selected database sizes in the context of the full sequence reservoirs.

One unlettered panel (called ``a`` in this repository): FASTA entries of every
acquisition-specific database (726 acquisitions) per cohort, cohort medians, the full
cohort reservoir sizes, and the 20,431-entry human reference proteome, on a log axis.

Source data
-----------
source_data/Fig1/acquisition_sizes_726.tsv
    One row per acquisition: cohort, sample, reservoir_entries, entries, ... This table is
    shared with Figure 1b and lives under source_data/Fig1/ (never duplicated).
source_data/SFig01/human_swissprot_canonical_2026_03.fasta.gz
    The reviewed human Swiss-Prot canonical proteome used as the size reference.

The plotting code is the SF1 block of the original producer (build_evidence.py, lines
8-27, as copied into build_sf1.py for the 726-row input) and the raster step of
reexport.py (as copied into render.py). It is moved unchanged except that input and
output paths are function arguments and that literal values are named constants below
(same values). Matplotlib default style is used, as in the producer; the house style
set by ``fl_style.init()`` is switched off inside ``panel_a``.

v30: ``panel_a_v30`` redraws the same values in the house style (Söhne, ink, Okabe-Ito by
meaning, 180 x 105.2 mm) and ``values_unchanged_v30`` proves the plotted arrays identical to
``panel_a``; moved unchanged from the v30 rebuild ``build_sfig01_style.py`` together with the
restyle helpers of ``restyle_common.py`` (fingerprint, compare, texts, capture_saved,
align_letters, draw_rel_label, render_png).
"""

import gzip
import hashlib
import io
import math
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---- choices (named constants) ------------------------------------------------------
HUMAN_REFERENCE_ENTRIES = 20431   # Methods, v28d line 673 (Swiss-Prot human canonical, 2026_03)
SIZES_SHA256 = "fd9ade3bb4ebd2ee6c129f69819b9697cf20abd7c10c77a12644bef25396e386"  # 726-row table
HUMAN_REFERENCE_SHA256 = "b205ed9970fee62413e2f192c6fb82daab86e4a5897ef521f907ecbfcd7b3e63"
JITTER_SEED = 20260928            # producer build_evidence.py line 15; not stated in the text
JITTER_HALF_WIDTH = 0.14          # producer line 17; horizontal jitter of points, display only
COHORT_ORDER = ["CAMPI", "HSTOOL", "PREDICT", "CapScan"]   # producer line 13; x-axis order
COHORT_NAMES = {"CAMPI": "CAMPI", "HSTOOL": "Healthy stool", "PREDICT": "PREDICT",
                "CapScan": "CapScan"}
Y_LIMITS = (150, 1.2e9)           # producer line 23; log axis (Methods, v28d line 675)
FIGSIZE_IN = (7.2, 4.2)           # producer line 14
PLACED_WIDTH_MM = 180.34          # placed width in the manuscript docx (reexport.py)
PLACED_HEIGHT_MM = 105.20         # placed height in the manuscript docx (reexport.py)
RASTER_DPI = 400                  # reexport.py TARGET_DPI (effective dpi at placed width)

_DEFAULT_RC = {k: v for k, v in mpl.rcParamsDefault.items() if k != "backend"}


def check_sha256(path, expected):
    """Raise if the sha256 of path differs from expected.

    Parameters
    ----------
    path : pathlib.Path
        File to check.
    expected : str
        Hex sha256.

    Returns
    -------
    str
        The sha256.
    """
    s = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if s != expected:
        raise ValueError(f"{Path(path).name}: sha256 {s} differs from expected {expected}")
    return s


def load_sizes(path):
    """Read the per-acquisition size table.

    Parameters
    ----------
    path : pathlib.Path
        ``source_data/SFig01/acquisition_sizes_726.tsv``.

    Returns
    -------
    pandas.DataFrame
        One row per acquisition (726 rows).
    """
    return pd.read_csv(path, sep="\t")


def panel_a(t, out_pdf):
    """Draw the figure exactly as the original producer and save it as a PDF.

    Parameters
    ----------
    t : pandas.DataFrame
        Output of ``load_sizes``.
    out_pdf : pathlib.Path
        PDF to write.

    Returns
    -------
    dict
        Per cohort ``(n, median, min, max)`` of FASTA entries, as the producer printed.
    """
    O = Path(out_pdf)
    with mpl.rc_context(_DEFAULT_RC):
        # ---- moved from build_sf1.py (SF1 block of build_evidence.py), unchanged logic ----
        assert (t.entries < HUMAN_REFERENCE_ENTRIES).all()
        order = COHORT_ORDER
        names = COHORT_NAMES
        fig, ax = plt.subplots(figsize=FIGSIZE_IN)
        rng = np.random.default_rng(JITTER_SEED)
        for i, c in enumerate(order):
            z = t[t.cohort == c]
            ax.scatter(i + rng.uniform(-JITTER_HALF_WIDTH, JITTER_HALF_WIDTH, len(z)), z.entries,
                       s=8, alpha=.35, color='#1675a9', linewidths=0)
            ax.plot([i - .2, i + .2], [z.entries.median()] * 2, color='#093e5a', lw=2)
            pool = int(z.reservoir_entries.iloc[0])
            assert z.reservoir_entries.nunique() == 1
            ax.scatter([i], [pool], marker='D', s=40, color='#c67621', zorder=4)
            ax.text(i, pool * 1.25, f'{pool/1e6:.1f} M', ha='center', fontsize=8)
        ax.axhline(HUMAN_REFERENCE_ENTRIES, ls='--', c='#555', lw=1,
                   label=f'Human reference: {HUMAN_REFERENCE_ENTRIES:,} entries')
        ax.set_yscale('log')
        ax.set_ylim(*Y_LIMITS)
        ax.set_xticks(range(4), [f'{names[c]}\nn={(t.cohort==c).sum()}' for c in order])
        ax.set_ylabel('Protein entries (log scale)')
        ax.spines[['top', 'right']].set_visible(False)
        ax.scatter([], [], s=20, color='#1675a9', label='Selected acquisition FASTAs (median line)')
        ax.scatter([], [], s=35, color='#c67621', marker='D', label='Full cohort reservoir')
        ax.legend(frameon=False, fontsize=8, loc='center left')
        fig.tight_layout()
        fig.savefig(O)
        plt.close(fig)
    return {c: (int((t.cohort == c).sum()), float(t[t.cohort == c].entries.median()),
                int(t[t.cohort == c].entries.min()), int(t[t.cohort == c].entries.max()))
            for c in COHORT_ORDER}


def render_composite(pdf, out_png, wmm=PLACED_WIDTH_MM, hmm=PLACED_HEIGHT_MM, dpi=RASTER_DPI):
    """Rasterise the figure PDF at ``dpi`` effective at the placed width (reexport.py step).

    Parameters
    ----------
    pdf : pathlib.Path
        PDF written by ``panel_a``.
    out_png : pathlib.Path
        PNG to write (the image placed in the manuscript).
    wmm, hmm : float
        Placed width and height in millimetres.
    dpi : int
        Target effective resolution.

    Returns
    -------
    dict
        Pixel size, effective dpi and aspect deviation from the placed extent.
    """
    import fitz
    from PIL import Image

    w_in = wmm / 25.4
    h_in = hmm / 25.4
    width_px = math.ceil(w_in * dpi)
    pg = fitz.open(pdf)[0]
    zoom = width_px / pg.rect.width
    pix = pg.get_pixmap(matrix=fitz.Matrix(zoom, zoom), alpha=False)
    im = Image.open(io.BytesIO(pix.tobytes('png')))
    im.save(out_png, dpi=(dpi, dpi), optimize=True)
    return dict(px=f'{im.width}x{im.height}', eff_dpi=round(min(im.width / w_in, im.height / h_in), 1),
                aspect_delta=round(im.width / im.height / (w_in / h_in) - 1, 5))


def stat_counts(t):
    """Number of acquisitions in total and per cohort.

    Statistical unit: acquisition; n = 726 (15 CAMPI, 58 HSTOOL, 306 PREDICT, 347 CapScan).

    Parameters
    ----------
    t : pandas.DataFrame
        Output of ``load_sizes``.

    Returns
    -------
    dict
        ``{"total": int, "<cohort>": int, ...}``.
    """
    d = {"total": int(len(t))}
    d.update({c: int((t.cohort == c).sum()) for c in COHORT_ORDER})
    return d


def stat_entries(t):
    """Median, minimum and maximum FASTA entries per cohort, overall maximum, and the
    number of databases below the human reference size.

    Statistical unit: acquisition; n = 726. Descriptive only, no test.

    Parameters
    ----------
    t : pandas.DataFrame
        Output of ``load_sizes``.

    Returns
    -------
    dict
        ``median``, ``min``, ``max`` (per cohort), ``max_all``, ``n_below_human``.
    """
    return dict(
        median={c: float(t[t.cohort == c].entries.median()) for c in COHORT_ORDER},
        min={c: int(t[t.cohort == c].entries.min()) for c in COHORT_ORDER},
        max={c: int(t[t.cohort == c].entries.max()) for c in COHORT_ORDER},
        max_all=int(t.entries.max()),
        n_below_human=int((t.entries < HUMAN_REFERENCE_ENTRIES).sum()),
    )


def stat_reservoirs(t):
    """Full reservoir size per cohort (one value per cohort, constant within cohort).

    Statistical unit: sequence (reservoir entries); one reservoir per cohort.

    Parameters
    ----------
    t : pandas.DataFrame
        Output of ``load_sizes``.

    Returns
    -------
    dict
        ``{cohort: reservoir entries}``.
    """
    out = {}
    for c in COHORT_ORDER:
        r = t.loc[t.cohort == c, "reservoir_entries"]
        assert r.nunique() == 1, c
        out[c] = int(r.iloc[0])
    return out


def stat_human_reference(fasta_gz):
    """Count entries of the human reference proteome FASTA.

    Statistical unit: sequence; n = entries in the file.

    Parameters
    ----------
    fasta_gz : pathlib.Path
        ``source_data/SFig01/human_swissprot_canonical_2026_03.fasta.gz``.

    Returns
    -------
    dict
        ``entries`` (header lines), ``sha256_uncompressed`` and whether it equals the
        sha256 recorded by the upstream analysis.
    """
    h = hashlib.sha256()
    n = 0
    with gzip.open(fasta_gz, "rb") as fh:
        for line in fh:
            h.update(line)
            n += line.startswith(b">")
    s = h.hexdigest()
    return dict(entries=n, sha256_uncompressed=s, sha256_matches_upstream=s == HUMAN_REFERENCE_SHA256)


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
# v30: house-style rebuild of the same panel
# =====================================================================================
V30_PLACED_W_MM = 180.0           # v30 build: drawn at the placed size (v29 docx extent 180.3 x 105.2 mm)
V30_PLACED_H_MM = 105.2
V30_PNG_DPI = 300                 # v30 build: PNG rasterised at 300 dpi from the PDF


def panel_a_v30(t, out_stem):
    """Draw Supplementary Figure 1 in the house style (v30) and write PDF and 300-dpi PNG.

    Moved unchanged from the v30 rebuild ``build_sfig01_style.py`` (step 2, house-style
    redraw, and step 3, write). Same data, seed and draw order as ``panel_a``, so the jitter
    is identical; only typeface, colours, sizes and the placement of the reservoir value
    labels (beside the diamonds instead of on top of them) differ.

    Parameters
    ----------
    t : pandas.DataFrame
        Output of ``load_sizes``.
    out_stem : pathlib.Path
        Output path without suffix; ``<stem>.pdf`` and ``<stem>.png`` are written.

    Returns
    -------
    tuple
        ``(fig, ax, pdf_path, png_path, png_px)``; the figure is left open for the
        values-unchanged proof (``fingerprint``).
    """
    from style import fl_style as fl
    mpl.rcdefaults()
    fl.init()
    assert (t.entries < HUMAN_REFERENCE_ENTRIES).all()
    order = COHORT_ORDER
    names = COHORT_NAMES
    fig, ax = plt.subplots(figsize=fl.mm(V30_PLACED_W_MM, V30_PLACED_H_MM))
    fig.subplots_adjust(left=0.085, right=0.985, top=0.97, bottom=0.115)
    rng = np.random.default_rng(JITTER_SEED)            # producer seed and draw order -> identical jitter
    for i, c in enumerate(order):
        z = t[t.cohort == c]
        ax.scatter(i + rng.uniform(-JITTER_HALF_WIDTH, JITTER_HALF_WIDTH, len(z)), z.entries, s=6, alpha=.35, color=fl.FL, linewidths=0)
        ax.plot([i - .2, i + .2], [z.entries.median()] * 2, color=fl.INK, lw=1.5, solid_capstyle='butt')
        pool = int(z.reservoir_entries.iloc[0]); assert z.reservoir_entries.nunique() == 1
        ax.scatter([i], [pool], marker='D', s=22, color=fl.ORANGE, linewidths=0, zorder=4)
        # value label beside the diamond (producer put it above, where it overlapped the marker)
        ax.text(i + 0.06, pool, f'{pool/1e6:.1f} M', ha='left', va='center', fontsize=7, color=fl.INK)
    ax.axhline(HUMAN_REFERENCE_ENTRIES, ls='--', c=fl.REF, lw=1, label=f'Human reference: {HUMAN_REFERENCE_ENTRIES:,} entries')
    ax.set_yscale('log'); ax.set_ylim(*Y_LIMITS)
    ax.set_xticks(range(4), [f'{names[c]}\nn={(t.cohort == c).sum()}' for c in order])
    ax.set_ylabel('Protein entries (log scale)')
    fl.style(ax)
    ax.tick_params(axis='y', labelsize=fl.LOG_TICK_PT)   # exponents at 0.7 x 7.2 = 5.04 pt
    ax.tick_params(axis='y', which='minor', length=1.2, width=fl.MIN_LINE_PT)
    ax.scatter([], [], s=14, color=fl.FL, linewidths=0, label='Selected acquisition FASTAs (median line)')
    ax.scatter([], [], s=18, color=fl.ORANGE, marker='D', linewidths=0, label='Full cohort reservoir')
    ax.legend(frameon=False, fontsize=7, loc='center left', bbox_to_anchor=(0.0, 0.47), handlelength=2.2)
    stem = Path(out_stem)
    pdf, _ = fl.save(fig, stem, formats=('pdf',))
    png = stem.with_suffix('.png')
    px = render_png(pdf, png, V30_PNG_DPI)
    return fig, ax, Path(pdf), png, px


def values_unchanged_v30(t, out_stem):
    """Prove that the v30 redraw plots exactly the producer's values.

    Runs ``panel_a`` (the original producer) with ``savefig`` intercepted and the v30
    redraw, then compares every plotted array, axis limit and scale exactly, and the label
    strings (step 1 and 3 of ``build_sfig01_style.py``).

    Parameters
    ----------
    t : pandas.DataFrame
        Output of ``load_sizes``.
    out_stem : pathlib.Path
        Passed to ``panel_a_v30``.

    Returns
    -------
    dict
        ``values_identical`` (bool), per-panel report, ``label_strings_lost`` and
        ``label_strings_added``, ``png_px``.
    """
    mpl.rcdefaults()
    with capture_saved() as cap:
        panel_a(t, Path(out_stem).with_name('_unused.pdf'))
    (old_fig,) = cap.figs
    OLD = fingerprint([old_fig.axes[0]])
    (old_text,) = cap.texts
    fig, ax, pdf, png, px = panel_a_v30(t, out_stem)
    NEW = fingerprint([ax])
    new_text = texts(fig)
    ok, rep = compare(OLD, NEW, ['SF1'])
    plt.close(fig)
    return dict(values_identical=ok, panels=rep, label_strings_lost=sorted(set(old_text) - set(new_text)),
                label_strings_added=sorted(set(new_text) - set(old_text)), png_px=px)
