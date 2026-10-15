"""Supplementary Figure 3: CAMPI member support and identity recovery (panels a-d).

Panels
------
a, b  Pooled member-specific peptide counts in the FastaLake (FL) arm at 5-50 and 9-50
      residues, labelled with the count and the number of acquisitions (of 15) that
      contribute at least one such peptide.
c     Within-acquisition recovery of reference identities by each Sage search arm.
d     Within-MSFragger recovery of reference and MG identities by FastaLake.

Source data (``source_data/SFig03/``)
-------------------------------------
CAMPI_member_support_recount.tsv     panels a, b; Results v28d line 83
CAMPI_member_acquisition_counts.tsv  per-acquisition member counts (Methods v28d line 685)
CAMPI_reference_recovery.tsv         panel c
CAMPI_within_MSFragger.tsv           panel d; Results v28d line 89
CAMPI_membership_<sample>.tsv        accepted peptides of all five Sage arms with their
                                     eight-member reference membership (15 files); used to
                                     recount the three tables above
second_engine_between_engine_identities.tsv  MSFragger versus Sage identities (Source Data
                                     only, not drawn)
CAMPI_MSFragger_search_receipts.tsv  one row per MSFragger + Percolator search (45): engine
                                     versions, accepted counts, input hashes, Slurm job id

The drawing code is moved unchanged from the original producer (``sf3.py`` and
``common.py`` of the round-5 figure build, see PROVENANCE.md). Only input and output
paths were made repository-relative, the style initialisation reads the Söhne file from
``FASTALAKE_FONT_DIR``, and the plotted-values table names repository files instead of
lab paths.

v30: panels c,d draw the nested 9-50-residue rows (``PLOTTED_MIN_LENGTH_CD = 9``, the len950
rebuild ``panels_len950.py``); ``build_figure_v30`` and ``build_v30`` redraw the figure in the
house style and prove every plotted value identical to ``build_figure``; moved unchanged from
the v30 rebuild ``build_sfig03_style.py`` together with the restyle helpers of
``restyle_common.py``.
"""

import collections
import hashlib
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import font_manager as fm  # noqa: E402

# ---------------------------------------------------------------------------------------
# Choices. Every value here is quoted in README.md, Choices table.
# ---------------------------------------------------------------------------------------
MIN_LENGTHS = (5, 9)          # Methods, v28d line 685: accepted 5-50 and nested 9-50 residues
MAX_LENGTH = 50               # Methods, v28d line 685
PLOTTED_MIN_LENGTH_CD = 9     # panels c, d draw the nested 9-50 population (len950 rebuild; the v28f producer sf3.py used 5)
N_ACQUISITIONS = 15           # Methods, v28d line 253: all 15 CAMPI acquisitions
SPECIFIC_N_MEMBERS = 1        # member-specific = assigned to exactly one of the eight reference proteomes (legend, v28d line 417)
ARMS = ["FL", "FL_R", "FL_R_M", "REF", "MG"]   # legend, v28d line 417
MSF_COMPARATORS = ["REF", "MG"]                # panel d; Results v28d line 89
ARMLAB = {"FL": "FL", "FL_R": "FL +\nReference", "FL_R_M": "FL +\nReference\n+ MG",
          "REF": "Reference", "MG": "MG"}     # label edit FAR15 (producer sf3.py)
COL = ["#2878a3", "#df8f44", "#619b68", "#b35872", "#8277ac"]   # producer palette
PLACED_W_MM = 160.0           # placed width in the manuscript (producer sf3.py)
PLACED_H_MM = 140.0           # figure height (producer sf3.py)
BASE_FONT_PT = 7.0            # producer common.init(7.0)
MIN_TEXT_PT = 5.0             # smallest text at placed size (producer assertion)
SLANT_SHEAR = 0.20            # synthetic slant for taxon labels, label edit FAR16
N_SLANTED_LABELS = 16         # 8 taxon labels in each of panels a and b
RENDER_DPI = 400              # resolution of the embedded manuscript image
XLIM_AB = (0.8, 1e7)          # producer sf3.py
YLIM_CD = (0, 105)            # producer sf3.py
INK = "#4e4f51"               # producer common.py
FONT_FILES = ("soehne-buch.ttf", "soehne-buch.otf")   # producer used soehne-buch.ttf
MM = 72 / 25.4

#: Crop boxes (left, top, right, bottom) in pixels of the 400-dpi composite (2520 x 2205),
#: used only to cut per-panel images for the pixel checks. Each box is the panel's axes
#: tight bounding box (ticks, labels, title, value labels) plus a 12-pixel margin, clipped
#: to the image; computed once with ``panel_crop_boxes`` on the producer layout and fixed here.
CROP_BOXES = {'a': (41, 21, 1257, 1043), 'b': (1310, 21, 2520, 1043),
              'c': (262, 1158, 1225, 2110), 'd': (1532, 1158, 2494, 2054)}

TABLES = {
    "support": "CAMPI_member_support_recount.tsv",
    "acq_counts": "CAMPI_member_acquisition_counts.tsv",
    "recovery": "CAMPI_reference_recovery.tsv",
    "msfragger": "CAMPI_within_MSFragger.tsv",
    "between_engine": "second_engine_between_engine_identities.tsv",
    "receipts": "CAMPI_MSFragger_search_receipts.tsv",
}
SAMPLES = ["S01", "S02", "S03", "S04", "S05", "S06", "S07", "S08", "S09", "S10",
           "S11_Fraction1", "S11_Fraction2", "S11_Fraction3", "S11_Fraction4", "S12"]


def sha(p):
    """Return the hex sha256 of a file."""
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


# ---------------------------------------------------------------------------------------
# Style (moved from producer common.py; font path read from FASTALAKE_FONT_DIR)
# ---------------------------------------------------------------------------------------
_FONT = {"path": None, "name": None}


def font_file():
    """Return the Söhne Buch file in ``FASTALAKE_FONT_DIR`` (or ``FL_STYLE_DIR``), or None.

    Returns
    -------
    pathlib.Path or None
        ``soehne-buch.ttf`` if present (the file the producer used), else
        ``soehne-buch.otf``, else None.
    """
    d = os.environ.get("FASTALAKE_FONT_DIR") or os.environ.get("FL_STYLE_DIR")
    if not d:
        return None
    for name in FONT_FILES:
        p = Path(d) / name
        if p.is_file():
            return p
    return None


def init(base=BASE_FONT_PT):
    """Apply the producer rcParams (common.init) and register Söhne Buch if available.

    Parameters
    ----------
    base : float
        Base font size in points.

    Returns
    -------
    bool
        True if the Söhne face was found and registered.
    """
    plt.rcdefaults()   # the producer started from matplotlib defaults; undo any earlier style
    f = font_file()
    _FONT.update(path=None, name=None)
    if f is not None:
        fm.fontManager.addfont(str(f))
        _FONT.update(path=f, name=fm.FontProperties(fname=str(f)).get_name())
    family = _FONT["name"] or "DejaVu Sans"
    plt.rcParams.update({
        'font.family': family, 'font.size': base, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
        'mathtext.fontset': 'custom', 'mathtext.rm': family, 'mathtext.it': family, 'mathtext.bf': family,
        'mathtext.default': 'rm', 'axes.unicode_minus': True,
        'text.color': INK, 'axes.labelcolor': INK, 'axes.edgecolor': INK, 'xtick.color': INK, 'ytick.color': INK,
        'axes.linewidth': .6, 'xtick.major.width': .6, 'ytick.major.width': .6,
        'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
        'axes.spines.top': False, 'axes.spines.right': False,
        'axes.titleweight': 'normal', 'figure.titleweight': 'normal', 'font.weight': 'normal',
    })
    return f is not None


# ---------------------------------------------------------------------------------------
# Data
# ---------------------------------------------------------------------------------------
def load(sd):
    """Read the four Supplementary Figure 3 tables.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig03``.

    Returns
    -------
    dict of pandas.DataFrame
        Keys ``support``, ``acq_counts``, ``recovery``, ``msfragger``, ``between_engine``,
        ``receipts``.
    """
    return {k: pd.read_csv(Path(sd) / v, sep="\t") for k, v in TABLES.items()}


def short_names(names):
    """Abbreviate binomials as in the producer ('Bacteroides thetaiotaomicron' -> 'B. thetaiotaomicron')."""
    return [n.split()[0][0] + '. ' + ' '.join(n.split()[1:]) for n in names]


# ---------------------------------------------------------------------------------------
# Panels (drawing code moved from producer sf3.py)
# ---------------------------------------------------------------------------------------
def new_figure():
    """Create the 2 x 2 figure at placed size with the producer's margins.

    Returns
    -------
    (matplotlib.figure.Figure, numpy.ndarray of Axes)
    """
    fig, axes = plt.subplots(2, 2, figsize=(PLACED_W_MM / 25.4, PLACED_H_MM / 25.4)); axes = axes.ravel()
    fig.subplots_adjust(left=.17, right=.985, top=.955, bottom=.10, wspace=.62, hspace=.52)
    return fig, axes


def _panel_ab(ax, p, i, lo, rows, src):
    names = list(p.member.unique())
    short = short_names(names)
    sub = p[p.arm.eq('FL') & p.min_length.eq(lo)].set_index('member').loc[names]
    ax.barh(np.arange(8), sub.pooled_specific_peptides, color=COL[i])
    ax.set_yticks(np.arange(8), short); ax.set_xscale('log'); ax.set_xlabel('Pooled member-specific\npeptide count'); ax.invert_yaxis()
    ax.tick_params(axis='x', labelsize=7.5)          # mathtext exponents print at 0.7 x 7.5 = 5.25 pt
    for k, x in enumerate(sub.itertuples()):
        lab = f'{x.pooled_specific_peptides:,}; {x.acquisitions_with_support}/15'
        ax.text(max(x.pooled_specific_peptides, 1) * 1.08, k, lab, va='center', fontsize=6.5)
        rows.append(dict(panel='ab'[i], group=x.Index, arm='FL', min_length=lo, sample='', x=k, value=x.pooled_specific_peptides,
                         acquisitions_with_support=x.acquisitions_with_support, label=lab, source=src))
    ax.set_xlim(*XLIM_AB); ax.set_title('ab'[i] + '  ' + f'FastaLake, {lo}–50 residues', loc='left', pad=6, fontsize=8)


def panel_a(ax, t, rows):
    """Draw panel a: pooled FL member-specific peptides, 5-50 residues.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    t : dict of pandas.DataFrame
        Output of ``load``.
    rows : list
        Plotted values are appended here (one dict per bar).
    """
    _panel_ab(ax, t["support"], 0, MIN_LENGTHS[0], rows, TABLES["support"])


def panel_b(ax, t, rows):
    """Draw panel b: pooled FL member-specific peptides, nested 9-50 residues (see ``panel_a``)."""
    _panel_ab(ax, t["support"], 1, MIN_LENGTHS[1], rows, TABLES["support"])


def dots(ax, groups, ylabel, panel, key, rows, src):
    """Producer ``dots``: one column of acquisition dots per group with a black median bar."""
    for i, (label, df) in enumerate(groups):
        v = df.reference_recovery.to_numpy(dtype=float); v = v[np.isfinite(v)] * 100
        x = i + np.linspace(-.16, .16, len(v)); ax.scatter(x, v, s=8, alpha=.48, color=COL[i % 5], lw=0)
        ax.hlines(np.median(v), i - .25, i + .25, color='black', lw=1.5)
        for s_, xx, vv in zip(df['sample'], x, v):
            rows.append(dict(panel=panel, group=label.replace('\n', ' '), arm=key[i], min_length=PLOTTED_MIN_LENGTH_CD, sample=s_, x=xx, value=vv,
                             acquisitions_with_support='', label='', source=src))
        rows.append(dict(panel=panel, group=label.replace('\n', ' '), arm=key[i], min_length=PLOTTED_MIN_LENGTH_CD, sample='median', x=i, value=float(np.median(v)),
                         acquisitions_with_support='', label='', source='median of the points above'))
    ax.set_xticks(range(len(groups)), [g[0] for g in groups], fontsize=6); ax.set_ylabel(ylabel)


def panel_c(ax, t, rows):
    """Draw panel c: Sage within-acquisition reference-identity recovery per arm (5-50).

    Statistical unit: acquisition; n = 15 per arm (the REF arm is 100% by construction).
    """
    rec = t["recovery"]
    g = rec[rec.min_length.eq(PLOTTED_MIN_LENGTH_CD)]
    dots(ax, [(ARMLAB[a], g[g.arm.eq(a)]) for a in ARMS], 'Reference identities recovered (%)', 'c', ARMS, rows, TABLES["recovery"])
    ax.set_ylim(*YLIM_CD); ax.set_title('c  Sage recovery', loc='left', pad=6, fontsize=8)


def panel_d(ax, t, rows):
    """Draw panel d: within-MSFragger recovery of REF and MG identities by FL (5-50).

    Statistical unit: acquisition; n = 15 per comparator.
    """
    ms = t["msfragger"]
    g = ms[ms.min_length.eq(PLOTTED_MIN_LENGTH_CD)]
    dots(ax, [(ARMLAB[c], g[g.comparator.eq(c)]) for c in MSF_COMPARATORS], 'Reference identities recovered (%)', 'd', MSF_COMPARATORS, rows, TABLES["msfragger"])
    ax.set_ylim(*YLIM_CD); ax.set_title('d  MSFragger recovery', loc='left', pad=6, fontsize=8)


def build_figure(t):
    """Draw all four panels in the producer's order.

    Returns
    -------
    (Figure, Axes array, pandas.DataFrame of plotted values)
    """
    rows = []
    fig, axes = new_figure()
    panel_a(axes[0], t, rows)
    panel_b(axes[1], t, rows)
    panel_c(axes[2], t, rows)
    panel_d(axes[3], t, rows)
    return fig, axes, pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------
# PDF post-processing (moved from producer common.py)
# ---------------------------------------------------------------------------------------
def spans(pdf):
    """All text spans on page 1 of ``pdf`` (PyMuPDF dict output)."""
    import fitz
    pg = fitz.open(pdf)[0]
    return [dict(s, direction=l['dir']) for b in pg.get_text('dict')['blocks'] for l in b.get('lines', []) for s in l['spans']]


def slant(pdf_in, pdf_out, labels, shear=SLANT_SHEAR):
    """Re-set the given labels in Söhne Buch with a synthetic slant (label edit FAR16).

    Removes only the label text and writes it again with a shear transform. Without the
    Söhne file the labels are re-set in PyMuPDF's built-in Helvetica.

    Returns
    -------
    int
        Number of labels re-set.
    """
    import fitz
    doc = fitz.open(pdf_in); page = doc[0]
    hits = [s for s in spans(pdf_in) if s['text'] in labels]
    for s in hits:
        page.add_redact_annot(fitz.Rect(s['bbox']), fill=None, cross_out=False)
    page.apply_redactions(images=fitz.PDF_REDACT_IMAGE_NONE, graphics=fitz.PDF_REDACT_LINE_ART_NONE, text=fitz.PDF_REDACT_TEXT_REMOVE)
    if _FONT["path"] is not None:
        page.insert_font(fontname='SB', fontfile=str(_FONT["path"])); fname = 'SB'
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


def render_png(pdf, png, dpi=RENDER_DPI):
    """Rasterise page 1 of ``pdf`` to ``png`` at ``dpi``; return (width, height) in pixels."""
    import fitz
    pix = fitz.open(pdf)[0].get_pixmap(dpi=dpi, alpha=False)
    pix.set_dpi(dpi, dpi); pix.save(png)
    return pix.width, pix.height


def min_text_pt(pdf, placed_width_mm=PLACED_W_MM):
    """Smallest text size (pt) at placed width, and the scale factor page -> placed."""
    import fitz
    pg = fitz.open(pdf)[0]; k = placed_width_mm * MM / pg.rect.width
    sz = [s['size'] * k for s in spans(pdf) if s['text'].strip()]
    return round(min(sz), 3), k


def save_figure(fig, t, out, stem='SFig03_len950', dpi=RENDER_DPI):
    """Save the figure as the producer did: raw PDF, slanted final PDF, 400-dpi PNG.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    t : dict of pandas.DataFrame
    out : pathlib.Path
        Output directory.
    stem : str
        Output file stem (``SFig03_len950``: the 9-50 producer step of the v30 figure).
    dpi : int
        PNG resolution (the v28f image was 400 dpi; the len950 PNG was rendered at 300).

    Returns
    -------
    dict
        Paths of ``raw_pdf``, ``pdf`` and ``png``, the PNG size, the number of slanted
        labels and the smallest text size at placed width.
    """
    out = Path(out)
    raw = out / 'SFig03_raw.pdf'; fig.savefig(raw); plt.close(fig)
    final = out / f'{stem}.pdf'
    n_sl = slant(raw, final, set(short_names(list(t["support"].member.unique()))))
    assert n_sl == N_SLANTED_LABELS, n_sl
    mt, k = min_text_pt(final)
    assert abs(k - 1) < 1e-6 and mt >= MIN_TEXT_PT, (mt, k)
    png = out / f'{stem}.png'; px = render_png(final, png, dpi)
    return dict(raw_pdf=raw, pdf=final, png=png, png_px=px, slanted=n_sl, smallest_text_pt=mt)


def panel_crop_boxes(fig, axes, dpi=RENDER_DPI, margin=12):
    """Pixel crop box of each panel in a ``dpi`` render: axes tight bbox plus ``margin``.

    Returns
    -------
    dict
        ``{'a': (left, top, right, bottom), ...}`` in pixels, origin top-left.
    """
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    H = fig.get_figheight() * dpi; s = dpi / fig.dpi; boxes = {}
    for letter, ax in zip('abcd', axes):
        bb = ax.get_tightbbox(r)
        boxes[letter] = (max(0, int(bb.x0 * s) - margin), max(0, int(H - bb.y1 * s) - margin),
                         int(np.ceil(bb.x1 * s)) + margin, int(np.ceil(H - bb.y0 * s)) + margin)
    return boxes


def crop(png, box, out):
    """Cut ``box`` (left, top, right, bottom) out of ``png`` and save it to ``out``."""
    from PIL import Image
    Image.open(png).crop(box).save(out)
    return out


def crop_pdf(pdf, box, out, dpi=RENDER_DPI):
    """Write a one-page PDF of ``pdf`` with its crop box set to the pixel ``box`` at ``dpi``.

    The vector content is unchanged; only the visible area is limited to the panel.
    """
    import fitz
    k = 72 / dpi
    doc = fitz.open(pdf); page = doc[0]
    page.set_cropbox(fitz.Rect(box[0] * k, box[1] * k, box[2] * k, box[3] * k) & page.mediabox)   # clip to the page
    doc.save(out, garbage=4, deflate=True)
    return out


# ---------------------------------------------------------------------------------------
# Statistics and printed numbers
# ---------------------------------------------------------------------------------------
_WORDS = ["zero", "one", "two", "three", "four", "five", "six", "seven", "eight", "nine",
          "ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen"]


def as_word(n):
    """Spell an integer 0-15 as the manuscript prints it (e.g. 12 -> 'twelve')."""
    return _WORDS[int(n)]


def read_membership(sd):
    """Read the 15 per-acquisition membership tables.

    Returns
    -------
    dict
        sample -> DataFrame with columns peptide, n_members, members, arms.
    """
    return {s: pd.read_csv(Path(sd) / f"CAMPI_membership_{s}.tsv", sep="\t") for s in SAMPLES}


def stat_recount_membership(sd, t):
    """Recount panels a-c tables from the accepted-peptide membership tables.

    Port of the member recount block of the producer ``build_revision.py`` and the
    paired-recovery block of ``CAMPI/aggregate.py``. For each acquisition, arm and minimum
    length (5 or 9, maximum 50), a peptide is member-specific when it maps to exactly one
    of the eight reference proteomes. Pooled counts are the size of the union over the 15
    acquisitions; support is the number of acquisitions with at least one such peptide.
    Reference recovery is |arm & REF| / |REF| per acquisition.

    Statistical unit: acquisition (n = 15 per arm and length); peptide for pooled counts.

    Parameters
    ----------
    sd : pathlib.Path
    t : dict of pandas.DataFrame
        The published tables, compared against the recount.

    Returns
    -------
    dict
        ``support``, ``acq_counts``, ``recovery`` (recounted DataFrames) and boolean
        ``equal_support``, ``equal_acq_counts``, ``equal_recovery`` (exact equality with the
        published tables; recovery ratios compared to 1e-12), and ``lengths_in_window``
        (all peptides between MIN_LENGTHS[0] and MAX_LENGTH residues).
    """
    mem = read_membership(sd)
    names = list(t["support"].member.unique())
    sets, recs, lens_ok = {}, [], True
    for s in SAMPLES:
        m = mem[s]
        lens_ok &= bool(m.peptide.str.len().between(MIN_LENGTHS[0], MAX_LENGTH).all())
        for arm in ARMS:
            a = m[m.arms.fillna('').str.split(';').map(lambda z: arm in z)]
            for lo in MIN_LENGTHS:
                a1 = a[a.peptide.str.len().between(lo, MAX_LENGTH)]
                sets[(s, arm, lo)] = set(a1.peptide)
                for member in names:
                    v = set(a1.loc[a1.n_members.eq(SPECIFIC_N_MEMBERS) & a1.members.eq(member), 'peptide'])
                    recs.append(dict(sample=s, arm=arm, min_length=lo, member=member, specific_peptides=len(v)))
                    sets[(s, arm, lo, member)] = v
    acq = pd.DataFrame(recs)
    pooled = []
    for arm in ARMS:
        for lo in MIN_LENGTHS:
            for member in names:
                g = acq[acq.arm.eq(arm) & acq.min_length.eq(lo) & acq.member.eq(member)]
                pooled.append(dict(arm=arm, min_length=lo, member=member,
                                   pooled_specific_peptides=len(set().union(*(sets[(x, arm, lo, member)] for x in g['sample']))),
                                   acquisitions_with_support=int((g.specific_peptides > 0).sum()), total_acquisitions=len(g)))
    sup = pd.DataFrame(pooled)
    rec = []
    for s in SAMPLES:
        for arm in ARMS:
            for lo in MIN_LENGTHS:
                p, r = sets[(s, arm, lo)], sets[(s, 'REF', lo)]
                rec.append(dict(sample=s, arm=arm, min_length=lo, peptides=len(p), reference=len(r), overlap=len(p & r),
                                gained=len(p - r), lost=len(r - p), count_ratio=len(p) / len(r), reference_recovery=len(p & r) / len(r)))
    rec = pd.DataFrame(rec)
    key = ['sample', 'arm', 'min_length']
    pub = t["recovery"].sort_values(key).reset_index(drop=True)
    new = rec.set_index(key).loc[pd.MultiIndex.from_frame(pub[key])].reset_index()
    ints = ['peptides', 'reference', 'overlap', 'gained', 'lost']
    eq_rec = (len(pub) == len(rec) and (new[ints].to_numpy() == pub[ints].to_numpy()).all()
              and np.allclose(new[['count_ratio', 'reference_recovery']], pub[['count_ratio', 'reference_recovery']], rtol=0, atol=1e-12))
    eq_sup = sup.reset_index(drop=True).equals(t["support"].reset_index(drop=True))
    eq_acq = acq.reset_index(drop=True).equals(t["acq_counts"].reset_index(drop=True))
    return dict(support=sup, acq_counts=acq, recovery=rec, equal_support=bool(eq_sup),
                equal_acq_counts=bool(eq_acq), equal_recovery=bool(eq_rec), lengths_in_window=lens_ok)


def stat_line83_member_support(t):
    """Numbers of Results v28d line 83 from the member-support table.

    Statistical unit: peptide (distinct peptides pooled over n = 15 acquisitions).

    Returns
    -------
    dict
        ``n_members`` (members with >= 1 FL member-specific peptide at 5-50),
        ``n_members_total``, ``lp_fl_5`` and ``lp_fl_9`` (L. plantarum pooled FL counts),
        ``mg_members_9`` (members with >= 1 MG member-specific peptide at 9-50),
        ``mg_zero_9`` (members with none), ``least3`` (the three members with the fewest FL
        member-specific peptides at 5-50), ``fold_ref_over_fl`` (REF / FL pooled count for
        those three), ``fold_min`` and ``fold_max``, and ``flr_exceeds_ref`` (FL + Reference
        exceeds REF for each of the three).
    """
    p = t["support"]
    def get(arm, lo):
        return p[p.arm.eq(arm) & p.min_length.eq(lo)].set_index('member').pooled_specific_peptides
    fl5, fl9, mg9, ref5, flr5 = get('FL', 5), get('FL', 9), get('MG', 9), get('REF', 5), get('FL_R', 5)
    lp = [m for m in p.member.unique() if m.endswith('plantarum')][0]
    least3 = list(fl5.sort_values().index[:3])
    fold = (ref5[least3] / fl5[least3]).to_dict()
    return dict(n_members=int((fl5 > 0).sum()), n_members_total=int(p.member.nunique()),
                lp_fl_5=int(fl5[lp]), lp_fl_9=int(fl9[lp]),
                mg_members_9=int((mg9 > 0).sum()), mg_zero_9=sorted(mg9[mg9 == 0].index),
                least3=least3, fold_ref_over_fl=fold, fold_min=min(fold.values()), fold_max=max(fold.values()),
                flr_exceeds_ref={m: bool(flr5[m] > ref5[m]) for m in least3})


def stat_line89_msfragger(t, lo=PLOTTED_MIN_LENGTH_CD):
    """Median paired FL/comparator peptide-count ratio and identity recovery in MSFragger.

    Descriptive medians over acquisitions of the within-MSFragger comparison (Results
    v28d line 89). count_ratio = |FL| / |comparator|; reference_recovery =
    |FL & comparator| / |comparator|; peptides accepted at Percolator peptide q <= 0.01.
    Medians use numpy.median (mean of the two middle values for even n; here n is odd).

    Statistical unit: acquisition; n = 15 per comparator.

    Returns
    -------
    dict
        ``n`` per comparator, ``ratio_REF``, ``ratio_MG``, ``recovery_REF``,
        ``recovery_MG`` (fractions) and ``recovery_pooled`` (median over all 30
        acquisition-comparator pairs), plus ``fl_more_REF`` / ``fl_more_MG`` (acquisitions
        in which FL gives more peptides).
    """
    m = t["msfragger"]; g = m[m.min_length.eq(lo)]
    r = {}
    for c in MSF_COMPARATORS:
        z = g[g.comparator.eq(c)]
        r['n_' + c] = int(len(z)); r['ratio_' + c] = float(np.median(z.count_ratio))
        r['recovery_' + c] = float(np.median(z.reference_recovery)); r['fl_more_' + c] = int((z.FL > z.reference).sum())
    r['recovery_pooled'] = float(np.median(g.reference_recovery))
    return r


def stat_panel_medians(plotted):
    """Median bars of panels c and d as drawn (percent), from the plotted-values table.

    Statistical unit: acquisition; n = 15 per group.

    Returns
    -------
    pandas.DataFrame
        panel, arm, n, median_percent.
    """
    pv = plotted[plotted.panel.isin(['c', 'd'])]
    pts = pv[pv['sample'].ne('median')]
    med = pv[pv['sample'].eq('median')].set_index(['panel', 'arm']).value
    out = pts.groupby(['panel', 'arm'], sort=False).value.agg(['size', 'median']).reset_index()
    out.columns = ['panel', 'arm', 'n', 'median_percent']
    assert np.allclose(out.median_percent.to_numpy(), med.loc[list(zip(out.panel, out.arm))].to_numpy())
    return out


def n_acquisitions(t):
    """Number of distinct CAMPI acquisitions in each table (expected 15 everywhere)."""
    return {k: int(v['sample'].nunique()) for k, v in t.items() if 'sample' in v}


def stat_msfragger_receipts(t):
    """Check panel d's table against the 45 search receipts and report engine versions.

    The within-MSFragger table holds, per acquisition and length, |FL| and |comparator|.
    Each must equal the accepted count written by the search itself (COMPLETE.json of the
    upstream run, exported to CAMPI_MSFragger_search_receipts.tsv). The FL arm is named
    ``FL_UNIPROT`` in the receipts.

    Statistical unit: acquisition; n = 15 per arm (45 searches).

    Returns
    -------
    dict
        ``n_searches``, ``n_acquisitions``, ``all_pass``, ``engine_strings`` (distinct
        version strings), ``msfragger_version``, ``percolator_version`` and
        ``counts_match`` (every FL and comparator count equals its receipt).
    """
    import re
    r = t["receipts"]; m = t["msfragger"]
    acc = {(x.sample, x.arm, 5): x.accepted_5_50 for x in r.itertuples()}
    acc.update({(x.sample, x.arm, 9): x.accepted_9_50 for x in r.itertuples()})
    ok = all(acc[(x.sample, 'FL_UNIPROT', x.min_length)] == x.FL and acc[(x.sample, x.comparator, x.min_length)] == x.reference
             for x in m.itertuples())
    eng = sorted(r.engine_and_confidence.unique())
    mf = sorted({re.search(r'MSFragger (\S+?);', e).group(1) for e in eng})
    pc = sorted({re.search(r'Percolator (\S+) ', e).group(1) for e in eng})
    return dict(n_searches=int(len(r)), n_acquisitions=int(r['sample'].nunique()), all_pass=bool(r.status.eq('PASS').all()),
                engine_strings=eng, msfragger_version=';'.join(mf), percolator_version=';'.join(pc), counts_match=bool(ok))


def recomputed_numbers(t, s83, s89, rcpt):
    """Every number of NUMBERS.tsv, in row order, as computed here.

    Parameters
    ----------
    t : dict of pandas.DataFrame
    s83, s89, rcpt : dict
        Outputs of ``stat_line83_member_support``, ``stat_line89_msfragger`` and
        ``stat_msfragger_receipts``.

    Returns
    -------
    list of (str, object)
        (row label, value). Word numbers are spelled as printed; percentages are in %.
    """
    lo5, lo9 = MIN_LENGTHS
    n_acq = sorted(set(t["support"].total_acquisitions))
    assert n_acq == [N_ACQUISITIONS] and rcpt["n_acquisitions"] == N_ACQUISITIONS
    return [
        ("83 members with FL member-specific peptides", as_word(s83["n_members"])),
        ("83 L. plantarum FL pooled 5-50", as_word(s83["lp_fl_5"])),
        ("83 lower length bound", lo5), ("83 upper length bound", MAX_LENGTH),
        ("83 L. plantarum FL pooled 9-50", as_word(s83["lp_fl_9"])),
        ("83 nested lower length bound", lo9), ("83 upper length bound", MAX_LENGTH),
        ("83 MG members with member-specific peptides at 9-50", as_word(s83["mg_members_9"])),
        ("83 nested lower length bound", lo9), ("83 upper length bound", MAX_LENGTH),
        ("83 REF/FL fold, least of three members (5-50)", s83["fold_min"]),
        ("83 REF/FL fold, most of three members (5-50)", s83["fold_max"]),
        ("89 median FL/REF count ratio", s89["ratio_REF"]),
        ("89 median FL/MG count ratio", s89["ratio_MG"]),
        ("89 median comparator identity recovery (%)", 100 * s89["recovery_pooled"]),
        ("253 CAMPI acquisitions searched", rcpt["n_acquisitions"]),
        ("253 MSFragger version", rcpt["msfragger_version"]),
        ("253 Percolator version", rcpt["percolator_version"]),
        ("417 acquisitions denominator", n_acq[0]),
        ("417 reference proteomes", as_word(s83["n_members_total"])),
        ("417 nested lower length bound", lo9), ("417 upper length bound", MAX_LENGTH),
        ("417 lower length bound", lo5), ("417 upper length bound", MAX_LENGTH),
        ("685 members", as_word(s83["n_members_total"])),
        ("685 lower length bound", lo5), ("685 upper length bound", MAX_LENGTH),
        ("685 nested lower length bound", lo9), ("685 upper length bound", MAX_LENGTH),
    ]


# =====================================================================================
# v30 house-style rebuild: shared helpers
# Moved unchanged from the v30 rebuild helper module (restyle_common.py, used for
# Supplementary Figures 1-3; see PROVENANCE.md). Only the lab paths and environment
# setup at the top of that module were dropped (sha and render_png are this module's own, same code); fl_style is imported from the
# repository's style package instead. Functions:
#   fingerprint()   plotted arrays + axis limits/scales of a list of axes (values-unchanged proof)
#   compare()       exact comparison of two fingerprints
#   texts()         every visible string (texts, tick labels, legend entries) of a figure
#   capture_saved() run a producer and capture the Figure objects it saves (nothing written)
#   align_letters() put the bold panel letters of one figure column on one x position
#   draw_rel_label() draw a label containing <= or >= in Söhne (the face has no such glyphs)
# =====================================================================================
import hashlib  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.collections import Collection, LineCollection, PathCollection  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402


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




# =====================================================================================
# v30: house-style rebuild (moved unchanged from build_sfig03_style.py)
# =====================================================================================
V30_PLACED_W_MM, V30_PLACED_H_MM = 160.0, 140.0   # v30 build: placed size (v29 docx extent 160.0 x 140.0 mm)
V30_PNG_DPI = 300                                 # v30 build: PNG rasterised at 300 dpi from the PDF
V30_TITLES = {'a  FastaLake, 5–50 residues', 'b  FastaLake, 9–50 residues', 'c  Sage recovery', 'd  MSFragger recovery',
              'Pooled member-specific\npeptide count'}   # removed in v30; their words are now in the x labels


def build_figure_v30(t):
    """Draw panels a-d in the house style (v30) from the same tables as ``build_figure``.

    Arm colours as in main Figure 2 (FL blue, FL + Reference green, FL + Reference + MG
    reddish purple, Reference grey, MG orange); both a and b and both d columns are
    FastaLake searches, so they are FL blue. Subpanel titles move into the x-axis labels;
    c is widened and d narrowed (axes positions only).

    Parameters
    ----------
    t : dict of pandas.DataFrame
        Output of ``load``.

    Returns
    -------
    tuple
        ``(fig, axes, plotted)``; ``plotted`` has one row per plotted value (panel, group,
        value, label) in drawing order.
    """
    from style import fl_style as fl
    matplotlib.rcdefaults()
    fl.init()
    ARM_COL = {'FL': fl.FL, 'FL_R': fl.GREEN, 'FL_R_M': fl.PURPLE, 'REF': fl.REF, 'MG': fl.METAG}   # = main Fig. 2
    fig, axs = plt.subplots(2, 2, figsize=fl.mm(V30_PLACED_W_MM, V30_PLACED_H_MM)); ax = axs.ravel()
    fig.subplots_adjust(left=.165, right=.975, top=.965, bottom=.115, wspace=.62, hspace=.50)
    # bottom row: c is widened (its arm labels collided at the producer's width); d narrows.
    _y0, _h = ax[2].get_position().y0, ax[2].get_position().height
    ax[2].set_position([0.085, _y0, 0.435, _h]); ax[3].set_position([0.64, _y0, 0.335, _h])
    rows = []
    p = t['support']
    names = list(p.member.unique()); short = short_names(names)
    for i, lo in enumerate(MIN_LENGTHS):
        a = ax[i]
        sub = p[p.arm.eq('FL') & p.min_length.eq(lo)].set_index('member').loc[names]
        a.barh(np.arange(8), sub.pooled_specific_peptides, color=fl.FL)
        a.set_yticks(np.arange(8), short); a.set_xscale('log'); a.invert_yaxis()
        a.set_xlabel(f'Pooled member-specific peptide count,\nFastaLake, {lo}–50 residues')
        for k, x in enumerate(sub.itertuples()):
            lab = f'{x.pooled_specific_peptides:,}; {x.acquisitions_with_support}/15'
            a.text(max(x.pooled_specific_peptides, 1) * 1.25, k, lab, va='center', fontsize=6.5, color=fl.INK)
            rows.append(dict(panel='ab'[i], group=x.Index, value=x.pooled_specific_peptides, label=lab))
        a.set_xlim(*XLIM_AB)

    def _dots(a, groups, cols, panel, ylabel='Reference identities recovered (%)'):
        for i, (label, df) in enumerate(groups):
            v = df.reference_recovery.to_numpy(dtype=float); v = v[np.isfinite(v)] * 100
            x = i + np.linspace(-.16, .16, len(v))
            a.scatter(x, v, s=7, alpha=.6, color=cols[i], lw=0)
            a.hlines(np.median(v), i - .25, i + .25, color=fl.INK, lw=1.5)
            rows.extend(dict(panel=panel, group=label.replace('\n', ' '), sample=s_, value=vv) for s_, vv in zip(df['sample'], v))
            rows.append(dict(panel=panel, group=label.replace('\n', ' '), sample='median', value=float(np.median(v))))
        a.set_xticks(range(len(groups)), [g[0] for g in groups], fontsize=6.5)
        a.set_ylabel(ylabel)

    rec = t['recovery']; g = rec[rec.min_length.eq(PLOTTED_MIN_LENGTH_CD)]
    _dots(ax[2], [(ARMLAB[k], g[g.arm.eq(k)]) for k in ARMS], [ARM_COL[k] for k in ARMS], 'c')
    ax[2].set_ylim(*YLIM_CD); ax[2].set_xlabel('Sage search arm')
    ms = t['msfragger']; g = ms[ms.min_length.eq(PLOTTED_MIN_LENGTH_CD)]
    # d: every dot is a FastaLake search (its recovery of the comparator's identities) -> FL blue
    _dots(ax[3], [(ARMLAB[c], g[g.comparator.eq(c)]) for c in MSF_COMPARATORS], [fl.FL, fl.FL], 'd',
          ylabel='Comparator identities recovered (%)')   # v35: d recovers the comparator's (Reference or MG) identities
    ax[3].set_ylim(*YLIM_CD); ax[3].set_xlabel('MSFragger comparator')
    for a in ax:
        fl.style(a)
    for a in ax[:2]:
        a.tick_params(axis='x', labelsize=fl.LOG_TICK_PT)
    ax[2].tick_params(axis='x', labelsize=6.5); ax[3].tick_params(axis='x', labelsize=6.5)
    for a, l in zip(ax, 'abcd'):
        fl.panel(a, l)
    fl.normalise_type(fig, 'pre-align', verbose=False)
    align_letters(fig, [[ax[0], ax[2]], [ax[1], ax[3]]])
    return fig, list(ax), pd.DataFrame(rows)


def build_v30(t, out):
    """Build the v30 figure and prove its plotted values equal the producer's.

    Steps of ``build_sfig03_style.py``: (1) run ``init`` and ``build_figure`` (the len950
    producer, c,d at 9-50 residues) and fingerprint its axes; (2) ``build_figure_v30``;
    (3) assert identical arrays, limits, scales, plotted-values table and a,b value labels,
    and that only the four titles and the old two-line x label are lost as strings;
    (4) save the PDF, re-set the 16 taxon labels with a slant (``slant``) and rasterise at
    300 dpi to ``SFig03_v30.png``.

    Parameters
    ----------
    t : dict of pandas.DataFrame
        Output of ``load``.
    out : pathlib.Path
        Output directory.

    Returns
    -------
    dict
        Comparison report, ``slanted`` label count, ``pdf``, ``png``, ``png_px`` and the
        v30 plotted-values table (``plotted``).
    """
    from style import fl_style as fl
    out = Path(out)
    have_font = init()
    old_fig, old_axes, old_plotted = build_figure(t)
    OLD = fingerprint(list(old_axes))
    old_text = texts(old_fig)
    plt.close(old_fig)
    fig, axes, newp = build_figure_v30(t)
    NEW = fingerprint(axes)
    new_text = texts(fig)
    ok, rep = compare(OLD, NEW, ['a', 'b', 'c', 'd'])
    lost = sorted(set(old_text) - set(new_text)); added = sorted(set(new_text) - set(old_text))
    assert ok, rep
    assert set(lost) <= V30_TITLES, lost
    oldv = old_plotted[['panel', 'value']].reset_index(drop=True); newv = newp[['panel', 'value']].reset_index(drop=True)
    plotted_equal = bool(oldv.panel.equals(newv.panel) and np.array_equal(oldv.value.to_numpy(float), newv.value.to_numpy(float)))
    assert plotted_equal
    labels_equal = list(old_plotted.loc[old_plotted.panel.isin(['a', 'b']), 'label']) == list(newp.loc[newp.panel.isin(['a', 'b']), 'label'])
    assert labels_equal
    raw, _ = fl.save(fig, out / '_raw_SFig03_v30', formats=('pdf',))
    plt.close(fig)
    pdf = out / 'SFig03_v30.pdf'
    short = short_names(list(t['support'].member.unique()))
    n_sl = slant(raw, pdf, set(short))
    assert n_sl == N_SLANTED_LABELS, n_sl
    Path(raw).unlink()
    png = out / 'SFig03_v30.png'
    px = render_png(pdf, png, V30_PNG_DPI)
    return dict(values_identical=ok, panels=rep, plotted_table_equal_to_producer=plotted_equal,
                value_labels_equal=labels_equal, label_strings_lost=lost, label_strings_added=added,
                slanted=n_sl, pdf=pdf, png=png, png_px=px, plotted=newp,
                producer_font_available=have_font)
