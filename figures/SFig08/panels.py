"""Supplementary Figure 8: clinical database comparators and evidence ablation.

Panels
------
a  accepted peptides in four PREDICT database arms, 296 acquisitions
   (``source_data/SFig08/PREDICT_curation_counts.tsv``)
b  LFQ intensity removed when FastaLake-only peptide evidence is dropped
c  KEGG orthologs (KOs) losing all support
d  Bray-Curtis distance between normalized all-evidence and shared-only KO profiles
   (b-d: ``source_data/SFig08/clinical_fixed_dictionary_ablation.tsv``, 296 PREDICT and
   295 CapScan acquisitions)

The original producer drew a-d as one matplotlib figure (block "Clinical comparator scope"
of ``build_revision.py``) and a second script removed the figure title from the PDF and
rasterised it at 300 dpi (``prepare_supplement_figures.py``). Both steps are moved here
unchanged: ``figure`` reproduces the first, ``redact_title`` the second. Original
locations are in ``PROVENANCE.md``.

The figure shows no hypothesis test. Its statistics are medians of per-acquisition
values and counts of acquisitions; see ``stat_a_paired_gain`` and ``stat_bcd_medians``.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

# ---------------------------------------------------------------------------------------
# Choices. Every value is the one the original producer used.
# ---------------------------------------------------------------------------------------

#: rcParams set by the producer on top of matplotlib defaults (build_revision.py, line 13).
PRODUCER_RC = {
    "font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "pdf.fonttype": 42, "svg.fonttype": "none",
    "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 180,
}
#: Producer palette (build_revision.py, line 14); colours cycle by group position.
COL = ["#2878a3", "#df8f44", "#619b68", "#b35872", "#8277ac"]
#: Figure size in inches (producer: plt.subplots(2, 2, figsize=(7.1, 6.5))).
FIGSIZE = (7.1, 6.5)
#: Figure title drawn by the producer and later removed by redact_title.
FIGURE_TITLE = "Comparator curation and the contribution of additional evidence"
#: Rectangle (PDF points) that prepare_supplement_figures.py redacted; recorded by the
#: producer in SUPPLEMENT_FIGURE_SOURCES.json entry 6 ("top_text" bbox).
TITLE_BBOX = (54.334373474121094, 7.500360488891602, 456.86846923828125, 20.31536102294922)
#: Rasterisation resolution of the composite (prepare_supplement_figures.py).
COMPOSITE_DPI = 300
#: In-figure panel titles, as drawn by the producer.
PANEL_TITLES = {"a": "PREDICT comparator curation", "b": "Removed peptide evidence",
                "c": "Functional support", "d": "Profile change"}
#: Ablation columns plotted in b-d and their y labels, as drawn by the producer.
ABLATION_COLUMNS = [("intensity_removed_percent", "LFQ intensity removed (%)"),
                    ("KO_support_lost_percent", "KO support lost (%)"),
                    ("BrayCurtis", "Normalized KO Bray–Curtis distance")]
#: Peptide acceptance and length window behind panel a counts: rank-1 target PSMs with
#: peptide q <= 0.01, 5-50 residues, modifications stripped, I and L merged
#: (Methods, v28d line 229; upstream collector, see UPSTREAM.json).
Q_CUTOFF = 0.01
LENGTH_WINDOW = (5, 50)
#: Arm names in the source table and their legend names (v28d line 441).
ARMS = {"FL_cohortMG": "full FastaLake including cohort metagenomes",
        "MG": "specimen MG",
        "FL_public": "public-sequence FastaLake",
        "MG_own_curated": "curated specimen MG"}
#: Expected number of rows per arm and per cohort (legend, v28d line 441).
N_PREDICT = 296
N_CAPSCAN = 295
#: Spelled-out counts, to compare with numbers the legend prints as words.
NUMBER_WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five"}


def read_tables(sd):
    """Read the two source tables.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig08``.

    Returns
    -------
    tuple of pandas.DataFrame
        ``(counts, ablation)``: long table (sample, arm, peptides), one row per
        acquisition and arm; ablation table, one row per acquisition.
    """
    counts = pd.read_csv(Path(sd) / "PREDICT_curation_counts.tsv", sep="\t")
    abl = pd.read_csv(Path(sd) / "clinical_fixed_dictionary_ablation.tsv", sep="\t")
    return counts, abl


def check_tables(counts, abl):
    """Run the producer's numerical checks and the roster checks; raise on failure.

    Parameters
    ----------
    counts, abl : pandas.DataFrame
        Tables from ``read_tables``.

    Returns
    -------
    dict
        Check name to True.
    """
    out = {}
    # Producer checks (build_revision.py, "Ablation intensity fractions" and "Ablation KO fraction").
    out["ablation_intensity_fraction"] = bool(np.allclose(
        abl.fraction_intensity_removed, abl.FL_only_LFQ_intensity / abl.all_LFQ_intensity))
    out["ablation_KO_fraction"] = bool(np.allclose(
        abl.KO_support_lost_percent, 100 * (abl.all_KOs - abl.shared_KOs) / abl.all_KOs))
    # Roster: four arms on the same 296 acquisitions; 296 PREDICT and 295 CapScan rows.
    per_arm = counts.groupby("arm")["sample"].agg(["size", "nunique"])
    out["four_arms"] = sorted(per_arm.index) == sorted(ARMS)
    out["296_per_arm"] = bool((per_arm["size"] == N_PREDICT).all() and (per_arm["nunique"] == N_PREDICT).all())
    out["same_acquisitions_all_arms"] = len({frozenset(g["sample"]) for _, g in counts.groupby("arm")}) == 1
    n = abl.groupby("cohort")["sample"].nunique().to_dict()
    out["ablation_roster"] = n == {"PREDICT": N_PREDICT, "CapScan": N_CAPSCAN} and abl["sample"].is_unique
    out["PREDICT_ablation_equals_curation_roster"] = (
        set(abl.loc[abl.cohort == "PREDICT", "sample"]) == set(counts["sample"]))
    bad = [k for k, v in out.items() if not v]
    if bad:
        raise AssertionError(f"failed checks: {bad}")
    return out


# ---------------------------------------------------------------------------------------
# Plotting, moved from the producer unchanged.
# ---------------------------------------------------------------------------------------

def dots(ax, groups, ylabel, percent=False):
    """Draw one evenly spread dot column per group with a black median bar (producer code).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    groups : list of (str, array-like)
        Tick label and values, in drawing order.
    ylabel : str
    percent : bool
        Multiply values by 100.
    """
    for i, (label, values) in enumerate(groups):
        v = np.asarray(values, dtype=float); v = v[np.isfinite(v)] * (100 if percent else 1)
        x = i + np.linspace(-.16, .16, len(v)); ax.scatter(x, v, s=8, alpha=.48, color=COL[i % 5], rasterized=len(v) > 100)
        if len(v): ax.hlines(np.median(v), i - .25, i + .25, color='black', lw=1.5)
    ax.set_xticks(range(len(groups)), [g[0] for g in groups]); ax.set_ylabel(ylabel)


def title(ax, l, t):
    """Left-aligned panel letter and in-figure title (producer code)."""
    ax.set_title(l + '  ' + t, loc='left', pad=9)


def panel_a(ax, counts):
    """Panel a: accepted peptides (thousands) per acquisition in each of four PREDICT arms.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    counts : pandas.DataFrame
        Long table (sample, arm, peptides); 296 acquisitions per arm.

    Returns
    -------
    None
        Arms are drawn in the producer's order (pandas groupby, alphabetical).
    """
    dots(ax, [(k.replace('_', '\n'), g.peptides / 1000) for k, g in counts.groupby('arm')], 'Accepted peptides (thousands)'); title(ax, 'a', PANEL_TITLES['a'])


def _panel_ablation(ax, abl, i):
    col, label = ABLATION_COLUMNS[i]
    dots(ax, [(co, g[col]) for co, g in abl.groupby('cohort')], label); title(ax, 'bcd'[i], PANEL_TITLES['bcd'[i]])


def panel_b(ax, abl):
    """Panel b: LFQ intensity removed (%) per acquisition, by cohort (CapScan 295, PREDICT 296).

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    abl : pandas.DataFrame
        Ablation table, one row per acquisition.
    """
    _panel_ablation(ax, abl, 0)


def panel_c(ax, abl):
    """Panel c: KOs losing all support (%) per acquisition, by cohort. See ``panel_b``."""
    _panel_ablation(ax, abl, 1)


def panel_d(ax, abl):
    """Panel d: normalized KO-profile Bray-Curtis distance per acquisition. See ``panel_b``."""
    _panel_ablation(ax, abl, 2)


def figure(counts, abl, out_stem):
    """Draw panels a-d exactly as the producer did and save PDF, SVG and PNG.

    Parameters
    ----------
    counts, abl : pandas.DataFrame
        Tables from ``read_tables``.
    out_stem : pathlib.Path
        Output path without suffix.

    Returns
    -------
    dict
        Suffix to written path; and the axes bounding boxes (inches) for cropping panels.
    """
    rc = {k: v for k, v in mpl.rcParamsDefault.items() if k not in ("backend", "backend_fallback")}
    rc.update(PRODUCER_RC)
    out_stem = Path(out_stem)
    with mpl.rc_context(rc):
        fig, axes = plt.subplots(2, 2, figsize=FIGSIZE); axes = axes.ravel()
        panel_a(axes[0], counts)
        for f, ax in zip([panel_b, panel_c, panel_d], axes[1:]):
            f(ax, abl)
        # Producer save(): bold suptitle, tight layout, then pdf, svg, png in that order.
        fig.suptitle(FIGURE_TITLE, fontsize=11, fontweight='bold')
        fig.tight_layout(rect=(0, 0, 1, .96), pad=1.6)
        written = {}
        for ext in ['pdf', 'svg', 'png']:
            p = out_stem.with_suffix('.' + ext); fig.savefig(p); written[ext] = p
        fig.canvas.draw()
        r = fig.canvas.get_renderer()
        boxes = {l: ax.get_tightbbox(r).transformed(fig.dpi_scale_trans.inverted()).bounds
                 for l, ax in zip('abcd', axes)}
        plt.close(fig)
    written["boxes_in"] = boxes
    return written


def redact_title(pdf_in, pdf_out, dpi=COMPOSITE_DPI):
    """Remove the figure title and rasterise, as prepare_supplement_figures.py did.

    The title rectangle is filled white with a PyMuPDF redaction (images and vector
    graphics untouched); the step asserts that every other number in the page text is
    unchanged, saves the PDF and writes a PNG at ``dpi``.

    Parameters
    ----------
    pdf_in, pdf_out : pathlib.Path
    dpi : int

    Returns
    -------
    pathlib.Path
        The PNG written next to ``pdf_out``.
    """
    import collections
    import re
    import fitz
    doc = fitz.open(pdf_in); p = doc[0]; before = p.get_text()
    p.add_redact_annot(fitz.Rect(TITLE_BBOX), fill=(1, 1, 1)); p.apply_redactions(images=0, graphics=0)
    assert FIGURE_TITLE not in p.get_text()
    before_nums = collections.Counter(re.findall(r'\d[\d,.]*', before.replace(FIGURE_TITLE, '')))
    after_nums = collections.Counter(re.findall(r'\d[\d,.]*', p.get_text()))
    assert before_nums == after_nums, 'numeric label changed'
    doc.save(pdf_out, garbage=4, deflate=True)
    png = Path(pdf_out).with_suffix('.png')
    fitz.open(pdf_out)[0].get_pixmap(dpi=dpi, alpha=False).save(png)
    return png


def title_bbox(pdf):
    """Return the bbox (PDF points) of the figure-title line in ``pdf``, for checking TITLE_BBOX."""
    import fitz
    with fitz.open(pdf) as doc:
        for b in doc[0].get_text("dict")["blocks"]:
            for line in b.get("lines", []):
                text = "".join(s["text"] for s in line["spans"])
                if text == FIGURE_TITLE:
                    return tuple(line["bbox"])
    return None


def crop_panel(png, box_in, dpi, out):
    """Crop one panel from a rendered composite for display.

    Parameters
    ----------
    png : pathlib.Path
        Raster of the full figure.
    box_in : tuple
        (x0, y0, width, height) in inches from the figure's lower-left corner.
    dpi : int
        Resolution of ``png``.
    out : pathlib.Path

    Returns
    -------
    pathlib.Path
    """
    from PIL import Image
    im = Image.open(png)
    x0, y0, w, h = box_in
    H = im.size[1]
    left, right = int(np.floor(x0 * dpi)), int(np.ceil((x0 + w) * dpi))
    top, bottom = int(np.floor(H - (y0 + h) * dpi)), int(np.ceil(H - y0 * dpi))
    im.crop((max(left, 0), max(top, 0), min(right, im.size[0]), min(bottom, H))).save(out)
    return Path(out)


# ---------------------------------------------------------------------------------------
# Statistics
# ---------------------------------------------------------------------------------------

def stat_a_paired_gain(counts, test_arm, reference_arm):
    """Median paired peptide gain of one arm over another, and sign counts.

    Statistical unit: acquisition (n = 296 PREDICT acquisitions; repeated acquisitions
    of one patient are not combined). Gain per acquisition is
    100 * (test / reference - 1); the reported value is the median over acquisitions.

    Parameters
    ----------
    counts : pandas.DataFrame
        Long table (sample, arm, peptides).
    test_arm, reference_arm : str
        Arm names in ``counts``.

    Returns
    -------
    dict
        n, median_gain_percent, n_higher, n_lower, n_equal, median_test, median_reference.
    """
    w = counts.pivot(index='sample', columns='arm', values='peptides')
    t, r = w[test_arm].astype(float), w[reference_arm].astype(float)
    assert t.notna().all() and r.notna().all() and (r > 0).all()
    g = 100 * (t / r - 1)
    return dict(n=int(len(w)), median_gain_percent=float(np.median(g)), n_higher=int((t > r).sum()),
                n_lower=int((t < r).sum()), n_equal=int((t == r).sum()),
                median_test=float(np.median(t)), median_reference=float(np.median(r)))


def stat_a_n_arms(counts):
    """Number of database arms and acquisitions per arm in panel a (unit: acquisition).

    Returns
    -------
    dict
        n_arms, n_arms_word (the count as the legend spells it) and n_acquisitions (296).
    """
    n = int(counts['arm'].nunique())
    return dict(n_arms=n, n_arms_word=NUMBER_WORDS[n], n_acquisitions=int(counts['sample'].nunique()))


def stat_bcd_medians(abl):
    """Per-cohort medians of the three ablation endpoints (unit: acquisition).

    n = 296 PREDICT and 295 CapScan acquisitions. Medians are those drawn as bars in b-d;
    the manuscript prints none of them.

    Returns
    -------
    pandas.DataFrame
        Index cohort; columns n and the median of each endpoint.
    """
    cols = [c for c, _ in ABLATION_COLUMNS]
    m = abl.groupby('cohort')[cols].median()
    m.insert(0, 'n', abl.groupby('cohort')['sample'].nunique())
    return m


# =======================================================================================
# v30: panel a at 9-50 residues and the v30 house-style restyle
# (moved from len950 build.py and build_sfig08_style.py; logic unchanged)
# =======================================================================================
import json  # noqa: E402

from matplotlib.collections import LineCollection, PathCollection  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from style import fl_style as fl  # noqa: E402

LENGTH_WINDOW_V30 = (9, 50)         # v30 legend a: "Accepted peptides of 9-50 residues"
V30_W_MM, V30_H_MM = 160.0, 134.0   # v30 placed size (v29 docx extent 160 x 146.5 mm minus the empty title band)
V30_TICK_PT, V30_LAB_PT = 6.5, 7.0  # type sizes (pt) at placed size
V30_DOT_S, V30_DOT_A = 9, 0.6       # dot size (pt^2) and alpha
ARM_COLOUR = {"FL_cohortMG": fl.FL, "FL_public": fl.FL, "MG": fl.METAG, "MG_own_curated": fl.METAG}  # colours by meaning (main Fig. 4)
#: v35: drawn arm names, short forms of the legend names in ``ARMS`` (order unchanged: groupby, alphabetical).
ARM_TICK = {"FL_cohortMG": "FastaLake +\ncohort MG", "FL_public": "Public\nFastaLake",
            "MG": "Specimen\nMG", "MG_own_curated": "Curated\nspecimen MG"}
#: producer tick label -> v35 drawn label, for the before/after tick check
TICK_RENAME = {k.replace("_", "\n"): v for k, v in ARM_TICK.items()}
COHORT_COLOUR = {"CapScan": fl.FL, "PREDICT": fl.FL}   # b-d: FastaLake-only evidence; cohort is the x axis


def read_tables_v30(sd):
    """Panel a counts at 9-50 residues in the shipped row order, and the ablation table.

    As ``len950/SFig08/build/build.py``: the 9-50 recount is merged onto the (sample, arm)
    order of ``PREDICT_curation_counts.tsv``; its 5-50 column must equal that table.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig08``.

    Returns
    -------
    tuple of pandas.DataFrame
        ``(counts_9_50, ablation)``; counts long table (sample, arm, peptides).
    """
    rc5 = pd.read_csv(Path(sd) / "PREDICT_curation_counts_by_length.tsv", sep="\t")
    shipped, abl = read_tables(sd)
    lo = LENGTH_WINDOW_V30[0]
    c = rc5[["sample", "arm", f"peptides_{lo}_50"]].rename(columns={f"peptides_{lo}_50": "peptides"})
    c9 = shipped[["sample", "arm"]].merge(c, on=["sample", "arm"], how="left")
    c5 = shipped[["sample", "arm"]].merge(rc5[["sample", "arm", "peptides_5_50"]].rename(columns={"peptides_5_50": "peptides"}),
                                          on=["sample", "arm"], how="left")
    assert (c5.peptides.values == shipped.peptides.values).all()
    check_tables(c9, abl)
    return c9, abl


def v30_collect(axes):
    """Every plotted quantity of each axes, colour-free, for the before/after check."""
    axes[0].figure.canvas.draw()
    out = []
    for ax in axes:
        out.append(dict(
            xlim=list(ax.get_xlim()), ylim=list(ax.get_ylim()),
            xticks=[float(v) for v in ax.get_xticks()], yticks=[float(v) for v in ax.get_yticks()],
            xticklabels=[t.get_text() for t in ax.get_xticklabels()],
            scatter=[np.asarray(c.get_offsets(), float) for c in ax.collections if isinstance(c, PathCollection)],
            segments=[np.asarray(s, float) for c in ax.collections if isinstance(c, LineCollection) for s in c.get_segments()],
            lines=[np.asarray(l.get_xydata(), float) for l in ax.lines],
            bars=[(p.get_x(), p.get_y(), p.get_width(), p.get_height()) for p in ax.patches if isinstance(p, Rectangle)],
            texts=[(tuple(map(float, t.get_position())), t.get_text()) for t in ax.texts],
            ylabel=ax.get_ylabel(), xlabel=ax.get_xlabel(), title=ax.get_title(loc="left")))
    return out


def v30_compare(before, after):
    """Exact comparison of the producer's and the restyled plotted arrays.

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
        ok("median bars identical", len(b["segments"]) == len(a["segments"]) and all(np.array_equal(x, y) for x, y in zip(b["segments"], a["segments"])), f"{len(a['segments'])} bars")
        ok("lines identical", len(b["lines"]) == len(a["lines"]) and all(np.array_equal(x, y) for x, y in zip(b["lines"], a["lines"])))
        ok("bars identical", b["bars"] == a["bars"])
        ok("texts identical", b["texts"] == a["texts"])
        ok("x and y limits identical", np.allclose(b["xlim"], a["xlim"], rtol=0, atol=1e-12) and np.allclose(b["ylim"], a["ylim"], rtol=0, atol=1e-12), f"x {a['xlim']} y {a['ylim']}")
        ok("x ticks and labels identical", b["xticks"] == a["xticks"] and [TICK_RENAME.get(t, t) for t in b["xticklabels"]] == a["xticklabels"], str(a["xticklabels"]))
        ok("y tick positions identical", np.allclose(b["yticks"], a["yticks"], rtol=0, atol=1e-12))
        ok("axis labels identical", (b["xlabel"], b["ylabel"]) == (a["xlabel"], a["ylabel"]), a["ylabel"])
    return rec


def v30_reference(c9, abl):
    """Run the producer drawing unchanged (panel a at 9-50) into a hidden figure; collect values."""
    rc = {k: v for k, v in mpl.rcParamsDefault.items() if k not in ("backend", "backend_fallback")}
    rc.update(PRODUCER_RC)
    with mpl.rc_context(rc):
        fig, axes = plt.subplots(2, 2, figsize=FIGSIZE)
        axes = axes.ravel()
        panel_a(axes[0], c9)
        for f, ax in zip([panel_b, panel_c, panel_d], axes[1:]):
            f(ax, abl)
        fig.tight_layout(rect=(0, 0, 1, .96), pad=1.6)
        vals = v30_collect(axes)
        plt.close(fig)
    return vals


def v30_dots(ax, groups, ylabel, colours):
    """Producer ``dots()`` in house style: same x (linspace jitter), values and median bars."""
    for i, (label, values) in enumerate(groups):
        v = np.asarray(values, dtype=float)
        v = v[np.isfinite(v)]
        x = i + np.linspace(-.16, .16, len(v))
        ax.scatter(x, v, s=V30_DOT_S, alpha=V30_DOT_A, color=colours[i], linewidths=0, rasterized=len(v) > 100)
        if len(v):
            ax.hlines(np.median(v), i - .25, i + .25, color=fl.INK, lw=1.5)
    ax.set_xticks(range(len(groups)), [g[0] for g in groups])
    ax.set_ylabel(ylabel)


def render_v30(sd, out):
    """Draw the v30 figure (a at 9-50 residues; b-d unchanged) and check values against the producer.

    Writes ``SupplementaryFigure8.pdf/.png`` (300 dpi) and ``values_check.json`` to ``out``.

    Returns
    -------
    dict
        ``checks``, ``pdf``, ``png``, ``counts`` (9-50 long table), ``abl``, ``axes_bbox_px``.
    """
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    c9, abl = read_tables_v30(sd)
    before = v30_reference(c9, abl)
    fl.init()
    fig, axes = plt.subplots(2, 2, figsize=fl.mm(V30_W_MM, V30_H_MM))
    axes = axes.ravel()
    fig.subplots_adjust(left=0.08, right=0.985, bottom=0.04, top=0.965, wspace=0.26, hspace=0.30)
    arms = [(k, g) for k, g in c9.groupby("arm")]
    v30_dots(axes[0], [(ARM_TICK[k], g.peptides / 1000) for k, g in arms], "Accepted peptides (thousands)",
             [ARM_COLOUR[k] for k, _ in arms])
    for i, ax in enumerate(axes[1:]):
        col, label = ABLATION_COLUMNS[i]
        grp = [(co, g[col]) for co, g in abl.groupby("cohort")]
        v30_dots(ax, grp, label, [COHORT_COLOUR[co] for co, _ in grp])
    for ax, L, b in zip(axes, "abcd", before):
        fl.style(ax)
        ax.tick_params(labelsize=V30_TICK_PT)
        ax.xaxis.label.set_size(V30_LAB_PT)
        ax.yaxis.label.set_size(V30_LAB_PT)
        ax.set_yticks(b["yticks"])           # keep the producer tick values (locator depends on pixels)
        ax.set_xlim(b["xlim"])
        ax.set_ylim(b["ylim"])
        fl.panel(ax, L)
    after = v30_collect(axes)
    rec = v30_compare(before, after)
    legend_home = {"a": "legend a 'four PREDICT database arms ... curated specimen MG'", "b": "legend b-d 'intensity removed'",
                   "c": "legend b-d 'KEGG orthologs (KOs) losing all support'", "d": "legend b-d 'Bray-Curtis distance between ... KO profiles'"}
    for L, b, a in zip("abcd", before, after):
        rec.append(dict(panel=L, check="producer title removed (content already in the legend)", passed=a["title"] in ("", L),
                        detail=f"{b['title']!r} -> {legend_home[L]}"))
    (out / "values_check.json").write_text(json.dumps(rec, indent=1) + "\n")
    bad = [r for r in rec if not r["passed"]]
    if bad:
        raise SystemExit("VALUE CHECK FAILED: " + json.dumps(bad, indent=1))
    pdf, png = fl.save(fig, out / "SupplementaryFigure8", formats=("pdf", "png"))
    fig.set_dpi(COMPOSITE_DPI)
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    H = fig.bbox.height
    boxes = {}
    for letter, ax in zip("abcd", axes):
        bb = ax.get_tightbbox(r)
        boxes[letter] = [max(0, int(np.floor(bb.x0))), max(0, int(np.floor(H - bb.y1))),
                         int(np.ceil(bb.x1)), int(np.ceil(H - bb.y0))]
    plt.close(fig)
    return dict(checks=rec, pdf=Path(pdf), png=Path(png), counts=c9, abl=abl, axes_bbox_px=boxes)


def crop_px(png, box, out):
    """Crop ``png`` to the pixel box (left, top, right, bottom) and write ``out``."""
    from PIL import Image
    with Image.open(png) as im:
        w, h = im.size
        l, t, r, b = box
        im.crop((l, t, min(r, w), min(b, h))).save(out)
