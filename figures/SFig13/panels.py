"""Supplementary Figure 13: histidine-boundary sequence context and PREDICT replication.

``figure_v30`` draws the figure as placed from manuscript v30 on: the same values in the house
style (port of the v30 restyle script ``build_sf13_style.py`` and ``sf_common.py``).

Panels
------
a-f  CapScan histidine-conditioned motif deviations (proximal fluid, distal fluid, stool;
     N- and C-terminal boundaries), one heat map per panel, dots at q < 0.05.
g,h  PREDICT histidine enrichment (de novo top 30%, non-specific Sage with both q gates),
     93 patient means per boundary side, medians, means with 95% patient-bootstrap intervals.

The figure is drawn as one matplotlib figure, exactly as the original producer drew it, and
then the in-image title is redacted and the page rasterised at 300 dpi.

Source data
-----------
``source_data/Fig6/capscan_histidine_statistical_families.tsv``
    A28 CapScan complete statistical families (347 acquisitions, 15 participants); owned by
    Figure 6, read here. Source of the 840 displayed motif cells and their q values.
``source_data/SFig13/capscan_histidine_motif_display.tsv``
    The 840 displayed cells, as exported by the producer (written by ``derive_display_table``).
``source_data/SFig13/capscan_participant_motif_differences.tsv.gz``
    Participant-region motif differences from which the 32,400-test motif family is recomputed.
``source_data/SFig13/predict_histidine_acquisition_enrichment.tsv.gz``
    PREDICT acquisition-level observed and source-expected residue counts, pseudonymised.
``source_data/SFig13/predict_histidine_statistical_families.tsv``
    A28 PREDICT complete statistical families (1,560 source tests and 12 reciprocal tests).

Original producers (moved here with paths made repository-relative, logic unchanged):
``build_sf13.py`` (filters of ``review_histidine_predict.py`` l.16 and
``plot_histidine_review.py`` l.40, and the HISTIDINE_CONTEXT_REPLICATION block of
``build_revision.py`` l.254-275 with its ``dots``, ``title``, ``check`` and ``save``), and
``redact_sf13.py`` (the logic of ``prepare_supplement_figures.py``). Statistics are moved
from the A28 code (``histidine_helpers.signflip``, ``histidine_helpers.bh``, and the motif
and source-enrichment loops of ``collect_revision4.py``).
"""

import collections
import hashlib
import itertools
import re
from pathlib import Path

import matplotlib as mpl
mpl.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# ----------------------------------------------------------------------------- files
F_CAPSCAN_FAMILIES = ("Fig6", "capscan_histidine_statistical_families.tsv")  # shared, owned by Fig6
F_CAPSCAN_DISPLAY = "capscan_histidine_motif_display.tsv"
F_CAPSCAN_PARTICIPANT_MOTIF = "capscan_participant_motif_differences.tsv.gz"
F_PREDICT_ACQUISITION = "predict_histidine_acquisition_enrichment.tsv.gz"
F_PREDICT_FAMILIES = "predict_histidine_statistical_families.tsv"

SHA256 = {
    F_CAPSCAN_FAMILIES[1]: "3c24a42eec725e0b43bdeecca3928f971925c25395dde7478d92a1db9e1c28d7",
    F_CAPSCAN_DISPLAY: "269fa990cd52029a6dd4961e7bf21e347436c374c89bee4fa5c9631dd87364ea",
    F_CAPSCAN_PARTICIPANT_MOTIF: "df84cc2720b0ed63448048ef0676ad33c44ef9fe46eb05414a1d98a3ecb5ae77",
    F_PREDICT_ACQUISITION: "b40ea94ad1ba55f9931bce83035728cd7429f4d4c5d20733c08733beb3bfb951",
    F_PREDICT_FAMILIES: "85c774cb604a1d9e6c9526b04750b10d6601a07855fe1552357f66be016d28d7",
}
#: sha256 of the producer's exported CapScan all-tests Source Data table (pandas re-export of
#: the families table); ``export_all_tests`` must reproduce it byte for byte.
SHA256_CAPSCAN_ALL_TESTS_EXPORT = "2adeeaabe990c2f5635e6fc58b53e6e9bac9745766e3302554a4f726ff200088"

# ----------------------------------------------------------------------------- choices
#: Routes of the primary CapScan flank table (review_histidine_predict.py l.16).
MOTIF_PRIMARY_ROUTES = ["denovo_top30", "sage_semi_bothq", "sage_unspecific_bothq"]
#: Displayed route and origin (plot_histidine_review.py l.40): de novo top-30% predictions,
#: all origins pooled (Supplementary Methods 9, v28d line 867).
DISPLAY_ROUTE = "denovo_top30"
DISPLAY_ORIGIN = "All"
#: Motif population: boundaries whose P1 residue is histidine (legend, v28d line 461).
DISPLAY_POPULATION = "H_only"
#: The fixed H-conditioning P1 position, omitted from display and tests (v28d lines 461, 873).
CONDITIONED_POSITION = -1
#: Significance threshold for the dots, BH q within the complete family (v28d lines 461, 889).
Q_CUTOFF = 0.05
#: Panel order a-f: compartment x boundary side (legend, v28d line 461).
COMPARTMENTS = ["proximal", "distal", "stool"]
SIDES = ["N", "C"]
#: Panels g,h routes (legend, v28d line 461; Supplementary Methods 9, v28d line 869).
PREDICT_ROUTES = ["denovo_top30", "sage_unspecific_bothq"]
PREDICT_RESIDUE = "H"
PREDICT_ORIGIN = "All"

#: Producer drawing settings (build_revision.py l.13-14 and the histidine block).
PRODUCER_RC = {"font.family": "DejaVu Sans", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
               "xtick.labelsize": 7, "ytick.labelsize": 7, "pdf.fonttype": 42, "svg.fonttype": "none",
               "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 180}
COL = ["#2878a3", "#df8f44", "#619b68", "#b35872", "#8277ac"]
FIGSIZE_IN = (7.1, 10.8)
HEIGHT_RATIOS = [1.15, 1.15, 1.15, .8]
SUPTITLE = "Histidine context and replication across cohorts"   # redacted before release
RASTER_DPI = 300                                               # prepare_supplement_figures.py
#: Margin added around each panel box when cutting panels out of the composite (display and
#: per-panel pixel checks only; the composite itself is checked whole).
PANEL_CROP_MARGIN_IN = 0.06

# ----------------------------------------------------------------------------- statistics
#: A28 histidine_helpers.signflip constants (Supplementary Methods 9, v28d line 869; line 295).
SIGNFLIP_MIN_N = 8               # n < 8: not estimable, p = 1
SIGNFLIP_EXACT_MAX_N = 16        # n <= 16: all 2^n sign assignments enumerated
SIGNFLIP_DRAWS = 200_000         # n > 16: random sign assignments, plus-one correction
SIGNFLIP_SEED_BASE = 20260914    # seed = base + n (not stated in the text)
BOOTSTRAP_DRAWS = 20_000         # whole-unit bootstrap resamples for the mean
BOOTSTRAP_SEED_BASE = 20261014   # seed = base + n (not stated in the text)
#: A28 constants for the families (collect_revision4.py).
AA = "ACDEFGHIKLMNPQRSTVWY"
REGIONS = ["duodenum", "jejunum", "ileum", "colon", "stool"]
MOTIF_ROUTES = ["denovo_score09", "denovo_top30", "sage_semi", "sage_semi_bothq",
                "sage_unspecific", "sage_unspecific_bothq"]
MOTIF_ORIGINS = ["All", "Human", "Microbial"]
MOTIF_POPULATIONS = ["all_nonKR", "H_only"]
MOTIF_POSITIONS = [-4, -3, -2, -1, 1, 2, 3, 4]
MOTIF_COMPARTMENTS = [("proximal", REGIONS[:2]), ("distal", REGIONS[2:4]), ("stool", ["stool"])]
PREDICT_SOURCE_ROUTES = ["score_0_0.5", "score_0.5_0.7", "score_0.7_0.9", "score_0.9_0.95",
                         "score_0.95_0.99", "score_0.99_plus", "denovo_score09",
                         "denovo_score09_eligible", "denovo_top30", "sage_tryptic",
                         "sage_tryptic_bothq", "sage_unspecific", "sage_unspecific_bothq"]
PREDICT_CONTRAST = "ALL_observed_specimens"

_DEFAULT_RC = {k: v for k, v in mpl.rcParamsDefault.items() if k != "backend"}


# ============================================================================ helpers
def sha256(path):
    """Return the hex sha256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def check_sha256(path, expected):
    """Assert that ``path`` has sha256 ``expected`` (a full hex digest or a prefix)."""
    got = sha256(path)
    assert got.startswith(expected), (Path(path).name, got, expected)
    return got


# ============================================================================ data
def derive_display_table(families):
    """Select the 840 displayed motif cells from the CapScan family table.

    Moved from ``review_histidine_predict.py`` l.16 (primary H-conditioned flank table) and
    ``plot_histidine_review.py`` l.40 (displayed route and origin). The producer wrote the
    intermediate table to disk and read it back; the round trip is kept, because it fixes the
    dtypes and float text of the exported Source Data table.

    Parameters
    ----------
    families : pathlib.Path
        ``capscan_histidine_statistical_families.tsv``.

    Returns
    -------
    flanks : pandas.DataFrame
        Motif rows, H-conditioned population, three primary routes, P1 omitted.
    display : pandas.DataFrame
        The 840 displayed cells (route de novo top 30%, origin All), read back from text.
    display_text : str
        The display table as the producer wrote it (tab-separated text).
    """
    import io
    mot = pd.read_csv(families, sep="\t")
    mot = mot[mot.endpoint.eq("motif") & mot.population.eq(DISPLAY_POPULATION)
              & mot.route.isin(MOTIF_PRIMARY_ROUTES) & mot.position.ne(CONDITIONED_POSITION)].copy()
    buf = io.StringIO(); mot.to_csv(buf, sep="\t", index=False); buf.seek(0)
    m = pd.read_csv(buf, sep="\t")
    m = m[m.route.eq(DISPLAY_ROUTE) & m.origin.eq(DISPLAY_ORIGIN)]
    buf = io.StringIO(); m.to_csv(buf, sep="\t", index=False); buf.seek(0)
    return mot, pd.read_csv(buf, sep="\t"), buf.getvalue()


def export_all_tests(families, out):
    """Write the producer's ``CapScan_histidine_all_tests.tsv`` export (read, then ``to_csv``).

    Parameters
    ----------
    families : pathlib.Path
        ``capscan_histidine_statistical_families.tsv``.
    out : pathlib.Path
        Output file (under ``build/``).

    Returns
    -------
    pandas.DataFrame
        The complete CapScan family table.
    """
    allhist = pd.read_csv(families, sep="\t")
    allhist.to_csv(out, sep="\t", index=False)
    return allhist


def predict_patient_means(acquisition):
    """PREDICT patient means of observed and source-expected residue fractions.

    Moved from ``build_revision.py`` l.258-259: acquisitions are averaged within specimen and
    specimens within patient; rows without a patient are dropped.

    Parameters
    ----------
    acquisition : pathlib.Path
        ``predict_histidine_acquisition_enrichment.tsv.gz``.

    Returns
    -------
    pandas.DataFrame
        One row per route, side, origin, residue and patient (statistical unit: patient,
        n = 93), columns ``O``, ``E`` and ``difference`` = O - E (fractions).
    """
    ef = pd.read_csv(acquisition, sep="\t")
    ef = ef[ef.participant.notna() & ef.participant.astype(str).ne("")].copy()
    key = ["route", "side", "origin", "residue"]
    ef["O"] = ef.observed / ef.denominator.replace(0, np.nan)
    ef["E"] = ef.expected / ef.denominator.replace(0, np.nan)
    sp = ef.groupby(key + ["cohort", "participant", "region", "specimen_id"], dropna=False)[["O", "E"]].mean().reset_index()
    pt = sp.groupby(key + ["cohort", "participant"], dropna=False)[["O", "E"]].mean().reset_index()
    pt["difference"] = pt.O - pt.E
    return pt


# ============================================================================ drawing
def check(name, ok, details):
    """Producer assertion helper (build_revision.py ``check``)."""
    assert ok, (name, details)


def dots(ax, groups, ylabel, percent=False):
    """Producer strip plot with median bars (build_revision.py ``dots``)."""
    for i, (label, values) in enumerate(groups):
        v = np.asarray(values, dtype=float); v = v[np.isfinite(v)] * (100 if percent else 1)
        x = i + np.linspace(-.16, .16, len(v)); ax.scatter(x, v, s=8, alpha=.48, color=COL[i % 5], rasterized=len(v) > 100)
        if len(v):
            ax.hlines(np.median(v), i - .25, i + .25, color="black", lw=1.5)
    ax.set_xticks(range(len(groups)), [g[0] for g in groups]); ax.set_ylabel(ylabel)


def title(ax, l, t):
    """Producer panel title (build_revision.py ``title``)."""
    ax.set_title(l + "  " + t, loc="left", pad=9)


def color_scale(display):
    """Symmetric colour limit: largest absolute displayed cell mean, in percentage points."""
    return float(np.abs(display["mean"]).max() * 100)


def _motif_panel(fig, ax, display, compartment, side, letter, scale):
    """Draw one motif heat map (build_revision.py l.264-267) and return its dot summary."""
    g = display[display.contrast.eq(compartment) & display.side.eq(side)]
    a = g.pivot(index="residue", columns="position", values="mean").sort_index()
    qmat = g.pivot(index="residue", columns="position", values="q").reindex_like(a)
    im = ax.imshow(a.values * 100, vmin=-scale, vmax=scale, cmap="RdBu_r", aspect="auto")
    ax.set_xticks(range(len(a.columns)), [f"{int(x):+}" for x in a.columns])
    ax.set_yticks(range(len(a)), a.index, fontsize=7)
    ax.set_xlabel("Position relative to boundary")
    title(ax, letter, compartment.capitalize() + "; " + side + " boundary")
    yy, xx = np.where(qmat.values < Q_CUTOFF)
    ax.scatter(xx, yy, s=3, c="black")
    fig.colorbar(im, ax=ax, fraction=.04, pad=.02, label="Observed − expected (pp)")
    return dict(region=compartment, side=side, cells=int(a.size), q_lt05=len(xx),
                n_min=int(g.n.min()), n_max=int(g.n.max()))


def panel_a(fig, ax, display, scale):
    """Panel a: CapScan proximal fluid, N-terminal boundary, H-conditioned motif deviations.

    Parameters
    ----------
    fig, ax : matplotlib Figure and Axes
        The figure and its first axes.
    display : pandas.DataFrame
        The 840 displayed cells (``derive_display_table``).
    scale : float
        Colour limit in percentage points (``color_scale``).

    Returns
    -------
    dict
        Cells, number of dots (q < 0.05) and participant n range (unit: participant, n = 15).
    """
    return _motif_panel(fig, ax, display, "proximal", "N", "a", scale)


def panel_b(fig, ax, display, scale):
    """Panel b: proximal fluid, C-terminal boundary. Parameters and returns as ``panel_a``."""
    return _motif_panel(fig, ax, display, "proximal", "C", "b", scale)


def panel_c(fig, ax, display, scale):
    """Panel c: distal fluid, N-terminal boundary. Parameters and returns as ``panel_a``."""
    return _motif_panel(fig, ax, display, "distal", "N", "c", scale)


def panel_d(fig, ax, display, scale):
    """Panel d: distal fluid, C-terminal boundary. Parameters and returns as ``panel_a``."""
    return _motif_panel(fig, ax, display, "distal", "C", "d", scale)


def panel_e(fig, ax, display, scale):
    """Panel e: stool, N-terminal boundary. Parameters and returns as ``panel_a``."""
    return _motif_panel(fig, ax, display, "stool", "N", "e", scale)


def panel_f(fig, ax, display, scale):
    """Panel f: stool, C-terminal boundary. Parameters and returns as ``panel_a``."""
    return _motif_panel(fig, ax, display, "stool", "C", "f", scale)


def _predict_panel(ax, pt, st, j):
    """Draw one PREDICT panel (build_revision.py l.269-274)."""
    route = PREDICT_ROUTES[j]
    g = pt[pt.route.eq(route) & pt.origin.eq(PREDICT_ORIGIN) & pt.residue.eq(PREDICT_RESIDUE)]
    dots(ax, [(s, g[g.side.eq(s)].difference) for s in SIDES], "H enrichment (percentage points)", True)
    out = {}
    for i, side in enumerate(SIDES):
        r = st[st.route.eq(route) & st.origin.eq(PREDICT_ORIGIN) & st.residue.eq(PREDICT_RESIDUE)
               & st.side.eq(side) & st.endpoint.eq("source_enrichment")].iloc[0]
        vals = g[g.side.eq(side)].difference.dropna()
        check("PREDICT H participant mean " + route + side, len(vals) == r.n and np.isclose(vals.mean(), r["mean"]),
              "Specimen means then equal participant means; n=" + str(len(vals)))
        ax.errorbar(i + .28, r["mean"] * 100, yerr=[[100 * (r["mean"] - r.CI95_low)], [100 * (r.CI95_high - r["mean"])]],
                    fmt="s", c="black", ms=3, capsize=3)
        out[side] = dict(n=len(vals), median_pp=float(np.median(vals) * 100), mean_pp=float(r["mean"] * 100),
                         ci_pp=(float(r.CI95_low * 100), float(r.CI95_high * 100)))
    ax.axhline(0, c="grey", lw=.6)
    title(ax, "gh"[j], "PREDICT: " + ("de novo top 30%" if j == 0 else "unspecific Sage; both q gates"))
    return out


def panel_g(ax, pt, st):
    """Panel g: PREDICT histidine enrichment, de novo top-30% predictions, N and C boundaries.

    Parameters
    ----------
    ax : matplotlib Axes
    pt : pandas.DataFrame
        Patient means (``predict_patient_means``); statistical unit patient, n = 93 per side.
    st : pandas.DataFrame
        PREDICT families (``predict_histidine_statistical_families.tsv``), source of the plotted
        means and 95% patient-bootstrap intervals.

    Returns
    -------
    dict
        Per side: n, median, mean and interval in percentage points.
    """
    return _predict_panel(ax, pt, st, 0)


def panel_h(ax, pt, st):
    """Panel h: as ``panel_g`` for non-specific Sage with peptide and spectrum q gates."""
    return _predict_panel(ax, pt, st, 1)


def draw_figure(display, pt, st, out_pdf):
    """Draw the complete figure as the producer did and save the unredacted PDF.

    Drawing runs under matplotlib defaults plus ``PRODUCER_RC`` (the producer did not use
    the house style).

    Parameters
    ----------
    display, pt, st : pandas.DataFrame
        Displayed motif cells, PREDICT patient means, PREDICT families.
    out_pdf : pathlib.Path
        Output PDF (with the in-image title, as the producer's
        ``HISTIDINE_CONTEXT_REPLICATION.pdf``).

    Returns
    -------
    dict
        ``scale`` (colour limit, pp), per-panel summaries ``a``..``h`` and ``panel_boxes_in``
        (each panel's bounding box in inches from the top-left, for panel crops).
    """
    with mpl.rc_context(_DEFAULT_RC):
        plt.rcParams.update(PRODUCER_RC)
        fig, axes = plt.subplots(4, 2, figsize=FIGSIZE_IN, gridspec_kw={"height_ratios": HEIGHT_RATIOS})
        axes = axes.ravel()
        scale = color_scale(display)
        summary = {"scale": scale}
        for i, f in enumerate([panel_a, panel_b, panel_c, panel_d, panel_e, panel_f]):
            summary["abcdef"[i]] = f(fig, axes[i], display, scale)
        summary["g"] = panel_g(axes[6], pt, st)
        summary["h"] = panel_h(axes[7], pt, st)
        # save() of build_revision.py
        fig.suptitle(SUPTITLE, fontsize=11, fontweight="bold")
        fig.tight_layout(rect=(0, 0, 1, .96), pad=1.6)
        fig.savefig(out_pdf)
        summary["panel_boxes_in"] = _panel_boxes(fig, axes)
        plt.close(fig)
    return summary


def _panel_boxes(fig, axes):
    """Bounding box of each panel (axes, title, labels and colour bar) in inches from top-left."""
    r = fig.canvas.get_renderer()
    H = fig.get_figheight()
    cbars = [a for a in fig.axes if a not in list(axes)]
    boxes = {}
    for i, ax in enumerate(axes):
        bbs = [ax.get_tightbbox(r), ax.yaxis.label.get_window_extent(r), ax.xaxis.label.get_window_extent(r),
               ax.title.get_window_extent(r), ax._left_title.get_window_extent(r)]
        bbs += [cbars[i].get_tightbbox(r)] if i < 6 else []
        x0 = min(b.x0 for b in bbs) / fig.dpi; x1 = max(b.x1 for b in bbs) / fig.dpi
        y0 = min(b.y0 for b in bbs) / fig.dpi; y1 = max(b.y1 for b in bbs) / fig.dpi
        boxes["abcdefgh"[i]] = (x0, H - y1, x1, H - y0)
    return boxes


def redact_and_rasterise(pdf_in, pdf_out, dpi=RASTER_DPI):
    """Remove the in-image title and render the page (logic of prepare_supplement_figures.py).

    The title box is located by text search; numeric labels are asserted unchanged.

    Parameters
    ----------
    pdf_in : pathlib.Path
        PDF from ``draw_figure``.
    pdf_out : pathlib.Path
        Redacted PDF; the PNG is written beside it with suffix ``.png``.
    dpi : int
        Raster resolution.

    Returns
    -------
    dict
        Title bounding box (pt) and page rectangle.
    """
    import fitz
    doc = fitz.open(pdf_in); p = doc[0]; before = p.get_text(); ttext = SUPTITLE
    spans = [sp for b in p.get_text("dict")["blocks"] for l in b.get("lines", []) for sp in l["spans"]
             if sp["text"].strip() == ttext]
    assert len(spans) == 1
    bbox = spans[0]["bbox"]; assert bbox[1] < 20
    p.add_redact_annot(fitz.Rect(bbox), fill=(1, 1, 1)); p.apply_redactions(images=0, graphics=0)
    assert ttext not in p.get_text()
    before_nums = collections.Counter(re.findall(r"\d[\d,.]*", before.replace(ttext, "")))
    after_nums = collections.Counter(re.findall(r"\d[\d,.]*", p.get_text()))
    assert before_nums == after_nums
    doc.save(pdf_out, garbage=4, deflate=True)
    p = fitz.open(pdf_out)[0]; p.get_pixmap(dpi=dpi, alpha=False).save(Path(pdf_out).with_suffix(".png"))
    return dict(title_bbox=[round(v, 3) for v in bbox], rect=list(p.rect))


def crop_panels(png, boxes_in, out_dir, prefix, dpi=RASTER_DPI):
    """Cut each panel out of a rendered composite PNG.

    Parameters
    ----------
    png : pathlib.Path
        Composite rendered at ``dpi``.
    boxes_in : dict
        Panel boxes in inches (``draw_figure``).
    out_dir : pathlib.Path
    prefix : str
        File name prefix, e.g. ``panel_`` or ``expected_panel_``.

    Returns
    -------
    dict
        Panel letter to written path.
    """
    from PIL import Image
    im = Image.open(png)
    out = {}
    for k, (x0, y0, x1, y1) in boxes_in.items():
        m = PANEL_CROP_MARGIN_IN
        box = (max(0, int(np.floor((x0 - m) * dpi))), max(0, int(np.floor((y0 - m) * dpi))),
               min(im.width, int(np.ceil((x1 + m) * dpi))), min(im.height, int(np.ceil((y1 + m) * dpi))))
        f = Path(out_dir) / f"{prefix}{k}.png"
        im.crop(box).save(f)
        out[k] = f
    return out


# ============================================================================ statistics
_null_cache = {}
_boot_cache = {}


def signflip(delta):
    """Two-sided paired sign-flip test of mean zero, with a whole-unit bootstrap interval.

    Moved unchanged from A28 ``histidine_helpers.signflip``. For n <= 16 all 2^n sign
    assignments are enumerated (exact P); for larger n, 200,000 random assignments
    (seed 20260914 + n) with a plus-one correction. The statistic is |mean|; P counts
    assignments with |mean| at least the observed |mean| (minus 1e-12). The 95% interval is
    the 2.5 and 97.5 percentiles of 20,000 multinomial bootstrap means (seed 20261014 + n).
    n < 8 returns p = 1, not estimable.

    Parameters
    ----------
    delta : array-like
        One value per statistical unit (participant or patient); non-finite values dropped.

    Returns
    -------
    dict
        n, mean, p, estimable, CI95_low, CI95_high, n_positive, n_negative.
    """
    d = np.asarray(delta, float); d = d[np.isfinite(d)]; n = len(d)
    if n < SIGNFLIP_MIN_N:
        return dict(n=n, mean=float(d.mean()) if n else np.nan, p=1., estimable=False)
    rng = np.random.default_rng(SIGNFLIP_SEED_BASE + n)
    if n not in _null_cache:
        if n <= SIGNFLIP_EXACT_MAX_N:
            _null_cache[n] = ((np.arange(2 ** n, dtype=np.uint32)[:, None] >> np.arange(n)) & 1).astype(np.int8) * 2 - 1
        else:
            _null_cache[n] = rng.choice(np.array([-1, 1], np.int8), size=(SIGNFLIP_DRAWS, n))
    signs = _null_cache[n]
    hits = int((np.abs(signs @ d / n) >= abs(d.mean()) - 1e-12).sum())
    p = hits / len(signs) if n <= SIGNFLIP_EXACT_MAX_N else (hits + 1) / (len(signs) + 1)
    if n not in _boot_cache:
        _boot_cache[n] = np.random.default_rng(BOOTSTRAP_SEED_BASE + n).multinomial(n, np.full(n, 1 / n), size=BOOTSTRAP_DRAWS)
    ci = np.quantile(_boot_cache[n] @ d / n, [.025, .975])
    return dict(n=n, mean=float(d.mean()), p=p, estimable=True, CI95_low=ci[0], CI95_high=ci[1],
                n_positive=int((d > 0).sum()), n_negative=int((d < 0).sum()))


def bh(p):
    """Benjamini-Hochberg adjusted P values (moved unchanged from A28 ``histidine_helpers.bh``)."""
    p = np.asarray(p, float); order = np.argsort(p); v = np.empty(len(p))
    v[order] = np.minimum(1, np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    return v


def stat_motif_family(participant_motif):
    """Recompute the complete CapScan motif family (32,400 tests) from participant values.

    Moved from ``collect_revision4.py`` (motif loop): for every route x side x origin x
    population x flank position x residue (P1 omitted for the H-conditioned population) and
    each compartment, the participant value is the mean of that participant's region values
    (proximal: duodenum and jejunum; distal: ileum and colon; stool), with a missing region
    making the participant missing. Each cell is tested with ``signflip`` (statistical unit:
    participant, n = 15 where all regions are present); BH over all 32,400 tests.
    Participant identifiers are read as text, as they were held in the original run, so
    participants are ordered "1", "10", ..., "15", "2", ...; this order fixes which bootstrap
    weight falls on which participant and so reproduces the published intervals exactly.

    Parameters
    ----------
    participant_motif : pathlib.Path
        ``capscan_participant_motif_differences.tsv.gz``.

    Returns
    -------
    pandas.DataFrame
        One row per test with key columns, n, mean, p, CI95_low, CI95_high, n_positive,
        n_negative and q.
    """
    pm = pd.read_csv(participant_motif, sep="\t", dtype={"participant": str})
    key = ["route", "side", "origin", "population", "position", "residue"]
    cells = {labels: g for labels, g in pm.groupby(key, sort=False)}
    tests = []
    for labels in itertools.product(MOTIF_ROUTES, SIDES, MOTIF_ORIGINS, MOTIF_POPULATIONS, MOTIF_POSITIONS, AA):
        if labels[3] == "H_only" and labels[4] == CONDITIONED_POSITION:
            continue
        g = cells.get(labels)
        pv = g.pivot(index="participant", columns="region", values="difference").reindex(columns=REGIONS) if g is not None else pd.DataFrame(columns=REGIONS, dtype=float)
        for compartment, columns in MOTIF_COMPARTMENTS:
            d = pv[columns].mean(axis=1, skipna=False)
            tests.append(dict(endpoint="motif", family="all_motif_positions_residues", **dict(zip(key, labels)),
                              contrast=compartment, **signflip(d)))
    tf = pd.DataFrame(tests)
    tf["q"] = bh(tf.p)
    return tf


def stat_predict_source_family(acquisition):
    """Recompute the PREDICT source-enrichment family (1,560 tests) from acquisition values.

    Moved from ``collect_revision4.py`` (``participant(..., all_specimens=True)`` and the
    source-enrichment loop): acquisition fractions O and E are averaged within specimen,
    specimens within patient; the patient value is O - E. Each route x side x origin x
    residue cell is tested with ``signflip`` (statistical unit: patient, n = 93 for the
    displayed cells); BH over the 1,560 tests. Planned cells without data enter with p = 1.

    Parameters
    ----------
    acquisition : pathlib.Path
        ``predict_histidine_acquisition_enrichment.tsv.gz``.

    Returns
    -------
    pandas.DataFrame
        One row per test with route, side, origin, residue, contrast, n, mean, p, CI95_low,
        CI95_high, n_positive, n_negative and q.
    """
    ef = pd.read_csv(acquisition, sep="\t", dtype={"participant": str, "region": str, "specimen_id": str},
                     keep_default_na=False, na_values={"observed": [""], "expected": [""], "denominator": [""]})
    for k in ["observed", "expected", "denominator"]:
        ef[k] = pd.to_numeric(ef[k], errors="coerce")
    ef["O"] = ef.observed / ef.denominator.replace(0, np.nan)
    ef["E"] = ef.expected / ef.denominator.replace(0, np.nan)
    key = ["route", "side", "origin", "residue"]

    def participant(frame, value):
        f = frame.copy(); f[value] = pd.to_numeric(f[value], errors="coerce"); f = f[f.participant.astype(str).ne("")]
        sp = f.groupby(key + ["cohort", "participant", "region", "specimen_id"], dropna=False)[value].mean().reset_index()
        return sp.groupby(key + ["cohort", "participant"], dropna=False)[value].mean().reset_index().assign(region=PREDICT_CONTRAST)

    pp = participant(ef, "O").merge(participant(ef, "E"), on=key + ["cohort", "participant", "region"], validate="1:1")
    pp["difference"] = pp.O - pp.E
    lookup = {labels: g for labels, g in pp.groupby(key + ["region"], sort=False)}
    tests = []
    for route, side, origin, residue in itertools.product(PREDICT_SOURCE_ROUTES, SIDES, MOTIF_ORIGINS, AA):
        g = lookup.get((route, side, origin, residue, PREDICT_CONTRAST), pd.DataFrame({"difference": []}))
        tests.append(dict(endpoint="source_enrichment", family="source_P1_all_residues", route=route, side=side,
                          origin=origin, residue=residue, contrast=PREDICT_CONTRAST, **signflip(g.difference)))
    tf = pd.DataFrame(tests)
    tf["q"] = bh(tf.p)
    return tf


def compare_family(recomputed, published, key, cols=("n", "mean", "p", "CI95_low", "CI95_high", "q")):
    """Compare a recomputed family with the published family table, row by row.

    Parameters
    ----------
    recomputed, published : pandas.DataFrame
    key : list of str
        Columns identifying a test.
    cols : sequence of str
        Columns compared.

    Returns
    -------
    dict
        Number of rows on each side, rows matched, and per column the number of exactly equal
        values and the largest absolute difference.
    """
    a = recomputed.set_index(key)
    b = published.set_index(key)
    common = a.index.intersection(b.index)
    out = dict(recomputed_rows=len(a), published_rows=len(b), matched_rows=len(common))
    for c in cols:
        x = pd.to_numeric(a.loc[common, c], errors="coerce").to_numpy(float)
        y = pd.to_numeric(b.loc[common, c], errors="coerce").to_numpy(float)
        both_nan = np.isnan(x) & np.isnan(y)
        diff = np.where(both_nan, 0, np.abs(x - y))
        out[c] = dict(equal=int(((x == y) | both_nan).sum()), max_abs_diff=float(np.nanmax(diff)) if len(diff) else 0.0)
    return out


def stat_display_counts(display, families):
    """Counts printed in or implied by the legend: family size, displayed cells, dots.

    Parameters
    ----------
    display : pandas.DataFrame
        The 840 displayed cells.
    families : pandas.DataFrame
        CapScan complete families (or the recomputed motif family).

    Returns
    -------
    dict
        motif_family_tests, motif_family_q_lt05, display_cells, display_q_lt05 and the
        participant n values among displayed cells (unit: participant).
    """
    mf = families[families.family.eq("all_motif_positions_residues")]
    return dict(motif_family_tests=len(mf), motif_family_q_lt05=int((mf.q < Q_CUTOFF).sum()),
                display_cells=len(display), display_q_lt05=int((display.q < Q_CUTOFF).sum()),
                display_n_values=sorted(display.n.unique().tolist()))


def stat_predict_all_positive(pt):
    """Patients with positive histidine enrichment, per route and boundary side.

    Statistical unit: patient (n = 93). Supports "all 93 patients showed positive enrichment
    with both de novo predictions and non-specific database searches" (v28d line 175).

    Parameters
    ----------
    pt : pandas.DataFrame
        Patient means (``predict_patient_means``).

    Returns
    -------
    dict
        (route, side) -> dict(n, n_positive, n_negative, min_pp).
    """
    out = {}
    h = pt[pt.origin.eq(PREDICT_ORIGIN) & pt.residue.eq(PREDICT_RESIDUE)]
    for route in PREDICT_ROUTES:
        for side in SIDES:
            v = h[h.route.eq(route) & h.side.eq(side)].difference.dropna()
            out[(route, side)] = dict(n=int(len(v)), n_positive=int((v > 0).sum()), n_negative=int((v < 0).sum()),
                                     min_pp=float(v.min() * 100))
    return out


# ============================================================================ v30 house style
# Ported from the v30 restyle helpers (LAB_ROOT/.../final/v30/figures_style/SFig1x/sf_common.py):
# fingerprint/capture/compare prove the restyled figure plots the same values as the producer;
# _v30_finish is the fl_style save sequence. Logic unchanged; fl_style comes from REPO/style.
V30_LETTER_PT = 7.0        # sf_common.LETTER_PT: house panel letters, 7 pt
V30_DPI = 300              # PNG = the PDF rasterised at 300 dpi (sf_common.finish)


def _v30_fl():
    """Return the repository fl_style module with the house rcParams applied (sf_common.house).

    Returns
    -------
    module
        ``style.fl_style`` after ``matplotlib.rcdefaults()`` and ``fl.init()``.
    """
    from style import fl_style as fl
    import matplotlib as _mpl
    _mpl.rcdefaults()
    fl.init()
    return fl


def _v30_is_cbar(ax):
    return hasattr(ax, '_colorbar') or ax.get_label() == '<colorbar>'


def _v30_arr(a):
    a = np.ma.asarray(a)
    return np.ma.filled(a.astype(float), np.nan)


def v30_fingerprint(fig, skip_line=lambda l: False, text_norm=lambda s: s, data_text=lambda t: True):
    """Every plotted array of a figure (sf_common.fingerprint).

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Drawn figure.
    skip_line, text_norm, data_text : callable
        Lines to ignore, text normaliser, and which axes texts are data labels.

    Returns
    -------
    list of dict
        One dict per non-colour-bar axes: limits, in-view ticks and tick labels, images,
        collection offsets, line-collection segments, lines and data texts.
    """
    from matplotlib.collections import LineCollection
    out = []
    for ax in fig.axes:
        if _v30_is_cbar(ax):
            continue
        xl, yl = sorted(ax.get_xlim()), sorted(ax.get_ylim())
        xin = [k for k, v in enumerate(ax.get_xticks()) if xl[0] - 1e-9 <= v <= xl[1] + 1e-9]
        yin = [k for k, v in enumerate(ax.get_yticks()) if yl[0] - 1e-9 <= v <= yl[1] + 1e-9]
        xt, yt = ax.get_xticklabels(), ax.get_yticklabels()
        d = {'xlim': list(ax.get_xlim()), 'ylim': list(ax.get_ylim()),
             'xticks': [float(ax.get_xticks()[k]) for k in xin], 'yticks': [float(ax.get_yticks()[k]) for k in yin],
             'xticklabels': [xt[k].get_text() for k in xin],
             'yticklabels': [yt[k].get_text() for k in yin],
             'images': [], 'offsets': [], 'segments': [], 'lines': [], 'texts': []}
        for im in ax.images:
            d['images'].append({'array': _v30_arr(im.get_array()), 'clim': list(im.get_clim()),
                                'extent': list(im.get_extent())})
        for c in ax.collections:
            if isinstance(c, LineCollection):
                d['segments'].append([np.asarray(s, float) for s in c.get_segments()])
            else:
                d['offsets'].append(np.asarray(c.get_offsets(), float))
        for l in ax.lines:
            if skip_line(l):
                continue
            d['lines'].append(np.column_stack([np.asarray(l.get_xdata(), float), np.asarray(l.get_ydata(), float)]))
        for t in ax.texts:
            if data_text(t):
                d['texts'].append((text_norm(t.get_text()), float(t.get_position()[0])))
        d['texts'].sort()
        out.append(d)
    return out


def v30_capture(run, **fpkw):
    """Run ``run()`` with ``Figure.savefig`` intercepted; return each figure's fingerprint.

    Parameters
    ----------
    run : callable
        Draws and saves one or more figures (the producer).
    **fpkw
        Passed to ``v30_fingerprint``.

    Returns
    -------
    list
        One fingerprint per captured figure; nothing is written by the intercepted saves.
    """
    from matplotlib.figure import Figure
    caught = []
    orig = Figure.savefig

    def fake(self, *a, **k):
        self.canvas.draw()
        if not any(c[0] is self for c in caught):
            caught.append((self, v30_fingerprint(self, **fpkw)))
    Figure.savefig = fake
    try:
        run()
    finally:
        Figure.savefig = orig
    return [c[1] for c in caught]


def _v30_eq(a, b):
    a = np.asarray(a, float); b = np.asarray(b, float)
    return a.shape == b.shape and np.array_equal(a, b, equal_nan=True)


def v30_compare(ref, new, ignore=()):
    """Exact comparison of two fingerprints (sf_common.compare); raises on any difference.

    Parameters
    ----------
    ref, new : list of dict
        Producer and restyled fingerprints.
    ignore : tuple of str
        Keys not compared.

    Returns
    -------
    list of str
        One report line per axes, ending ``IDENTICAL``.
    """
    rep = []
    assert len(ref) == len(new), ('data axes', len(ref), len(new))
    for i, (r, n) in enumerate(zip(ref, new)):
        for k in ['xlim', 'ylim']:
            if k in ignore:
                continue
            assert _v30_eq(r[k], n[k]), (i, k, r[k], n[k])
        for k in ['xticks', 'yticks', 'xticklabels', 'yticklabels']:
            if k in ignore:
                continue
            assert r[k] == n[k], (i, k, r[k], n[k])
        assert len(r['images']) == len(n['images']), (i, 'images')
        for a, b in zip(r['images'], n['images']):
            assert _v30_eq(a['array'], b['array']) and _v30_eq(a['clim'], b['clim']) and _v30_eq(a['extent'], b['extent']), (i, 'image')
        assert len(r['offsets']) == len(n['offsets']), (i, 'n offsets', len(r['offsets']), len(n['offsets']))
        for a, b in zip(r['offsets'], n['offsets']):
            assert _v30_eq(a, b), (i, 'offsets')
        assert len(r['segments']) == len(n['segments']), (i, 'n segment sets')
        for a, b in zip(r['segments'], n['segments']):
            assert len(a) == len(b) and all(_v30_eq(x, y) for x, y in zip(a, b)), (i, 'segments')
        assert len(r['lines']) == len(n['lines']), (i, 'n lines', len(r['lines']), len(n['lines']))
        for a, b in zip(r['lines'], n['lines']):
            assert _v30_eq(a, b), (i, 'line')
        if 'texts' not in ignore:
            assert r['texts'] == n['texts'], (i, 'texts', r['texts'][:5], n['texts'][:5])
        npts = sum(x['array'].size for x in r['images']) + sum(len(x) for x in r['offsets']) \
            + sum(len(s) for x in r['segments'] for s in x) + sum(len(x) for x in r['lines'])
        rep.append(f"axes {i}: images {len(r['images'])}, point sets {len(r['offsets'])}, segment sets "
                   f"{len(r['segments'])}, lines {len(r['lines'])}, data texts {len(r['texts'])}, "
                   f"values compared {npts}: IDENTICAL")
    return rep


def _v30_finish(fl, fig, base, letters, keep_sizes=()):
    """fl_style save sequence of the v30 rebuild (sf_common.finish).

    enforce_nbt, normalise_type, panel letters re-set to ``V30_LETTER_PT``, optional text
    sizes restored, clip_check, PDF save, then the PDF rasterised at ``V30_DPI`` with PyMuPDF.

    Parameters
    ----------
    fl : module
        ``style.fl_style``.
    fig : matplotlib.figure.Figure
        Restyled figure.
    base : pathlib.Path
        Output path without suffix; writes ``.pdf`` and ``.png``.
    letters : list
        Axes carrying panel letters.
    keep_sizes : sequence
        Text artists whose font size must survive ``normalise_type``.

    Returns
    -------
    bool
        ``clip_check`` result.
    """
    import fitz
    base = str(base)
    fl.enforce_nbt(fig)
    sizes = [(t, t.get_fontsize()) for t in keep_sizes]
    fl.normalise_type(fig, Path(base).name, verbose=False)
    for ax in letters:
        ax._left_title.set_fontsize(V30_LETTER_PT)
    for t, s in sizes:
        t.set_fontsize(s)
    ok = fl.clip_check(fig, Path(base).name)
    fig.savefig(f'{base}.pdf', dpi=V30_DPI, bbox_inches=None, pad_inches=0, facecolor='white', edgecolor='none')
    plt.close(fig)
    pix = fitz.open(f'{base}.pdf')[0].get_pixmap(dpi=V30_DPI, alpha=False); pix.set_dpi(V30_DPI, V30_DPI); pix.save(f'{base}.png')
    return ok


# v30 house-style figure (LAB_ROOT/.../final/v30/figures_style/SFig13/build_sf13_style.py)
V30_SIZE_MM = (160, 243.4)                 # placed size in the v29b docx
V30_HEIGHT_RATIOS = [1.1, 1.1, 1.1, .95]
V30_W_PAD, V30_H_PAD, V30_HSPACE, V30_WSPACE = 2 / 72, 4 / 72, 0.03, 0.04   # constrained layout
V30_REGION_LABEL = {'proximal': 'Regions 1–2', 'distal': 'Regions 3–4'}     # v29 relabel; stool -> "Stool"
V30_ROUTE_LAB = {'denovo_top30': 'PREDICT: de novo top 30%',
                 'sage_unspecific_bothq': 'PREDICT: unspecific Sage;\nboth q gates'}   # v29 titles, verbatim


def _v30_dots(ax, groups, ylabel, col, ink, percent=False):
    """Producer strip plot (``dots``) with one route colour and ink median bars."""
    for i, (label, values) in enumerate(groups):
        v = np.asarray(values, dtype=float); v = v[np.isfinite(v)] * (100 if percent else 1)
        x = i + np.linspace(-.16, .16, len(v)); ax.scatter(x, v, s=8, alpha=.48, color=col, linewidths=0, rasterized=len(v) > 100)
        if len(v):
            ax.hlines(np.median(v), i - .25, i + .25, color=ink, lw=1.5)
    ax.set_xticks(range(len(groups)), [g[0] for g in groups]); ax.set_ylabel(ylabel)


def figure_v30(mot, pt, st, out_stem, ref_pdf):
    """Draw Supplementary Figure 13 in the v30 house style and prove its values unchanged.

    Port of ``build_sf13_style.py``: the producer (``draw_figure``) is run with savefig
    intercepted (reference fingerprint); the figure is redrawn with fl_style (Söhne,
    160 x 243.4 mm, ``fl.div_cmap`` at the same symmetric limits, panel titles a-f moved to the
    y-axis labels with the v29 "Regions 1–2"/"Regions 3–4" wording, g,h titles moved to the
    legend, route colours of main Figure 6: de novo top 30% blue, non-specific Sage green);
    every plotted array is compared exactly with the reference.

    Parameters
    ----------
    mot : pandas.DataFrame
        The 840 displayed motif cells (``capscan_histidine_motif_display.tsv``).
    pt : pandas.DataFrame
        PREDICT patient means (``predict_patient_means``), n = 93 patients per side.
    st : pandas.DataFrame
        PREDICT statistical families.
    out_stem : pathlib.Path
        Output path without suffix; ``.pdf`` and 300 dpi ``.png`` are written.
    ref_pdf : pathlib.Path
        Path passed to the intercepted producer (nothing is written there).

    Returns
    -------
    dict
        ``panels`` (comparison report lines), ``clip_ok``.
    """
    from matplotlib.lines import Line2D
    ref = v30_capture(lambda: draw_figure(mot, pt, st, ref_pdf))
    assert len(ref) == 1
    ref = ref[0]
    with mpl.rc_context():
        fl = _v30_fl(); INK = fl.INK
        ROUTE_COL = {'denovo_top30': fl.BLUE, 'sage_unspecific_bothq': fl.GREEN}   # main Figure 6 legend
        fig, axes = plt.subplots(4, 2, figsize=fl.mm(*V30_SIZE_MM), layout='constrained',
                                 gridspec_kw={'height_ratios': V30_HEIGHT_RATIOS})
        fig.get_layout_engine().set(w_pad=V30_W_PAD, h_pad=V30_H_PAD, hspace=V30_HSPACE, wspace=V30_WSPACE)
        axes = axes.ravel(); scale = float(np.abs(mot['mean']).max() * 100)
        for i, (region, side) in enumerate((r, s) for r in COMPARTMENTS for s in SIDES):
            ax = axes[i]; fl.style(ax)
            g = mot[mot.contrast.eq(region) & mot.side.eq(side)]; a = g.pivot(index='residue', columns='position', values='mean').sort_index(); qmat = g.pivot(index='residue', columns='position', values='q').reindex_like(a)
            im = ax.imshow(a.values * 100, vmin=-scale, vmax=scale, cmap=fl.div_cmap, aspect='auto')
            ax.set_xticks(range(len(a.columns)), [f'{int(x):+}' for x in a.columns]); ax.set_yticks(range(len(a)), a.index)
            ax.tick_params(axis='y', labelsize=6, length=1.5, pad=1.5)
            ax.set_xlabel('Position relative to boundary')
            ax.set_ylabel(V30_REGION_LABEL.get(region, region.capitalize()) + '; ' + side + ' boundary')
            yy, xx = np.where(qmat.values < Q_CUTOFF); ax.scatter(xx, yy, s=3, c=INK, linewidths=0)
            cb = fig.colorbar(im, ax=ax, fraction=.04, pad=.02, label='Observed − expected (pp)')
            cb.outline.set_linewidth(fl.MIN_LINE_PT); cb.ax.tick_params(labelsize=7, width=fl.MIN_LINE_PT, length=2.6)
            fl.panel(ax, 'abcdef'[i])
        for j, route in enumerate(PREDICT_ROUTES):
            g = pt[pt.route.eq(route) & pt.origin.eq(PREDICT_ORIGIN) & pt.residue.eq(PREDICT_RESIDUE)]; ax = axes[j + 6]; fl.style(ax)
            _v30_dots(ax, [(s, g[g.side.eq(s)].difference) for s in ['N', 'C']], 'H enrichment (percentage points)', ROUTE_COL[route], INK, True)
            for i, side in enumerate(['N', 'C']):
                r = st[st.route.eq(route) & st.origin.eq(PREDICT_ORIGIN) & st.residue.eq(PREDICT_RESIDUE) & st.side.eq(side) & st.endpoint.eq('source_enrichment')].iloc[0]; vals = g[g.side.eq(side)].difference.dropna()
                assert len(vals) == r.n and np.isclose(vals.mean(), r['mean'])
                ax.errorbar(i + .28, r['mean'] * 100, yerr=[[100 * (r['mean'] - r.CI95_low)], [100 * (r.CI95_high - r['mean'])]], fmt='s', c=INK, ms=3, capsize=3, elinewidth=fl.MIN_LINE_PT, capthick=fl.MIN_LINE_PT)
            ax.axhline(0, c=fl.GREY, lw=fl.MIN_LINE_PT)
            h = Line2D([], [], ls='none', marker='o', ms=4, mfc=ROUTE_COL[route], mec='none', alpha=.7)
            ax.legend([h], [V30_ROUTE_LAB[route]], loc='upper left', frameon=False, handletextpad=0.3, borderaxespad=0.2, fontsize=7)
            fl.panel(ax, 'gh'[j])
            yl = ax.get_ylim(); ax.set_yticks([t for t in ref[6 + j]['yticks']]); ax.set_ylim(yl)   # producer tick positions
        fig.canvas.draw()
        rep = v30_compare(ref, v30_fingerprint(fig))
        ok = _v30_finish(fl, fig, out_stem, letters=list(axes))
    assert ok
    return dict(panels=rep, clip_ok=ok)
