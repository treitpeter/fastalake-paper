"""Figure 6: histidine-rich internal peptide termini in intestinal fluid (panels a-f).

Plot code is moved from the original producer chain without change of logic:
``build_figure6.py`` (all six panels on one 180 x 166 mm canvas), then
``postprocess_fig6.py`` (panel-letter replacement, blank-band trim, 183 mm export),
``relabel_regions.py`` (Duodenum..Colon to Region 1-4 as PDF text overlays) and a final
PyMuPDF re-save. Statistics functions are moved from the upstream analyses that wrote
the tables the figure reads: A24 terminal inference (``search.py``, sign-flip tests,
bootstrap intervals, Benjamini-Hochberg) and A28 histidine inference
(``histidine_helpers.signflip``, ``collect_revision4.py``).

source_data/Fig6 files read here (data dictionary in source_data/Fig6/README.md):

* ``capscan_terminal_participant_region.tsv``  panels a, b, A24 terminal-fraction family
* ``capscan_residue_participant_region.tsv``   panel c, A24 residue-composition family
* ``capscan_terminal_statistical_families.tsv`` panels a-c legends (q, n)
* ``capscan_acquisition_design.tsv``            acquisition count (Methods)
* ``capscan_histidine_statistical_families.tsv`` A28 families, source of d-f CapScan rows
* ``capscan_histidine_primary_source_enrichment.tsv``, ``capscan_histidine_primary_reciprocal_tests.tsv``
                                               panels d-f, CapScan rows
* ``capscan_histidine_participant_region_enrichment.tsv``, ``capscan_histidine_participant_reciprocal.tsv``
                                               recomputation of the A28 CapScan families
* ``predict_histidine_primary_source_enrichment.tsv``, ``predict_histidine_primary_reciprocal_tests.tsv``
                                               panels d-f, PREDICT rows
* ``predict_histidine_patient_enrichment.tsv``  PREDICT patient values (pseudonymised)
* ``predict_histidine_statistical_families.tsv`` complete PREDICT A28 families (1,560 + 12 tests)
"""

import hashlib
import io
import itertools
import math
import re
from pathlib import Path

import matplotlib as mpl

mpl.use("Agg")
import matplotlib.font_manager as fm  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

# ---------------------------------------------------------------------------------------
# Source data file names
# ---------------------------------------------------------------------------------------
F_A24_PARTICIPANT = "capscan_terminal_participant_region.tsv"
F_A24_RESIDUE = "capscan_residue_participant_region.tsv"
F_A24_FAMILIES = "capscan_terminal_statistical_families.tsv"
F_A24_DESIGN = "capscan_acquisition_design.tsv"
F_A28_FAMILIES = "capscan_histidine_statistical_families.tsv"
F_A28_PARTICIPANT = "capscan_histidine_participant_region_enrichment.tsv"
F_A28_RECIPROCAL = "capscan_histidine_participant_reciprocal.tsv"
F_CAPSCAN_SOURCE = "capscan_histidine_primary_source_enrichment.tsv"
F_CAPSCAN_RECIPROCAL = "capscan_histidine_primary_reciprocal_tests.tsv"
F_PREDICT_SOURCE = "predict_histidine_primary_source_enrichment.tsv"
F_PREDICT_RECIPROCAL = "predict_histidine_primary_reciprocal_tests.tsv"
F_PREDICT_PATIENT = "predict_histidine_patient_enrichment.tsv"
F_PREDICT_FAMILIES = "predict_histidine_statistical_families.tsv"

# ---------------------------------------------------------------------------------------
# Choices (each a constant of the original code; v28d line of the Methods statement)
# ---------------------------------------------------------------------------------------
MIN_LENGTH = 9                 # terminus analyses use 9-50 residues, v28d lines 285, 855
TERMINUS = "C_last"            # C-terminal residue, v28d lines 169, 855
REGIONS = ["duodenum", "jejunum", "ileum", "colon", "stool"]  # regions 1-4 + stool, line 169
PRIMARY_ORIGINS = ["All", "Human", "Microbial"]  # conservative origin rule, lines 287, 861
SEED_PLOT = 6                  # builder rng: point jitter (a, b) and panel c bootstrap; not stated in v28d
N_BOOT_PANEL_C = 2_000         # Fig. 6c participant-bootstrap resamples, v28d line 869
CI_QUANTILES = (0.025, 0.975)  # 95% percentile intervals (all bootstraps), v28d line 169
# A24 terminal inference (search.py)
A24_SEED_SALT = "20260914"     # seed = first 4 bytes of sha256("|".join(key)+"|20260914"); not stated in v28d
A24_N_BOOT_TERMINAL = 50_000   # non-K/R fraction intervals, v28d line 869
A24_N_BOOT_COMPOSITION = 2_000 # residue-composition intervals, v28d line 869
EXACT_MAX_N = 16               # exact sign enumeration for n <= 16, v28d line 869
N_SIGN_DRAWS = 200_000         # random sign assignments for n > 16 (plus-one), v28d line 869
MIN_PAIRED_N = 8               # fewer paired units: test kept in family with p = 1 (A24/A28 REPORT.md); not stated in v28d
SIGN_TOL = 1e-12               # |null mean| >= |observed mean| - 1e-12 counts as a hit; not stated in v28d
# A28 histidine inference (histidine_helpers.signflip, collect_revision4.ratio_interval)
A28_SEED_NULL = 20260914       # rng(20260914 + n) for n > 16 sign draws; not stated in v28d
A28_SEED_BOOT = 20261014       # rng(20261014 + n) for the mean-difference bootstrap; not stated in v28d
A28_SEED_RATIO = 20261114      # rng(20261114 + n) for the enrichment-ratio bootstrap; not stated in v28d
A28_N_BOOT = 20_000            # internal-boundary and complementary-product intervals, v28d line 869
A28_ROUTES_CAPSCAN = ["score_0_0.5", "score_0.5_0.7", "score_0.7_0.9", "score_0.9_0.95", "score_0.95_0.99",
                      "score_0.99_plus", "denovo_score09", "denovo_score09_eligible", "denovo_top30",
                      "sage_tryptic", "sage_tryptic_bothq", "sage_semi", "sage_semi_bothq", "sage_unspecific",
                      "sage_unspecific_bothq"]  # complete planned grid, family size line 893
A28_PRIMARY = ["denovo_score09", "denovo_top30", "sage_semi", "sage_semi_bothq", "sage_unspecific",
               "sage_unspecific_bothq"]  # complementary-product routes, family size line 893
DISPLAY_ROUTES = ["denovo_top30", "sage_semi_bothq", "sage_unspecific_bothq"]  # d-f routes, lines 285, 867-869
AA = "ACDEFGHIKLMNPQRSTVWY"
# Composite geometry and post-processing (build_figure6.py, postprocess_fig6.py)
W_MM, H_MM = 180, 166          # builder canvas
FINAL_WIDTH_MM = 183           # journal double-column width after export
TRIM_THRESHOLD_PX, TRIM_PAD_PX = 15, 5   # blank-band trim at 72 dpi
RENDER_DPI = 300
REGION_LABELS = {"Duodenum": "Region 1", "Jejunum": "Region 2", "Ileum": "Region 3", "Colon": "Region 4",
                 "duodenum": "region 1", "jejunum": "region 2", "ileum": "region 3", "colon": "region 4"}
#: Panel crop boxes in builder millimetres (x0, y0, x1, y1), used only for per-panel pixel checks.
PANEL_BOXES_MM = {"a": (0, 0, 92, 60), "b": (92, 0, 180, 60), "c": (0, 60, 180, 100),
                  "d": (0, 100, 76, 157), "e": (76, 100, 126, 157), "f": (126, 100, 180, 157)}
NUMBER_WORDS = {14: "fourteen", 15: "fifteen"}
#: Panel b key word for human-origin peptides (v30 legend term, "human-assigned", v30 paragraph 87).
#: The builder legend writes it directly; the v28f producer wrote "Host".
KEY_B_HUMAN = "Human-assigned"
# Former v30 cover-and-retype (v29 figures/_tools/flatten.py --relabel); kept only for the
# optional ``relabel_and_flatten(old=..., new=...)`` path, not used by the notebook.
KEY_B_OLD = "Host: capsule > stool in 15/15, q < 0.001"
KEY_B_NEW = "Human-assigned: capsule > stool in 15/15, q < 0.001"
FLATTEN_DPI = 600              # v30 PDF holds one 600 dpi raster only (covered text removed)
COVER_PAD_PT = 0.3             # white cover box padding around the old span, points

# Palette and ink exactly as resolved by the producer's style module (Okabe-Ito).
INK, GREY = "#4e4f51", "#7f7f7f"
FL, METAG, GREEN, HOST, SKY = "#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9"


# =======================================================================================
# Style: the producer's settings, reproduced exactly
# =======================================================================================
def producer_style(font_dir):
    """Set matplotlib rcParams exactly as the original producer had them.

    The producer imported the lab styleguide ``init_plotting`` (registers ``soehne-buch.otf``
    as family ``soehne``, ink #4e4f51, 0.6 pt lines), then the project ``fl_style.init``
    of that date (no extra weights, 0.6 pt lines), then the builder's own overrides. The
    repository ``style/fl_style.py`` has since moved to 1 pt lines, so it is not used for
    this figure; this function restores the producer's values instead.

    Parameters
    ----------
    font_dir : str or pathlib.Path or None
        Directory with ``soehne-buch.otf``. Without it text falls back to the
        matplotlib default face.

    Returns
    -------
    bool
        True if the Söhne face was registered.
    """
    mpl.rcdefaults()
    mpl.rcParams["figure.dpi"] = 300
    mpl.rcParams["pdf.fonttype"] = 42
    mpl.rcParams["ps.fonttype"] = 42
    mpl.rcParams["svg.fonttype"] = "none"
    font = Path(font_dir) / "soehne-buch.otf" if font_dir else None
    have = bool(font) and font.is_file()
    if have:
        fe = fm.FontEntry(fname=str(font), name="soehne")
        fm.fontManager.ttflist.insert(0, fe)
        mpl.rcParams["font.family"] = fe.name
    else:
        mpl.rcParams["font.family"] = "sans-serif"
    for k in ("text.color", "axes.labelcolor", "axes.edgecolor", "xtick.color", "ytick.color"):
        mpl.rcParams[k] = INK
    for k in ("lines.linewidth", "axes.linewidth", "xtick.major.width", "ytick.major.width",
              "xtick.minor.width", "ytick.minor.width"):
        mpl.rcParams[k] = 0.6
    mpl.rcParams.update({"font.size": 7, "axes.titleweight": "normal", "axes.titlesize": 7,
                         "legend.fontsize": 6.5, "xtick.major.size": 2.6, "ytick.major.size": 2.6,
                         "figure.dpi": 300, "savefig.dpi": 300, "savefig.bbox": "tight",
                         "savefig.pad_inches": 0.02, "svg.fonttype": "none"})
    mpl.rcParams.update({"mathtext.fontset": "custom", "mathtext.rm": "soehne", "mathtext.it": "soehne",
                         "mathtext.bf": "soehne", "mathtext.sf": "soehne", "mathtext.tt": "soehne",
                         "mathtext.cal": "soehne", "mathtext.default": "regular"})
    # build_figure6.py overrides
    plt.rcParams["savefig.bbox"] = None
    plt.rcParams.update({"axes.labelcolor": INK, "text.color": INK, "xtick.color": INK, "ytick.color": INK,
                         "font.size": 6.5, "axes.labelsize": 6.5, "xtick.labelsize": 6.5,
                         "ytick.labelsize": 6.5, "legend.fontsize": 6})
    return have


def _style(ax):
    """Producer ``fl_style.style``: no top/right spines, 2.6 pt ticks 0.6 pt wide, 7 pt."""
    for sp in ("top", "right"):
        ax.spines[sp].set_visible(False)
    ax.tick_params(length=2.6, width=0.6, labelsize=7)
    ax.xaxis.label.set_size(7)
    ax.yaxis.label.set_size(7)
    return ax


def _panel(ax, letter):
    """Producer ``panel()`` in build_figure6.py: style the axes and draw the letter."""
    _style(ax)
    ax.set_title(letter, loc="left", fontsize=8, color=INK, pad=4)


# =======================================================================================
# Loading
# =======================================================================================
def load(sd):
    """Read every Figure 6 table.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/Fig6``.

    Returns
    -------
    dict of pandas.DataFrame
        Keys ``pp`` (A24 participant-region values; participant ids read as integers, as
        in the figure builder), ``fam`` (A24 families), ``design``, ``a28``, ``pe``, ``pr``, ``S`` and
        ``R`` (panel d-f source and reciprocal rows, CapScan + PREDICT, origin All),
        ``predict_patient``.
    """
    rd = lambda f, **k: pd.read_csv(sd / f, sep="\t", **k)  # noqa: E731
    # The builder read participant ids as integers (pandas default); the order of points
    # within a region, and so the jitter and the panel c bootstrap draws, follow that order.
    pp = pd.concat([rd(F_A24_PARTICIPANT), rd(F_A24_RESIDUE)], ignore_index=True)
    d = dict(pp=pp, fam=rd(F_A24_FAMILIES),
             design=rd(F_A24_DESIGN, dtype={"participant": str}), a28=rd(F_A28_FAMILIES),
             pe=rd(F_A28_PARTICIPANT, keep_default_na=False, na_values=[""]),
             pr=rd(F_A28_RECIPROCAL, keep_default_na=False, na_values=[""]),
             predict_patient=rd(F_PREDICT_PATIENT), predict_fam=rd(F_PREDICT_FAMILIES))
    S = pd.concat([rd(F_CAPSCAN_SOURCE).assign(cohort="CAPSCAN"), rd(F_PREDICT_SOURCE).assign(cohort="PREDICT")])
    R = pd.concat([rd(F_CAPSCAN_RECIPROCAL).assign(cohort="CAPSCAN"),
                   rd(F_PREDICT_RECIPROCAL).assign(cohort="PREDICT")])
    d["S"], d["R"] = S[S.origin.eq("All")], R[R.origin.eq("All")]
    return d


def derive_histidine_inputs(a28):
    """Select the panel d-f CapScan rows from the A28 family table.

    Moved from ``prep_histidine_inputs.py`` (same filter lines as the earlier
    ``review_histidine_capscan.py``).

    Parameters
    ----------
    a28 : pandas.DataFrame
        A28 CapScan ``complete_statistical_families`` table.

    Returns
    -------
    tuple of pandas.DataFrame
        (source-enrichment rows for residue H, matched-reciprocal rows), display routes only.
    """
    stats = a28
    routes = DISPLAY_ROUTES
    primary = stats[(stats.endpoint == "source_enrichment") & stats.residue.eq("H") & stats.route.isin(routes)].copy()
    rec = stats[(stats.endpoint == "matched_reciprocal") & stats.route.isin(routes)]
    return primary, rec


# =======================================================================================
# Panels (build_figure6.py, unchanged logic)
# =======================================================================================
def new_figure():
    """Open the 180 x 166 mm builder canvas.

    Returns
    -------
    tuple
        ``(fig, ax_mm, state)`` where ``ax_mm(x, y, w, h)`` adds axes in millimetres from
        the top-left corner and ``state`` carries the shared random generator (seed
        ``SEED_PLOT``) and the dictionary of plotted numbers.
    """
    mmw = 1 / 25.4
    fig = plt.figure(figsize=(W_MM * mmw, H_MM * mmw))

    def ax_mm(x, y, w, h):
        """Add axes at ``x, y`` mm from the top-left corner, ``w`` x ``h`` mm in size.

        Returns
        -------
        matplotlib.axes.Axes
        """
        return fig.add_axes([x / W_MM, 1 - (y + h) / H_MM, w / W_MM, h / H_MM])

    return fig, ax_mm, dict(rng=np.random.default_rng(SEED_PLOT), nums={})


def _vals(pp, route, unit, origin, metric):
    """Participant x region table of one endpoint at ``MIN_LENGTH`` and ``TERMINUS``.

    Parameters
    ----------
    pp : pandas.DataFrame
        ``capscan_terminal_participant_region.tsv``.
    route, unit, origin, metric : str
        Values of the columns of the same name.

    Returns
    -------
    pandas.DataFrame
        Index participant, columns region, values ``value`` (fractions, not percent).
    """
    d = pp[pp.route.eq(route) & pp.min_length.eq(MIN_LENGTH) & pp.unit.eq(unit) & pp.origin.eq(origin)
           & pp.terminus.eq(TERMINUS) & pp.metric.eq(metric)]
    return d.pivot(index="participant", columns="region", values="value")


def _qrow(fam, route, unit, origin, metric, contrast, family):
    """The single A24 test row for one endpoint, contrast and family (asserts exactly one).

    Parameters
    ----------
    fam : pandas.DataFrame
        A24 test table (``capscan_terminal_statistical_families.tsv`` or its recomputation).
    route, unit, origin, metric, contrast, family : str
        Values of the columns of the same name; ``min_length`` and ``terminus`` are fixed
        to ``MIN_LENGTH`` and ``TERMINUS``.

    Returns
    -------
    pandas.Series
    """
    d = fam[fam.route.eq(route) & fam.min_length.eq(MIN_LENGTH) & fam.unit.eq(unit) & fam.origin.eq(origin)
            & fam.terminus.eq(TERMINUS) & fam.metric.eq(metric) & fam.contrast.eq(contrast) & fam.family.eq(family)]
    assert len(d) == 1, (route, unit, origin, metric, contrast, len(d))
    return d.iloc[0]


def _qtxt(q):
    """Legend text of a Benjamini-Hochberg q value: "q < 0.001" below 0.001, else three decimals.

    Returns
    -------
    str
    """
    return "q < 0.001" if q < 0.001 else f"q = {q:.3f}"


REGL = ["Duodenum", "Jejunum", "Ileum", "Colon", "Stool"]
ROUTES = [("denovo", "De novo (top 30%)", FL, "o"), ("semi", "Sage semi-tryptic", METAG, "s"),
          ("unspecific", "Sage non-specific", GREEN, "^")]


def panel_a(fig, ax_mm, state, pp, fam):
    """Panel a: non-K/R C-terminal peptide fraction by region and route.

    Parameters
    ----------
    fig, ax_mm, state
        From :func:`new_figure`.
    pp, fam : pandas.DataFrame
        A24 participant-region values and A24 families.

    Returns
    -------
    matplotlib.axes.Axes
        Points are participant-region means (unit participant, n = 15 per region); the
        legend reports the paired sign-flip test of capsule mean minus stool, BH-adjusted
        within ``terminal_fraction__length_9``.
    """
    rng, nums = state["rng"], state["nums"]
    axA = ax_mm(15, 5, 72, 34); _panel(axA, "a"); off = [-.26, 0, .26]; nums["a"] = {}
    tw = _vals(pp, "tryptic", "canonical_peptide", "All", "nonKR_fraction_all") * 100
    for (route, lab, col, mk), o in zip(ROUTES, off):
        w = _vals(pp, route, "canonical_peptide", "All", "nonKR_fraction_all") * 100; nums["a"][route] = {}
        for i, r in enumerate(REGIONS):
            v = w[r].dropna(); axA.scatter(i + o + rng.uniform(-.05, .05, len(v)), v, s=3, color=col, alpha=.5, lw=0, zorder=2)
            axA.plot([i + o - .1, i + o + .1], [v.mean()] * 2, color="k", lw=1, zorder=3); nums["a"][route][r] = float(v.mean())
            nums["a"][route][r + "_n"] = int(len(v))
        q = _qrow(fam, route, "canonical_peptide", "All", "nonKR_fraction_all", "capsule__minus__stool", "terminal_fraction__length_9")
        nums["a"][route]["capsule_minus_stool"] = dict(mean_pp=100 * q.mean_difference, pos=int(q.n_positive), neg=int(q.n_negative), n=int(q.n_participants), p=float(q.p), q=float(q.q))
        axA.scatter([], [], s=8, color=col, marker="o", label=f"{lab}: capsule > stool in {int(q.n_positive)}/{int(q.n_participants)}, {_qtxt(q.q)}")
    axA.axhline(float(np.nanmean(tw.values)), color=GREY, lw=.6, ls="--"); axA.text(4.55, float(np.nanmean(tw.values)) + .3, "tryptic search", fontsize=5.5, color=GREY, ha="right", va="bottom")
    axA.set_xticks(range(5), REGL); axA.tick_params(axis="x", length=0); axA.set_ylabel("Peptides with a non-K/R\nC-terminus (%)"); axA.set_ylim(0, 14); axA.set_yticks([0, 4, 8, 12]); axA.set_xlim(-.55, 4.6)
    axA.legend(frameon=False, loc="upper left", bbox_to_anchor=(-.02, -.13), ncol=1, handletextpad=.2, labelspacing=.2, borderpad=0)
    nums["a"]["tryptic_mean_pct"] = float(np.nanmean(tw.values))
    return axA


def panel_b(fig, ax_mm, state, pp, fam):
    """Panel b: histidine share of non-K/R C termini, host versus microbial, non-specific route.

    Returns
    -------
    matplotlib.axes.Axes
        Bars are means of participant-region values (unit participant); legend tests as
        in panel a, BH within ``terminal_fraction__length_9``.
    """
    rng, nums = state["rng"], state["nums"]
    axB = ax_mm(106, 5, 70, 34); _panel(axB, "b"); nums["b"] = {}; hB = []
    for origin, col, o in (("Human", HOST, -.18), ("Microbial", FL, .18)):
        w = _vals(pp, "unspecific", "canonical_peptide", origin, "H_fraction_nonKR") * 100; nums["b"][origin] = {}
        for i, r in enumerate(REGIONS):
            v = w[r].dropna(); axB.bar(i + o, v.mean(), width=.34, color=col, alpha=.9, zorder=1); axB.scatter(i + o + rng.uniform(-.07, .07, len(v)), v, s=3, facecolor="white", edgecolor=INK, lw=.3, zorder=2); nums["b"][origin][r] = float(v.mean())
        q = _qrow(fam, "unspecific", "canonical_peptide", origin, "H_fraction_nonKR", "capsule__minus__stool", "terminal_fraction__length_9")
        nums["b"][origin]["capsule_minus_stool"] = dict(mean_pp=100 * q.mean_difference, pos=int(q.n_positive), neg=int(q.n_negative), n=int(q.n_participants), p=float(q.p), q=float(q.q), ci_pp=(100 * q.CI95_low, 100 * q.CI95_high))
        hB.append(plt.Rectangle((0, 0), 1, 1, color=col, label=f'{KEY_B_HUMAN if origin == "Human" else "Microbial"}: capsule > stool in {int(q.n_positive)}/{int(q.n_participants)}, {_qtxt(q.q)}'))
    axB.set_xticks(range(5), REGL); axB.tick_params(axis="x", length=0); axB.set_ylabel("Histidine among non-K/R\nC-termini (%)"); axB.set_xlim(-.55, 4.55)
    axB.legend(handles=hB, frameon=False, loc="upper left", bbox_to_anchor=(-.02, -.13), handletextpad=.4, labelspacing=.2, borderpad=0, handlelength=1)
    return axB


def panel_c(fig, ax_mm, state, pp, fam):
    """Panel c: residue composition of non-K/R C termini (non-specific route).

    Bars are means over participants of the region-group value (regions 1-2 = mean of
    duodenum and jejunum, regions 3-4 = mean of ileum and colon, both required); error
    bars are 95% percentile intervals of ``N_BOOT_PANEL_C`` participant bootstrap
    resamples drawn from the shared generator (seed ``SEED_PLOT``).

    Returns
    -------
    matplotlib.axes.Axes
    """
    rng, nums = state["rng"], state["nums"]
    AAc = [a for a in "ACDEFGHILMNPQSTVWY"]; axC = ax_mm(15, 64, 161, 28); _panel(axC, "c"); nums["c"] = {}
    comp = {}
    for a in AAc:
        w = _vals(pp, "unspecific", "assigned_bare_sequence", "All", f"residue_{a}_fraction_nonKR") * 100
        prox = w[["duodenum", "jejunum"]].mean(axis=1, skipna=False); dist = w[["ileum", "colon"]].mean(axis=1, skipna=False); comp[a] = {"proximal": prox.dropna(), "distal": dist.dropna(), "stool": w["stool"].dropna()}
    order = sorted(AAc, key=lambda a: -(comp[a]["proximal"].mean() + comp[a]["distal"].mean()) / 2)
    cols = {"proximal": SKY, "distal": FL, "stool": METAG}; offs = {"proximal": -.27, "distal": 0, "stool": .27}
    labs = {"proximal": "Proximal fluid (duodenum + jejunum)", "distal": "Distal fluid (ileum + colon)", "stool": "Stool"}
    iH = order.index("H"); axC.axvspan(iH - .45, iH + .45, color="#f2f2f2", zorder=0)
    nums["c"]["ci"] = {}
    for k in ["proximal", "distal", "stool"]:
        m = [comp[a][k].mean() for a in order]; lo = []; hi = []
        for a in order:
            v = comp[a][k].values; bs = rng.choice(v, (N_BOOT_PANEL_C, len(v))).mean(1); nums["c"]["draws"] = len(bs); lo.append(np.quantile(bs, CI_QUANTILES[0])); hi.append(np.quantile(bs, CI_QUANTILES[1]))
        x = np.arange(len(order)) + offs[k]; axC.bar(x, m, width=.26, color=cols[k], label=labs[k], zorder=1)
        axC.errorbar(x, m, yerr=[np.array(m) - np.array(lo), np.array(hi) - np.array(m)], fmt="none", ecolor=INK, elinewidth=.5, capsize=0, zorder=2)
        nums["c"][k] = {a: float(comp[a][k].mean()) for a in order}; nums["c"][k + "_n"] = int(len(comp["H"][k]))
        nums["c"]["ci"][k] = {a: (float(l_), float(h_)) for a, l_, h_ in zip(order, lo, hi)}
    nums["c"]["order"] = order
    axC.set_xticks(range(len(order)), order); axC.tick_params(axis="x", length=0); axC.set_ylabel("Share of non-K/R\nC-termini (%)"); axC.set_xlim(-.6, len(order) - .4)
    axC.legend(frameon=False, loc="upper right", ncol=3, handlelength=1, columnspacing=1)
    for k in ("proximal", "distal"):
        q = _qrow(fam, "unspecific", "assigned_bare_sequence", "All", "residue_H_fraction_nonKR", f"{k}__minus__stool", "all18_nonKR_composition__length_9"); nums["c"][f"H_{k}_minus_stool"] = dict(mean_pp=100 * q.mean_difference, pos=int(q.n_positive), neg=int(q.n_negative), p=float(q.p), q=float(q.q))
    axC.text(iH + .5, 43, "Histidine", ha="left", va="center", fontsize=6)
    return axC


FR = [("denovo_top30", "De novo (top 30%)", FL, "o"), ("sage_semi_bothq", "Sage semi-tryptic", METAG, "s"),
      ("sage_unspecific_bothq", "Sage non-specific", GREEN, "^")]
ROWS = [("CAPSCAN", "duodenum", "Duodenum"), ("CAPSCAN", "jejunum", "Jejunum"), ("CAPSCAN", "ileum", "Ileum"),
        ("CAPSCAN", "colon", "Colon"), ("CAPSCAN", "stool", "Stool"), ("PREDICT", "ALL_observed_specimens", "PREDICT stool")]
_DEF = [("C", "d", 13, "Histidine at C-boundaries,\nobserved − expected (points)"),
        ("N", "e", 70, "Histidine at N-boundaries,\nobserved − expected (points)"),
        ("rec", "f", 127, "Complementary products,\nhistidine − control (points)")]


def _boundary_panel(k, fig, ax_mm, state, S, R):
    """Draw one of panels d, e, f: route means with 95% bootstrap intervals per cohort row.

    Each symbol is the mean over participants (CapScan, n = 15) or patients (PREDICT,
    n = 93) of the per-unit value, with the A28 percentile bootstrap interval
    (``A28_N_BOOT`` resamples of whole participants or patients). Values are read from
    the A28 tables, not recomputed here.

    Parameters
    ----------
    k : int
        0 = panel d (C-boundaries), 1 = panel e (N-boundaries), 2 = panel f
        (complementary products); selects the entry of ``_DEF``.
    fig, ax_mm, state
        From :func:`new_figure`; plotted values are stored in ``state["nums"]["def"]``.
    S : pandas.DataFrame
        Internal-boundary observed-minus-expected summaries (``side`` C or N).
    R : pandas.DataFrame
        Complementary-product summaries (used when ``k == 2``).

    Returns
    -------
    matplotlib.axes.Axes
    """
    nums = state["nums"]; nums.setdefault("def", {})
    key, lab, x0, xl = _DEF[k]
    ax = ax_mm([30, 82, 132][k], 104, [42, 42, 44][k], 42); _panel(ax, lab); d = S[S.side.eq(key)] if key != "rec" else R
    ax.axvline(0, color=GREY, lw=.6); ax.axhline(4.5, color="#bdbdbd", lw=.5, ls="--"); nums["def"][key] = {}
    for ri, (route, rl, col, mk) in enumerate(FR):
        for j, (co, ct, _) in enumerate(ROWS):
            x = d[d.route.eq(route) & d.contrast.eq(ct) & d.cohort.eq(co)]
            if x.empty:
                assert route == "sage_semi_bothq" and co == "PREDICT"; continue
            assert len(x) == 1; v = x.iloc[0]; y = j + (ri - 1) * .24
            ax.plot([100 * v.CI95_low, 100 * v.CI95_high], [y, y], color=col, lw=.9); ax.scatter([100 * v["mean"]], [y], s=9, color=col, marker=mk, lw=0, zorder=3)
            nums["def"][key][f"{route}:{co}:{ct}"] = [100 * v["mean"], 100 * v.CI95_low, 100 * v.CI95_high, float(v.q), int(v.n)]
    ax.set_ylim(5.6, -.6); ax.set_yticks(range(6), [r[2] for r in ROWS] if k == 0 else [""] * 6); ax.tick_params(axis="y", length=0); ax.set_xlabel(xl)
    lo, hi = ax.get_xlim(); ax.set_xlim(min(lo, -.5), hi)
    return ax


def panel_d(fig, ax_mm, state, S, R):
    """Panel d: histidine at internal C-terminal boundaries, observed minus source-expected.

    Symbols are participant (CapScan, n = 15) or patient (PREDICT, n = 93) means with 95%
    bootstrap intervals (``A28_N_BOOT`` resamples of whole participants or patients).

    Returns
    -------
    matplotlib.axes.Axes
    """
    return _boundary_panel(0, fig, ax_mm, state, S, R)


def panel_e(fig, ax_mm, state, S, R):
    """Panel e: as panel d for internal N-terminal boundaries.

    Returns
    -------
    matplotlib.axes.Axes
    """
    return _boundary_panel(1, fig, ax_mm, state, S, R)


def panel_f(fig, ax_mm, state, S, R):
    """Panel f: complementary-product rate at histidine minus matched non-histidine boundaries.

    Also draws the shared route legend at the bottom of the figure (drawn last in the
    producer).

    Returns
    -------
    matplotlib.axes.Axes
    """
    ax = _boundary_panel(2, fig, ax_mm, state, S, R)
    for route, rl, col, mk in FR:
        fig.axes[-1].scatter([], [], s=12, color=col, marker=mk, label=rl)
    fig.legend(*fig.axes[-1].get_legend_handles_labels(), frameon=False, loc="lower center", ncol=3, bbox_to_anchor=(.55, -.005), handletextpad=.2, columnspacing=1.2)
    return ax


def save_builder(fig, out):
    """Save the builder composite as ``panels/Figure6.pdf`` and ``panels/Figure6.png`` (300 dpi).

    Returns
    -------
    tuple of pathlib.Path
    """
    p = Path(out) / "panels"; p.mkdir(parents=True, exist_ok=True)
    fig.savefig(p / "Figure6.pdf"); fig.savefig(p / "Figure6.png", dpi=300)
    return p / "Figure6.pdf", p / "Figure6.png"


def panel_crop(fig, letter, out):
    """Render the figure as it stands at 300 dpi and save the crop of one panel.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
    letter : str
        Key of ``PANEL_BOXES_MM``.
    out : pathlib.Path
        Writes ``out / f"panel_{letter}.png"``.

    Returns
    -------
    PIL.Image.Image
    """
    from PIL import Image
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=300); buf.seek(0)
    im = crop_image(Image.open(buf), letter)
    im.save(Path(out) / f"panel_{letter}.png")
    return im


def crop_image(im, letter):
    """Crop a 300 dpi builder-canvas image to the panel box ``PANEL_BOXES_MM[letter]``.

    Parameters
    ----------
    im : PIL.Image.Image
        Builder canvas rendered at 300 dpi.
    letter : str
        Panel letter, a key of ``PANEL_BOXES_MM``.

    Returns
    -------
    PIL.Image.Image
        The crop; the box is clipped to the image edge.
    """
    x0, y0, x1, y1 = (round(v / 25.4 * 300) for v in PANEL_BOXES_MM[letter])
    return im.crop((x0, y0, min(x1, im.width), min(y1, im.height)))


# =======================================================================================
# Post-processing (postprocess_fig6.py, relabel_regions.py, final re-save)
# =======================================================================================
def _trim_blank_bands(doc, threshold=TRIM_THRESHOLD_PX, pad=TRIM_PAD_PX):
    """Remove horizontal all-white bands taller than ``threshold`` px (72 dpi) from page 1.

    A row is white when every pixel channel is above 249. ``pad`` px are kept at each
    side of a removed band. Producer step of ``postprocess_fig6.py``.

    Parameters
    ----------
    doc : pymupdf.Document
    threshold, pad : int

    Returns
    -------
    tuple
        ``(document, cuts)``; ``cuts`` lists the removed ``(start, end)`` bands in points.
        The input document is returned unchanged when no band qualifies.
    """
    import pymupdf as f
    p = doc[0]; pix = p.get_pixmap(dpi=72, alpha=False); a = np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, 3)
    white = np.all(a > 249, axis=(1, 2)); cuts = []; start = None
    for y, v in enumerate(list(white) + [False]):
        if v and start is None: start = y
        elif not v and start is not None:
            if y - start > threshold: cuts.append((start + pad, y - pad))
            start = None
    if not cuts: return doc, []
    out = f.open(); dest = out.new_page(width=p.rect.width, height=p.rect.height - sum(b - a for a, b in cuts)); srcy = 0; dsty = 0
    for a, b in cuts + [(p.rect.height, p.rect.height)]:
        if a > srcy: rect = f.Rect(0, srcy, p.rect.width, a); dest.show_pdf_page(f.Rect(0, dsty, p.rect.width, dsty + rect.height), doc, 0, clip=rect); dsty += rect.height
        srcy = b
    return out, cuts


def _labels_and_small_fixes(doc, font_dir):
    """Rewrite single-letter panel labels (a-i, size >= 5 pt) in Söhne Halbfett 8 pt.

    The original spans are redacted (text only, no images or vector graphics) and the
    letters re-inserted at the same origin. Producer step of ``postprocess_fig6.py``.

    Parameters
    ----------
    doc : pymupdf.Document
    font_dir : str or pathlib.Path
        Must contain ``soehne-halbfett.ttf`` and ``soehne-buch.ttf``.

    Returns
    -------
    pymupdf.Document
        The same document, modified in place.
    """
    import pymupdf as f
    p = doc[0]; ss = [s for b in p.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"]]
    changes = [s for s in ss if re.fullmatch("[a-i]", s["text"]) and s["size"] >= 5]
    for s in changes: p.add_redact_annot(f.Rect(s["bbox"]), fill=False, cross_out=False)
    if changes: p.apply_redactions(images=0, graphics=0)
    for key, file in [("flbold", "soehne-halbfett.ttf"), ("flreg", "soehne-buch.ttf")]: p.insert_font(fontname=key, fontfile=str(Path(font_dir) / file))
    for s in changes:
        p.insert_text(list(s["origin"]), s["text"], fontsize=8, fontname="flbold", color=(.15, .15, .15))
    return doc


def postprocess(builder_pdf, outdir, font_dir):
    """Panel letters to Söhne Halbfett 8 pt, trim blank bands, export at 183 mm width, 300 dpi PNG.

    Moved from ``postprocess_fig6.py`` (k = 6 branch of the original r12 loop).

    Parameters
    ----------
    builder_pdf : pathlib.Path
    outdir : pathlib.Path
    font_dir : pathlib.Path
        Must contain ``soehne-halbfett.ttf`` and ``soehne-buch.ttf``.

    Returns
    -------
    dict
        Paths written, blank bands removed, page size and minimum font size.
    """
    import fitz
    import pymupdf as f
    O = Path(outdir); O.mkdir(parents=True, exist_ok=True)
    d = f.open(builder_pdf); d = _labels_and_small_fixes(d, font_dir); d, cuts = _trim_blank_bands(d)
    p = O / "Figure6_r12.pdf"
    if p.exists(): p.unlink()
    d.save(p, garbage=4, deflate=True)
    doc = fitz.open(p); pg = doc[0]
    nat = O / "Figure6_native.pdf"
    if nat.exists(): nat.unlink()
    doc.save(nat, garbage=4, deflate=True)
    scaled = fitz.open(); width = FINAL_WIDTH_MM / 25.4 * 72; factor = width / pg.rect.width
    page = scaled.new_page(width=width, height=pg.rect.height * factor); page.show_pdf_page(page.rect, doc, 0)
    dest = O / "Figure6.pdf"; scaled.save(dest, garbage=4, deflate=True)
    final = fitz.open(dest); page = final[0]
    page.get_pixmap(dpi=300, alpha=False).save(O / "Figure6.png")
    spans = [s for b in page.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"] if s["text"].strip()]
    smallest = min(s["size"] for s in spans); assert smallest >= 5, smallest
    return dict(pdf=dest, png=O / "Figure6.png", blank_bands_removed=cuts, page_pt=list(page.rect)[2:], min_font_pt=smallest)


def relabel_regions(src_pdf, outdir, font_dir):
    """Replace region names by Region 1-4 as opaque PDF text overlays (``relabel_regions.py``).

    The original vector text stays underneath; every original span is asserted unchanged.

    Returns
    -------
    dict
        ``pdf`` and ``png`` written to ``outdir/region_figures`` and the list of edits.
    """
    import fitz
    H = Path(outdir); out = H / "region_figures"; out.mkdir(parents=True, exist_ok=True)
    font = Path(font_dir) / "soehne-buch.ttf"
    d = fitz.open(src_pdf); edits = []
    for page in d:
        spans = []; before = page.get_text()
        for block in page.get_text("rawdict")["blocks"]:
            for line in block.get("lines", []):
                for sp in line["spans"]:
                    sp["text"] = "".join(c["c"] for c in sp["chars"]); old = sp["text"]; new = old
                    for a, b in REGION_LABELS.items(): new = re.sub(r"\b" + a + r"\b", b, new)
                    if new != old:
                        new = new.replace("Proximal fluid (region 1 + region 2)", "Fluid regions 1–2").replace("Distal fluid (region 3 + region 4)", "Fluid regions 3–4")
                        spans.append((sp, line["dir"], new))
                        quad = fitz.recover_quad(line["dir"], sp)
                        shape = page.new_shape(); shape.draw_quad(quad); shape.finish(color=None, fill=(1, 1, 1)); shape.commit()
        page.insert_font(fontname="KCRegion", fontfile=str(font))
        for sp, direction, new in spans:
            origin = fitz.Point(sp["origin"]); angle = math.degrees(math.atan2(-direction[1], direction[0]))
            color = tuple(((sp["color"] >> shift) & 255) / 255 for shift in [16, 8, 0])
            page.insert_text(origin, new, fontname="KCRegion", fontsize=sp["size"], color=color, morph=(origin, fitz.Matrix(angle)) if abs(angle) > 1e-3 else None)
            edits.append(dict(before=sp["text"], after=new, font_size=sp["size"]))
        after = page.get_text()
        assert after.startswith(before), "Original vector text changed unexpectedly"
    target = out / "Figure6.pdf"
    if target.exists(): target.unlink()
    d.save(target, garbage=4, deflate=True)
    nd = fitz.open(target); nd[0].get_pixmap(dpi=300, alpha=False).save(out / "Figure6.png")
    return dict(pdf=target, png=out / "Figure6.png", edits=edits)


def final_resave(region_pdf, region_png, out):
    """Final step of ``run_all.sh``: PyMuPDF re-save of the PDF (garbage=4, deflate) and PNG copy.

    Returns
    -------
    tuple of pathlib.Path
        ``(out/Figure6.pdf, out/Figure6.png)``, the v28f composite.
    """
    import shutil
    import fitz
    d = fitz.open(region_pdf); d.save(Path(out) / "Figure6.pdf", garbage=4, deflate=True)
    shutil.copyfile(region_png, Path(out) / "Figure6.png")
    return Path(out) / "Figure6.pdf", Path(out) / "Figure6.png"


def relabel_and_flatten(v28f_pdf, out, font_dir, old=None, new=None, dpi=FLATTEN_DPI):
    """v30 step: flatten the composite to one raster (optionally relabel a text span first).

    The panel b key is written as "Human-assigned" by :func:`panel_b`, so the notebook calls
    this with ``old=None`` (no relabel). With ``old`` and ``new`` given it reproduces the
    former cover-and-retype of the v29 tool.

    Moved from the v29 figure tool ``flatten.py --relabel`` without change of logic. The
    visible span ``old`` is covered with its background colour (sampled just above the span
    at 600 dpi) and ``new`` is typed in the span's Söhne face, size, colour and origin.
    The page is then written as a PDF holding only one ``dpi`` raster (no text, no vectors),
    and a 300 dpi PNG is rendered from the relabelled vector page.

    Parameters
    ----------
    v28f_pdf : path-like
        The v28f composite (output of :func:`final_resave`).
    out : path-like
        Output directory; writes ``Figure6_v30.pdf`` and ``Figure6_v30.png``.
    font_dir : path-like
        Directory with ``soehne-buch.ttf``, ``soehne-halbfett.ttf`` and ``soehne-kraftig.ttf``.
    old, new : str
        Key text before and after.
    dpi : int
        Raster resolution of the flattened PDF.

    Returns
    -------
    dict
        ``pdf``, ``png`` paths and a receipt (old span, cover colour, widths, changed rectangle).
    """
    import fitz
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    src = fitz.open(v28f_pdf); p = src[0]; W, H = p.rect.width, p.rect.height
    rec = dict(page_pt=[W, H])
    work = fitz.open(); work.insert_pdf(src); wp = work[0]
    if old is None:
        pm6 = wp.get_pixmap(dpi=dpi, alpha=False)
        rec.update(relabel=None)
        o = fitz.open(); op = o.new_page(width=W, height=H); op.insert_image(op.rect, pixmap=pm6)
        pdf = out / "Figure6_v30.pdf"; png = out / "Figure6_v30.png"
        if pdf.exists(): pdf.unlink()
        o.save(pdf, garbage=4, deflate=True)
        chk = fitz.open(pdf)[0]
        assert chk.get_text().strip() == "" and len(chk.get_drawings()) == 0 and len(chk.get_images()) == 1
        pn = wp.get_pixmap(dpi=RENDER_DPI, alpha=False); pn.set_dpi(RENDER_DPI, RENDER_DPI); pn.save(png)
        rec["flat_pdf_raster_px"] = [pm6.width, pm6.height]
        return dict(pdf=pdf, png=png, receipt=rec)
    sp = [s for b in p.get_text("dict")["blocks"] for l in b.get("lines", []) for s in l["spans"] if s["text"].strip() == old]
    assert len(sp) >= 1, old
    rec["old_spans"] = [dict(bbox=list(s["bbox"]), font=s["font"], size=s["size"], color="%06x" % s["color"]) for s in sp]
    s = sp[-1]                                        # last drawn = visible copy
    r = fitz.Rect(s["bbox"])
    pm = p.get_pixmap(dpi=600, clip=fitz.Rect(r.x0, r.y0 - 1, r.x1 + 1, r.y0 - 0.3))
    bgc = tuple(v / 255 for v in pm.pixel(pm.width // 2, 0))
    face = "soehne-halbfett.ttf" if "Halb" in s["font"] else ("soehne-kraftig.ttf" if "Kraft" in s["font"] else "soehne-buch.ttf")
    font = fitz.Font(fontfile=str(Path(font_dir) / face))
    nw = font.text_length(new, s["size"])
    wp.draw_rect(fitz.Rect(r.x0 - COVER_PAD_PT, r.y0, max(r.x1, s["origin"][0] + nw) + COVER_PAD_PT, r.y1), color=None, fill=bgc, width=0)
    tw = fitz.TextWriter(wp.rect); tw.append(s["origin"], new, font=font, fontsize=s["size"])
    tw.write_text(wp, color=tuple(((s["color"] >> k) & 255) / 255 for k in (16, 8, 0)))
    rec.update(relabel=[old, new], face=face, cover_bg=[round(v * 255) for v in bgc], new_width_pt=nw, old_width_pt=r.width,
               changed_rect_pt=[r.x0 - COVER_PAD_PT, r.y0, max(r.x1, s["origin"][0] + nw) + COVER_PAD_PT, r.y1])
    pm6 = wp.get_pixmap(dpi=dpi, alpha=False)
    o = fitz.open(); op = o.new_page(width=W, height=H); op.insert_image(op.rect, pixmap=pm6)
    pdf = out / "Figure6_v30.pdf"; png = out / "Figure6_v30.png"
    if pdf.exists(): pdf.unlink()
    o.save(pdf, garbage=4, deflate=True)
    chk = fitz.open(pdf)[0]
    assert chk.get_text().strip() == "" and len(chk.get_drawings()) == 0 and len(chk.get_images()) == 1
    pn = wp.get_pixmap(dpi=RENDER_DPI, alpha=False); pn.set_dpi(RENDER_DPI, RENDER_DPI); pn.save(png)
    rec["flat_pdf_raster_px"] = [pm6.width, pm6.height]
    return dict(pdf=pdf, png=png, receipt=rec)


# =======================================================================================
# Statistics: A24 terminal inference (moved from search.py)
# =======================================================================================
def _a24_seed(key):
    """Seed of an A24 random generator, derived from a key.

    The seed is the first 4 bytes, read little-endian, of
    ``sha256("|".join(map(str, key)) + "|" + A24_SEED_SALT)``. A24 uses two keys:
    ``["null", n]`` for the 200,000 random sign vectors of a test with n > ``EXACT_MAX_N``
    paired participants, and ``["bootstrap", n, draws]`` for the multinomial bootstrap
    weights. Both depend only on n (participants with a finite paired difference) and the
    number of draws, so every test with the same n reuses the same sign vectors and
    bootstrap weights. Statistical unit: participant. In Figure 6, n = 15, so the tests
    are exact and only the bootstrap seed is used.

    Parameters
    ----------
    key : list
        ``["null", n]`` or ``["bootstrap", n, draws]``.

    Returns
    -------
    int
        Seed in [0, 2**32) for ``numpy.random.default_rng``.
    """
    return int.from_bytes(hashlib.sha256(("|".join(map(str, key)) + "|" + A24_SEED_SALT).encode()).digest()[:4], "little")


def _a24_assess(delta, key, bootstrap, cache):
    """A24 ``assess``: paired sign-flip test and participant bootstrap interval of the mean."""
    delta = np.asarray(delta, float); delta = delta[np.isfinite(delta)]; n = len(delta)
    result = dict(n_participants=n, mean_difference=float(delta.mean()) if n else np.nan,
                  median_difference=float(np.median(delta)) if n else np.nan, n_positive=int((delta > 0).sum()),
                  n_negative=int((delta < 0).sum()), p=1., estimable=n >= MIN_PAIRED_N, CI95_low=np.nan, CI95_high=np.nan)
    if n < MIN_PAIRED_N:
        return result
    nc, bc = cache.setdefault("null", {}), cache.setdefault("boot", {})
    if n not in nc:
        if n <= EXACT_MAX_N:
            z = np.arange(2 ** n, dtype=np.uint32)[:, None]; nc[n] = ((z >> np.arange(n)) & 1).astype(np.int8) * 2 - 1
        else:
            nc[n] = np.random.default_rng(_a24_seed(["null", n])).choice(np.array([-1, 1], np.int8), size=(N_SIGN_DRAWS, n))
    means = nc[n] @ delta / n; hits = int((abs(means) >= abs(delta.mean()) - SIGN_TOL).sum())
    result["p"] = hits / len(means) if n <= EXACT_MAX_N else (hits + 1) / (len(means) + 1)
    result["test"] = "exact_sign_flip" if n <= EXACT_MAX_N else "200000_sign_flips_plus_one"
    result["null_size"] = len(means); result["hits"] = hits
    if (n, bootstrap) not in bc:
        bc[n, bootstrap] = np.random.default_rng(_a24_seed(["bootstrap", n, bootstrap])).multinomial(n, np.full(n, 1 / n), size=bootstrap)
    ci = np.quantile(bc[n, bootstrap] @ delta / n, CI_QUANTILES); result.update(CI95_low=float(ci[0]), CI95_high=float(ci[1]))
    result["draws"] = bootstrap
    return result


def bh(p):
    """Benjamini-Hochberg adjusted q values (the A24/A28 implementation).

    Parameters
    ----------
    p : array-like of float

    Returns
    -------
    numpy.ndarray
    """
    p = np.asarray(p, float); order = np.argsort(p); v = np.empty(len(p))
    v[order] = np.minimum(1, np.minimum.accumulate((p[order] * len(p) / np.arange(1, len(p) + 1))[::-1])[::-1])
    return v


def stat_terminal_families(pp):
    """Recompute the A24 CapScan families ``terminal_fraction__length_9`` and ``all18_nonKR_composition__length_9``.

    Statistical unit: participant (CapScan, n = 15 paired participants where both sides
    exist). For every route x unit x origin x terminus x metric cell and each of 14
    contrasts (10 region pairs, regions 1-2, regions 3-4 and capsule against stool, and
    regions 1-2 against regions 3-4), the paired participant differences are tested by a
    two-sided sign-flip test (exact enumeration of 2^n signs for n <= 16) and given a 95%
    percentile bootstrap interval of the mean (50,000 resamples for terminal fractions,
    2,000 for residue composition; weights drawn once per n from a seed derived from
    sha256). Cells with fewer than eight participants stay in the family with p = 1.
    q values are Benjamini-Hochberg within each family.

    Parameters
    ----------
    pp : pandas.DataFrame
        Both A24 participant-region tables. Participant ids are converted to text as in
        ``search.py``, so that bootstrap weights meet participants in the original order
        (lexicographic).

    Returns
    -------
    pandas.DataFrame
        One row per test, columns as the original family table plus ``null_size``, ``hits`` and ``draws``.
    """
    keys = ["cohort", "route", "min_length", "unit", "origin", "terminus", "metric"]
    pp = pp.assign(participant=pp.participant.astype(str))  # search.py: f.participant.astype(str)
    cap = pp[pp.origin.isin(PRIMARY_ORIGINS) & pp.cohort.eq("CAPSCAN") & pp.min_length.eq(MIN_LENGTH)]
    tests = []; cache = {}
    for labels, gg in cap.groupby(keys, sort=False):
        info = dict(zip(keys, labels)); is_residue = info["metric"].startswith("residue_")
        if info["metric"] == "other_residue_fraction_nonKR": continue
        if is_residue and (info["unit"] != "assigned_bare_sequence" or info["min_length"] != 9): continue
        if not is_residue and info["unit"] == "assigned_bare_sequence": continue
        pv = gg.pivot(index="participant", columns="region", values="value").reindex(columns=REGIONS)
        pv["proximal"] = pv[["duodenum", "jejunum"]].mean(axis=1, skipna=False); pv["distal"] = pv[["ileum", "colon"]].mean(axis=1, skipna=False)
        pv["capsule"] = pv[["proximal", "distal"]].mean(axis=1, skipna=False)
        contrasts = list(itertools.combinations(REGIONS, 2)) + [("proximal", "stool"), ("distal", "stool"), ("capsule", "stool"), ("proximal", "distal")]
        for a, b in contrasts:
            pair = pv[[a, b]].dropna(); delta = pair[a] - pair[b]
            family = ("all18_nonKR_composition" if is_residue else "terminal_fraction") + "__length_" + str(info["min_length"])
            tests.append(dict(info, contrast=a + "__minus__" + b, family=family,
                              **_a24_assess(delta, [*labels, a, b], A24_N_BOOT_COMPOSITION if is_residue else A24_N_BOOT_TERMINAL, cache)))
    tf = pd.DataFrame(tests); tf["q"] = tf.groupby("family").p.transform(bh)
    return tf


def compare_families(new, ref, keys, cols, atol=0.0):
    """Join recomputed and shipped family rows and report the largest absolute difference per column.

    Returns
    -------
    dict
        ``rows_new``, ``rows_ref``, ``rows_joined`` and ``max_abs_diff`` per column.
    """
    m = new.merge(ref, on=keys, suffixes=("_new", "_ref"), how="outer", indicator=True)
    out = dict(rows_new=len(new), rows_ref=len(ref), rows_joined=int((m["_merge"] == "both").sum()), max_abs_diff={})
    for c in cols:
        a = pd.to_numeric(m[c + "_new"], errors="coerce").to_numpy(float); b = pd.to_numeric(m[c + "_ref"], errors="coerce").to_numpy(float)
        both = np.isfinite(a) & np.isfinite(b)
        nan_mismatch = int((np.isfinite(a) != np.isfinite(b)).sum())
        out["max_abs_diff"][c] = (float(np.max(np.abs(a[both] - b[both]))) if both.any() else 0.0, nan_mismatch)
    return out


# =======================================================================================
# Statistics: A28 histidine inference (moved from histidine_helpers.py, collect_revision4.py)
# =======================================================================================
def signflip(delta, cache):
    """A28 ``signflip``: two-sided paired sign-flip test and 20,000-draw bootstrap interval of the mean.

    Exact enumeration for n <= 16; otherwise 200,000 random sign vectors from
    ``rng(A28_SEED_NULL + n)`` with plus-one correction. Bootstrap weights from
    ``rng(A28_SEED_BOOT + n)``, drawn once per n.

    Returns
    -------
    dict
    """
    d = np.asarray(delta, float); d = d[np.isfinite(d)]; n = len(d)
    if n < MIN_PAIRED_N: return dict(n=n, mean=float(d.mean()) if n else np.nan, p=1., estimable=False)
    rng = np.random.default_rng(A28_SEED_NULL + n)
    nc, bc = cache.setdefault("null", {}), cache.setdefault("boot", {})
    if n not in nc:
        if n <= EXACT_MAX_N:
            nc[n] = ((np.arange(2 ** n, dtype=np.uint32)[:, None] >> np.arange(n)) & 1).astype(np.int8) * 2 - 1
        else: nc[n] = rng.choice(np.array([-1, 1], np.int8), size=(N_SIGN_DRAWS, n))
    signs = nc[n]
    hits = int((np.abs(signs @ d / n) >= abs(d.mean()) - SIGN_TOL).sum())
    p = hits / len(signs) if n <= EXACT_MAX_N else (hits + 1) / (len(signs) + 1)
    if n not in bc: bc[n] = np.random.default_rng(A28_SEED_BOOT + n).multinomial(n, np.full(n, 1 / n), size=A28_N_BOOT)
    ci = np.quantile(bc[n] @ d / n, CI_QUANTILES)
    return dict(n=n, mean=float(d.mean()), p=p, estimable=True, CI95_low=ci[0], CI95_high=ci[1],
                n_positive=int((d > 0).sum()), n_negative=int((d < 0).sum()), null_size=len(signs), hits=hits,
                draws=int(bc[n].shape[0]))


def _ratio_interval(g):
    """Enrichment ratio (mean observed / mean expected fraction) with a bootstrap interval.

    Statistical unit: participant within one CapScan region (one row of ``g`` per
    participant, n up to 15). The point estimate is the ratio of the participant means of
    ``O`` (observed residue fraction at internal boundaries) and ``E`` (source-protein
    expectation), a ratio of means rather than a mean of ratios. Rows with a missing or
    infinite ``O`` or ``E`` are dropped first. If n >= ``MIN_PAIRED_N``, the interval is a
    percentile bootstrap: ``A28_N_BOOT`` (20,000) multinomial resamples of whole
    participants from ``rng(A28_SEED_RATIO + n)``; each resample gives a weighted
    observed and expected mean and their ratio; resamples with expected mean 0 are
    undefined and dropped; the interval is the ``CI_QUANTILES`` quantiles of the defined
    ratios. Used only in the A28 family ``source_P1_all_residues``; Figure 6 does not
    plot it.

    Parameters
    ----------
    g : pandas.DataFrame
        Participant rows of one route, side, origin, residue and region.

    Returns
    -------
    dict
        Empty if ``O``/``E`` are absent or no row remains; otherwise
        ``observed_fraction_mean``, ``expected_fraction_mean``, ``enrichment_ratio`` and,
        for n >= ``MIN_PAIRED_N``, ``ratio_bootstrap_defined`` and the interval bounds
        ``enrichment_CI95_low``/``enrichment_CI95_high``.
    """
    if not {"O", "E"} <= set(g.columns): return {}
    a = g[["O", "E"]].replace([np.inf, -np.inf], np.nan).dropna().values; n = len(a)
    if not n: return {}
    means = a.mean(axis=0); r = dict(observed_fraction_mean=means[0], expected_fraction_mean=means[1], enrichment_ratio=means[0] / means[1] if means[1] else np.nan)
    if n >= MIN_PAIRED_N:
        weights = np.random.default_rng(A28_SEED_RATIO + n).multinomial(n, np.full(n, 1 / n), size=A28_N_BOOT)
        v = weights @ a / n; ratios = np.divide(v[:, 0], v[:, 1], out=np.full(len(v), np.nan), where=v[:, 1] > 0)
        valid = ratios[np.isfinite(ratios)]; r["ratio_bootstrap_defined"] = len(valid)
        if len(valid): r.update(enrichment_CI95_low=np.quantile(valid, CI_QUANTILES[0]), enrichment_CI95_high=np.quantile(valid, CI_QUANTILES[1]))
    return r


def stat_histidine_source_family(pe):
    """Recompute the A28 CapScan family ``source_P1_all_residues`` (9,000 planned tests).

    Statistical unit: participant within region (CapScan, up to n = 15). Each cell
    (15 routes x 2 sides x 3 origins x 20 residues x 5 regions) tests whether the mean
    participant difference (observed residue fraction at internal boundaries minus the
    source-protein expectation) is zero by :func:`signflip`; unavailable cells are kept
    with p = 1. q values are Benjamini-Hochberg over the whole family.

    Parameters
    ----------
    pe : pandas.DataFrame
        ``capscan_histidine_participant_region_enrichment.tsv`` (file order preserved).

    Returns
    -------
    pandas.DataFrame
    """
    key = ["route", "side", "origin", "residue"]; cache = {}; tests = []
    lookup = {labels: g for labels, g in pe.groupby(key + ["region"], sort=False)}
    for route, side, origin, residue, region in itertools.product(A28_ROUTES_CAPSCAN, ["N", "C"], PRIMARY_ORIGINS, AA, REGIONS):
        g = lookup.get((route, side, origin, residue, region), pd.DataFrame({"difference": []}))
        tests.append(dict(endpoint="source_enrichment", family="source_P1_all_residues", route=route, side=side, origin=origin,
                          residue=residue, contrast=region, **signflip(g.difference, cache), **_ratio_interval(g)))
    tf = pd.DataFrame(tests); tf["q"] = bh(tf.p)
    return tf


def stat_histidine_reciprocal_family(pr):
    """Recompute the A28 CapScan family ``all_reciprocal_controls`` (90 tests).

    Statistical unit: participant within region (CapScan). For 6 routes x 3 origins x
    5 regions the participant difference (complementary-product rate at histidine
    boundaries minus matched non-histidine boundaries) is tested by :func:`signflip`;
    q values are Benjamini-Hochberg over the 90 tests.

    Parameters
    ----------
    pr : pandas.DataFrame
        Participant-level reciprocal-control rows, ``load(sd)["pr"]`` (``F_A28_RECIPROCAL``).

    Returns
    -------
    pandas.DataFrame
    """
    cells = {labels: g for labels, g in pr.groupby(["route", "origin", "region"], sort=False)}
    routes = [r for r in A28_PRIMARY if r in A28_ROUTES_CAPSCAN]; tests = []; cache = {}
    for labels in itertools.product(routes, PRIMARY_ORIGINS, REGIONS):
        g = cells.get(labels, pd.DataFrame({"difference": []}))
        tests.append(dict(endpoint="matched_reciprocal", family="all_reciprocal_controls", route=labels[0], origin=labels[1], contrast=labels[2], **signflip(g.difference, cache)))
    tf = pd.DataFrame(tests); tf["q"] = bh(tf.p)
    return tf


def stat_predict_patient(pt):
    """PREDICT patient-level histidine source enrichment (de novo and non-specific routes).

    Statistical unit: patient (n = 93; specimen means averaged within patient upstream).
    Mean, sign-flip P (200,000 random signs, plus-one) and 20,000-draw bootstrap interval
    by :func:`signflip`; ``n_positive`` counts patients with a positive difference. The
    q values need the full PREDICT family (1,560 tests over all residues); they are
    recomputed by :func:`bh` from the P values of ``predict_histidine_statistical_families.tsv``
    (the notebook checks them), not from patient values of other residues.

    Parameters
    ----------
    pt : pandas.DataFrame
        Patient-level PREDICT values, ``load(sd)["predict_patient"]`` (``F_PREDICT_PATIENT``).

    Returns
    -------
    pandas.DataFrame
        One row per route x side x origin, residue H.
    """
    rows = []; cache = {}
    for (route, side, origin, residue), g in pt.groupby(["route", "side", "origin", "residue"], sort=False):
        rows.append(dict(route=route, side=side, origin=origin, residue=residue, **signflip(g.difference, cache)))
    return pd.DataFrame(rows)


def predict_family_q(pfam):
    """Benjamini-Hochberg q values recomputed within each PREDICT A28 family from its P values.

    Parameters
    ----------
    pfam : pandas.DataFrame
        ``predict_histidine_statistical_families.tsv`` (1,560 source tests, 12 reciprocal tests).

    Returns
    -------
    pandas.DataFrame
        The table with a ``q_recomputed`` column.
    """
    return pfam.assign(q_recomputed=pfam.groupby("family").p.transform(bh))


def stat_fluid_observed_expected(a28):
    """Observed and source-expected histidine fractions at internal C boundaries in intestinal fluid.

    Statistical unit: participant (region means of participant values, from the A28
    family table columns ``observed_fraction_mean`` and ``expected_fraction_mean``).

    Parameters
    ----------
    a28 : pandas.DataFrame
        A28 family table, ``load(sd)["a28"]`` (``F_A28_FAMILIES``).

    Returns
    -------
    pandas.DataFrame
        Display routes x regions 1-4, origin All, residue H, side C, in percent.
    """
    x = a28[a28.endpoint.eq("source_enrichment") & a28.residue.eq("H") & a28.side.eq("C") & a28.origin.eq("All")
            & a28.route.isin(DISPLAY_ROUTES) & a28.contrast.isin(REGIONS[:4])]
    return x[["route", "contrast", "observed_fraction_mean", "expected_fraction_mean"]].assign(
        observed_pct=lambda d: 100 * d.observed_fraction_mean, expected_pct=lambda d: 100 * d.expected_fraction_mean)


# =======================================================================================
# Numbers and statistics tables (NUMBERS.tsv / STATS.tsv order)
# =======================================================================================
def collect_numbers(nums, d, tf, sf, rf, pt):
    """Every number of NUMBERS.tsv, in file order, as computed here.

    Parameters
    ----------
    nums : dict
        Plotted values collected by the panel functions.
    d : dict
        Output of :func:`load`.
    tf, sf, rf, pt : pandas.DataFrame
        Outputs of :func:`stat_terminal_families`, :func:`stat_histidine_source_family`,
        :func:`stat_histidine_reciprocal_family` and :func:`stat_predict_patient`.

    Returns
    -------
    list of tuple
        ``(row_id, value)``; values are full precision (participant counts in words where
        the v30 text prints words, digits where it prints digits: v30 paragraphs 85 and 89).
    """
    a, b, c = nums["a"], nums["b"], nums["c"]
    un = [a["unspecific"][r] for r in REGIONS[:4]]
    def row(fam_, **k):
        """The single row of ``fam_`` matching every column=value in ``k`` (asserts one)."""
        x = fam_
        for kk, v in k.items(): x = x[x[kk].eq(v)]
        assert len(x) == 1, k; return x.iloc[0]
    t_un = row(tf, route="unspecific", unit="canonical_peptide", origin="All", terminus=TERMINUS, metric="nonKR_fraction_all",
               contrast="capsule__minus__stool", family="terminal_fraction__length_9")
    t_h = [row(tf, route="unspecific", unit="canonical_peptide", origin=o, terminus=TERMINUS, metric="H_fraction_nonKR",
               contrast="capsule__minus__stool", family="terminal_fraction__length_9") for o in ("Human", "Microbial")]
    assert int(t_un.n_positive) == int(t_un.n_participants)
    oe = stat_fluid_observed_expected(d["a28"]); oe_u = oe[oe.route.eq("sage_unspecific_bothq")]
    pc = pt[pt.side.eq("C") & pt.origin.eq("All")]
    n_pairs = int(t_un.n_participants)
    words = lambda n: NUMBER_WORDS[int(n)]  # noqa: E731
    pred_n = int(d["S"][d["S"].cohort.eq("PREDICT")].n.unique()[0])
    stool_n = int(d["pp"].pipe(lambda x: x[x.route.eq("unspecific") & x.unit.eq("canonical_peptide") & x.origin.eq("All")
                                        & x.terminus.eq(TERMINUS) & x.metric.eq("nonKR_fraction_all") & x.region.eq("stool")]).value.notna().sum())
    fam_sizes = tf.groupby("family").size().to_dict(); a28_sizes = d["a28"].groupby("family").size().to_dict()
    pfam = d["predict_fam"].groupby("family").size().to_dict()
    n_pred_rec = int(pfam["all_reciprocal_controls"])
    lengths = sorted(d["pp"].min_length.unique())
    return [
        ("165_nonKR_min", float(min(un))), ("165_nonKR_max", float(max(un))), ("165_stool", float(a["unspecific"]["stool"])),
        ("165_n", int(t_un.n_positive)), ("165_tryptic", a["tryptic_mean_pct"]),
        ("169_min_length", int(lengths[0])), ("169_max_length", None), ("169_predict_n", pred_n),
        ("169_ci_level", 100 * (CI_QUANTILES[1] - CI_QUANTILES[0])), ("169_n", n_pairs),
        ("173_H_r12", float(c["proximal"]["H"])), ("173_H_r34", float(c["distal"]["H"])), ("173_H_stool", float(c["stool"]["H"])),
        ("173_n", min(int(x.n_positive) for x in t_h)),
        ("173_human_pp", float(100 * t_h[0].mean_difference)), ("173_microbial_pp", float(100 * t_h[1].mean_difference)),
        ("175_observed", float(oe_u.observed_pct.mean())), ("175_expected", float(oe_u.expected_pct.mean())),
        ("175_predict_positive", int(pc.n_positive.min())),
        ("285_capscan_acq", int(len(d["design"]))), ("285_predict_acq", None),
        ("287_n_words", words(n_pairs)),
        ("855_capscan_acq", int(len(d["design"]))), ("855_predict_acq", None),
        ("863_n_words", words(n_pairs)),
        ("869_boot_def", int(sf.draws.max())), ("869_exact_n", None),
        ("869_sign_draws", int(pt.null_size.max())), ("869_predict_n", pred_n),
        ("869_boot_c", int(c["draws"])), ("869_boot_terminal", int(tf[tf.family.eq("terminal_fraction__length_9")].draws.max())),
        ("869_n_words", words(n_pairs)), ("869_stool_words", words(stool_n)),
        ("893_terminal", int(fam_sizes["terminal_fraction__length_9"])), ("893_composition", int(fam_sizes["all18_nonKR_composition__length_9"])),
        ("893_source_capscan", len(sf)), ("893_source_predict", int(pfam["source_P1_all_residues"])), ("893_rec_capscan", len(rf)),
        ("893_rec_predict", n_pred_rec), ("893_spatial", int(a28_sizes["H_spatial_all_routes"])),
        ("893_motif", int(a28_sizes["all_motif_positions_residues"])),
    ]


def _r(v):
    """Full-precision text of a number (Python float repr)."""
    return repr(float(v))


STATS_COLUMNS = ["claim", "v28d_line", "test", "sided", "unit", "n", "distinct_or_repeated", "covariates", "correction",
                 "family", "family_size", "estimate", "uncertainty", "uncertainty_construction", "statistic", "df",
                 "exact_P", "effect_size", "effect_size_how", "draws", "seed", "software_version", "code"]


def stats_table(nums, d, tf, sf, rf, pt, software):
    """Rows of STATS.tsv (every test and interval Figure 6 shows or the text cites), from the code.

    Parameters
    ----------
    nums, d, tf, sf, rf, pt
        As in :func:`collect_numbers`.
    software : str
        Interpreter and package versions of this run.

    Returns
    -------
    pandas.DataFrame
        Columns ``STATS_COLUMNS``.
    """
    rows = []
    fam_size = tf.groupby("family").size().to_dict()
    def a24(route, origin, metric, unit, contrast, family, claim, line, code):
        """Append the STATS.tsv row of one A24 paired participant test (endpoint selected by the arguments)."""
        x = tf[tf.route.eq(route) & tf.origin.eq(origin) & tf.metric.eq(metric) & tf.unit.eq(unit) & tf.terminus.eq(TERMINUS)
               & tf.contrast.eq(contrast) & tf.family.eq(family)].iloc[0]
        n = int(x.n_participants); draws = int(x.draws)
        rows.append(dict(claim=claim, v28d_line=line, test="paired sign-flip test of the mean difference (exact enumeration of all 2^n sign vectors)",
                         sided="two-sided", unit="participant", n=str(n),
                         distinct_or_repeated="repeated: sampling regions and stool within participant; one paired difference per participant",
                         covariates="none", correction="Benjamini-Hochberg", family=family, family_size=str(fam_size[family]),
                         estimate=f"mean paired difference {_r(100 * x.mean_difference)} percentage points",
                         uncertainty=f"95% CI {_r(100 * x.CI95_low)} to {_r(100 * x.CI95_high)} percentage points",
                         uncertainty_construction=f"percentile bootstrap, {draws:,} draws, resampling participants (multinomial weights drawn once per n)",
                         statistic=f"mean paired difference; {int(x.hits)} of {int(x.null_size):,} sign vectors as extreme (|mean| >= observed - 1e-12)",
                         df="NA", exact_P=_r((x.p)), effect_size=f"{_r(100 * x.mean_difference)} percentage points; q = {_r(float(x.q))}; {int(x.n_positive)} of {n} positive",
                         effect_size_how="mean over participants of (capsule or region-group value minus stool value)",
                         draws=str(draws), seed=f"test: none (exact); bootstrap: {_a24_seed(['bootstrap', n, draws])} (sha256 of 'bootstrap|{n}|{draws}|{A24_SEED_SALT}', first 4 bytes little-endian)",
                         software_version=software, code=code))
    for route, lab in [("denovo", "de novo"), ("semi", "semi-tryptic"), ("unspecific", "non-specific")]:
        a24(route, "All", "nonKR_fraction_all", "canonical_peptide", "capsule__minus__stool", "terminal_fraction__length_9",
            f"non-K/R C-terminal fraction higher in intestinal fluid than stool, {lab} route (Fig. 6a)", "165;169",
            "figures/Fig6/panels.py:stat_terminal_families")
    for origin, lab in [("Human", "human"), ("Microbial", "microbial")]:
        a24("unspecific", origin, "H_fraction_nonKR", "canonical_peptide", "capsule__minus__stool", "terminal_fraction__length_9",
            f"histidine share among non-K/R C termini higher in intestinal fluid than stool, {lab}-assigned (Fig. 6b)", "169;173",
            "figures/Fig6/panels.py:stat_terminal_families")
    for k, lab in [("proximal", "regions 1-2"), ("distal", "regions 3-4")]:
        a24("unspecific", "All", "residue_H_fraction_nonKR", "assigned_bare_sequence", f"{k}__minus__stool", "all18_nonKR_composition__length_9",
            f"histidine share of non-K/R C termini, {lab} minus stool (Fig. 6c; not printed)", "173;893",
            "figures/Fig6/panels.py:stat_terminal_families")
    c = nums["c"]
    for k, lab in [("proximal", "regions 1-2"), ("distal", "regions 3-4"), ("stool", "stool")]:
        for aa in c["order"]:
            lo, hi = c["ci"][k][aa]
            rows.append(dict(claim=f"Fig. 6c residue {aa} share of non-K/R C termini, {lab}", v28d_line="169;173;869", test="none (interval only)",
                             sided="NA", unit="participant", n=str(c[k + "_n"]), distinct_or_repeated="distinct participants within group",
                             covariates="none", correction="none", family="NA", family_size="NA",
                             estimate=f"mean {_r(c[k][aa])}%", uncertainty=f"95% CI {_r(lo)} to {_r(hi)}%",
                             uncertainty_construction=f"percentile bootstrap (np.quantile at CI_QUANTILES 0.025/0.975), {c['draws']:,} resamples of participants with replacement",
                             statistic="NA", df="NA", exact_P="NA", effect_size="NA", effect_size_how="NA", draws=str(c["draws"]),
                             seed=f"{SEED_PLOT} (one numpy default_rng shared with the point jitter of panels a and b; draws taken in panel order)",
                             software_version=software, code="figures/Fig6/panels.py:panel_c"))
    def a28(row, fam, fam_n, claim, line, code, cohort, sided="two-sided"):
        """Append the STATS.tsv row of one A28 test.

        ``row`` is the test row (mean, CI, p, q, n); ``fam``/``fam_n`` its family and
        size; ``cohort`` CAPSCAN (unit participant) or PREDICT (unit patient).
        """
        n = int(row.n); exact = n <= EXACT_MAX_N
        rows.append(dict(claim=claim, v28d_line=line,
                         test="one-sample sign-flip test of participant values against zero (" + ("exact enumeration)" if exact else f"{N_SIGN_DRAWS:,} random sign vectors, plus-one correction)"),
                         sided=sided, unit="participant" if cohort == "CAPSCAN" else "patient", n=str(n),
                         distinct_or_repeated=("repeated: acquisitions within specimen and specimens within participant-region averaged first" if cohort == "CAPSCAN"
                                               else "repeated: specimens within patient averaged first"),
                         covariates="none", correction="Benjamini-Hochberg", family=fam, family_size=str(fam_n),
                         estimate=f"mean {_r(100 * row['mean'])} percentage points", uncertainty=f"95% CI {_r(100 * row.CI95_low)} to {_r(100 * row.CI95_high)} percentage points",
                         uncertainty_construction=f"percentile bootstrap, {A28_N_BOOT:,} draws, resampling whole {'participants' if cohort == 'CAPSCAN' else 'patients'} (multinomial weights drawn once per n)",
                         statistic=(f"mean; {int(row.hits)} of {int(row.null_size):,} sign vectors as extreme" if "hits" in row and pd.notna(row.get("hits")) else "mean"),
                         df="NA", exact_P=_r((row.p)), effect_size=f"{_r(100 * row['mean'])} percentage points; q = {_r(float(row.q))}" if "q" in row and pd.notna(row.get("q")) else f"{_r(100 * row['mean'])} percentage points",
                         effect_size_how=("observed histidine indicator minus source-protein expected fraction, averaged within specimen, then participant-region" if fam.startswith("source")
                                          else "complementary-product rate at histidine boundaries minus matched non-histidine boundaries"),
                         draws=str(A28_N_BOOT), seed=(f"test: none (exact)" if exact else f"test: {A28_SEED_NULL + n}") + f"; bootstrap: {A28_SEED_BOOT + n}",
                         software_version=software, code=code))
    regl = dict(zip(REGIONS, ["region 1", "region 2", "region 3", "region 4", "stool"]))
    for side, lab, panel in [("C", "C-terminal", "d"), ("N", "N-terminal", "e")]:
        for route in DISPLAY_ROUTES:
            for reg in REGIONS:
                x = sf[sf.route.eq(route) & sf.side.eq(side) & sf.origin.eq("All") & sf.residue.eq("H") & sf.contrast.eq(reg)].iloc[0]
                a28(x, "source_P1_all_residues (CapScan)", len(sf), f"Fig. 6{panel}: histidine enrichment at internal {lab} boundaries, CapScan {regl[reg]}, {route}",
                    "169;175;869;893", "figures/Fig6/panels.py:stat_histidine_source_family", "CAPSCAN")
        for route in ["denovo_top30", "sage_unspecific_bothq"]:
            x = pt[pt.route.eq(route) & pt.side.eq(side) & pt.origin.eq("All")].iloc[0].copy()
            pq = predict_family_q(d["predict_fam"])
            q = pq[pq.endpoint.eq("source_enrichment") & pq.route.eq(route) & pq.side.eq(side) & pq.origin.eq("All") & pq.residue.eq("H")].iloc[0].q_recomputed
            x["q"] = q
            a28(x, "source_P1_all_residues (PREDICT)", str(int((d["predict_fam"].family == "source_P1_all_residues").sum())), f"Fig. 6{panel}: histidine enrichment at internal {lab} boundaries, PREDICT stool, {route} (q: Benjamini-Hochberg over the shipped PREDICT family P values)",
                "169;175;869;893", "figures/Fig6/panels.py:stat_predict_patient", "PREDICT")
    for route in DISPLAY_ROUTES:
        for reg in REGIONS:
            x = rf[rf.route.eq(route) & rf.origin.eq("All") & rf.contrast.eq(reg)].iloc[0]
            a28(x, "all_reciprocal_controls (CapScan)", len(rf), f"Fig. 6f: complementary products at histidine minus matched boundaries, CapScan {regl[reg]}, {route}",
                "169;177;869;893", "figures/Fig6/panels.py:stat_histidine_reciprocal_family", "CAPSCAN")
    pr = d["R"][d["R"].cohort.eq("PREDICT")]
    for route in ["denovo_top30", "sage_unspecific_bothq"]:
        x = pr[pr.route.eq(route)].iloc[0]
        r = dict(claim=f"Fig. 6f: complementary products, PREDICT stool, {route} (values from the shipped PREDICT table; patient-level values not in this repository)",
                 v28d_line="169;177;869;893", test=f"one-sample sign-flip test ({N_SIGN_DRAWS:,} random sign vectors, plus-one correction)", sided="two-sided",
                 unit="patient", n=str(int(x.n)), distinct_or_repeated="repeated: specimens within patient averaged first", covariates="none",
                 correction="Benjamini-Hochberg", family="all_reciprocal_controls (PREDICT)", family_size=str(int((d["predict_fam"].family == "all_reciprocal_controls").sum())),
                 estimate=f"mean {_r(100 * x['mean'])} percentage points", uncertainty=f"95% CI {_r(100 * x.CI95_low)} to {_r(100 * x.CI95_high)} percentage points",
                 uncertainty_construction=f"percentile bootstrap, {A28_N_BOOT:,} draws, resampling whole patients", statistic="mean (sign-vector hit count not in the shipped table)",
                 df="NA", exact_P=_r((x.p)), effect_size=f"{_r(100 * x['mean'])} percentage points; q = {_r(float(x.q))}",
                 effect_size_how="complementary-product rate at histidine boundaries minus matched non-histidine boundaries",
                 draws=str(A28_N_BOOT), seed=f"test: {A28_SEED_NULL + int(x.n)}; bootstrap: {A28_SEED_BOOT + int(x.n)} (A28 code; not re-run here)",
                 software_version="not recomputed; value from the shipped PREDICT table (A28 PREDICT run, interpreter versions not recorded in this repository)", code="figures/Fig6/panels.py:load (table read, not recomputed)")
        rows.append(r)
    return pd.DataFrame(rows, columns=STATS_COLUMNS)
