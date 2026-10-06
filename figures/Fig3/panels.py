"""Figure 3: database transfer, pooling and quantitative agreement in healthy stool.

Panels
------
a      schematic (``artwork/panel_a_design.pdf``); numbers traced in ``ARTWORK_NUMBERS.md``.
b-e    ``panels_bcde``: transfer (b, c) and pooling (d, e), drawn on one canvas exactly as the
       original producer drew them.
f      ``panel_f``: pairwise agreement heatmaps for protein groups and KEGG orthologs.

Source data (``source_data/Fig3/``)
-----------------------------------
transfer_pool_comparisons.tsv   one row per target acquisition, experiment, arm and length
                                floor (panels b-e; Results lines 99, 103; Methods line 739).
                                v30 plots the min_length 9 rows (9-50 residues).
protein_group_pairwise.tsv      per acquisition pair Pearson r, CV and feature count for
                                protein groups (panel f, upper row; Results line 109)
ko_pairwise.tsv                 per acquisition pair r, CV and feature count for KEGG
                                orthologs (panel f, lower row; Results line 109; Methods 747)
design.tsv                      acquisition to sampling-design assignment and order (panel f)

The plotting code is moved from the original producers with the logic unchanged; only
input and output paths were made repository-relative and the producers' checks against
lab files were replaced by checks against the receipts in ``expected/``. The v30 producer
(``build_fig3_len950.py`` with ``panels_len950.py``) differs from the v28f producer only in
``MIN_LENGTH`` (9 instead of 5); its composite step is ported as ``assemble_v30``. See
``PROVENANCE.md`` for the producers.
"""

import hashlib
import os
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib as mpl  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib import font_manager  # noqa: E402

# ---------------------------------------------------------------------------------------
# Choices (README "Choices" table). Every value is the one the original producer used.
# ---------------------------------------------------------------------------------------
MIN_LENGTH = 9                     # v30 (len950): b-e use 9-50-residue peptides, the accepted 5-50 set
                                   # filtered to >= 9 residues, q not re-estimated (v30 legend; Supp.
                                   # Methods "reported at 9-50 residues"); v28d line 739 used 5
PEPTIDE_Q = 0.01                   # Methods, v28d line 739: peptide q <= 0.01 (applied upstream, A15)
TRANSFER_ARMS = ("T1_within_group", "T4_different_donor")   # Methods, v28d line 739
POOL_ARMS = ("own", "pool_all12", "pool_other11")            # Results, v28d line 103; legend 107
DESIGNS = ("TechnicalReplicates", "WorkflowReplicates", "OnePiece", "Intra", "Interindividual")
DESIGN_LABELS = ("Technical\nreplicates", "Workflow\nreplicates", "Spatial\nsampling",
                 "Intra-\nindividual", "Inter-\nindividual")   # Results, v28d line 99 (five series)
JITTER_SEED_BC = 3                 # producer build_figure3_cdef.py: rng(3 + design index), b and c
JITTER_SEED_DE = 9                 # producer build_figure3_cdef.py: rng(9 + arm index), d and e
JITTER_HALFWIDTH_BC = 0.05         # producer: uniform(-0.05, 0.05) horizontal jitter, b and c
JITTER_HALFWIDTH_DE = 0.08         # producer: uniform(-0.08, 0.08) horizontal jitter, d and e
PG_CLASSES = ("host", "microbial")                 # legend v28d line 107: host above, microbial below
KO_PROFILES = ("new_all", "new_other_unresolved")  # legend v28d line 107: all origins above,
                                                   # human-assigned excluded below (Methods 747)
HEATMAP_SCALE = (0.0, 1.0)         # producer build_fig3f_cr003.py: Normalize(0, 1)
HEATMAP_COLOURS = ("#FFFFFF", "#0072B2")   # producer: white to Okabe-Ito blue
SIZE_BCDE_MM = (180, 97)           # producer build_figure3_cdef.py canvas (mm)
SIZE_F_PT = (518.7401733398438, 211.8685302734375)  # producer build_fig3f_cr003.py canvas:
                                   # page width x (page height - 344 pt) of the FAR21 Figure 3
FONT_F_REGULAR = "soehne-buch.ttf"       # producer build_fig3f_cr003.py
FONT_F_BOLD = "soehne-halbfett.ttf"      # producer build_fig3f_cr003.py
FONT_BCDE = "soehne-buch.otf"            # producer style (fl_style 2026-08-02 via init_plotting)

_FONT_DIR = Path(os.environ.get("FASTALAKE_FONT_DIR") or os.environ.get("FL_STYLE_DIR") or "")


def fonts_present():
    """Return True if every font file the producers used is in ``FASTALAKE_FONT_DIR``."""
    return bool(str(_FONT_DIR)) and all(
        (_FONT_DIR / f).is_file() for f in (FONT_BCDE, FONT_F_REGULAR, FONT_F_BOLD))


def sha256_file(path):
    """Return the hex sha256 of a file."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def load(sd):
    """Read the four source tables.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/Fig3``.

    Returns
    -------
    dict of pandas.DataFrame
        Keys ``comparisons``, ``pg``, ``ko``, ``design``.
    """
    sd = Path(sd)
    return dict(
        comparisons=pd.read_csv(sd / "transfer_pool_comparisons.tsv", sep="\t"),
        pg=pd.read_csv(sd / "protein_group_pairwise.tsv", sep="\t"),
        ko=pd.read_csv(sd / "ko_pairwise.tsv", sep="\t"),
        design=pd.read_csv(sd / "design.tsv", sep="\t"),
    )


# ---------------------------------------------------------------------------------------
# Producer style for panels b-e. The b-e producer imported the FastaLake figure style as
# it stood on 2026-08-02 (sha256 58b64c34...; 0.6 pt lines, tight bounding box on save),
# which differs from the later style in style/fl_style.py (1 pt lines, no tight box).
# The functions it called (init, style, save and the lab helpers init_plotting and
# save_figure) are moved here unchanged so the panels render as they did.
# ---------------------------------------------------------------------------------------
class _ProducerStyle:
    """The 2026-08-02 FastaLake figure style used by the panel b-e producer (moved unchanged)."""

    INK = "#4e4f51"
    GREY = "#7f7f7f"
    GRID = "#E9E9E9"
    BLUE, ORANGE, GREEN, VERM = "#0072B2", "#E69F00", "#009E73", "#D55E00"
    FL = BLUE

    @staticmethod
    def init_plotting(font_path):
        """Lab ``shared_functions_plotting.init_plotting`` (moved unchanged)."""
        mpl.rcParams["figure.dpi"] = 300
        mpl.rcParams["pdf.fonttype"] = 42
        mpl.rcParams["ps.fonttype"] = 42
        mpl.rcParams["svg.fonttype"] = "none"
        if os.path.isfile(font_path):
            fe = font_manager.FontEntry(fname=font_path, name="soehne")
            font_manager.fontManager.ttflist.insert(0, fe)
            mpl.rcParams["font.family"] = fe.name
        else:
            warnings.warn("Söhne absent: panels b-e use the default typeface", stacklevel=2)
            mpl.rcParams["font.family"] = "sans-serif"
        for key in ("text.color", "axes.labelcolor", "axes.edgecolor", "xtick.color", "ytick.color"):
            mpl.rcParams[key] = "#4e4f51"
        for key in ("lines.linewidth", "axes.linewidth", "xtick.major.width",
                    "ytick.major.width", "xtick.minor.width", "ytick.minor.width"):
            mpl.rcParams[key] = 0.6

    @classmethod
    def init(cls):
        """``fl_style.init`` of 2026-08-02 (moved unchanged)."""
        cls.init_plotting(str(_FONT_DIR / FONT_BCDE))
        mpl.rcParams.update({"font.size": 7, "axes.titleweight": "normal",
                             "axes.titlesize": 7, "legend.fontsize": 6.5,
                             "xtick.major.size": 2.6, "ytick.major.size": 2.6,
                             "figure.dpi": 300, "savefig.dpi": 300,
                             "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
                             "svg.fonttype": "none"})
        mpl.rcParams.update({"mathtext.fontset": "custom",
                             "mathtext.rm": "soehne", "mathtext.it": "soehne",
                             "mathtext.bf": "soehne", "mathtext.sf": "soehne",
                             "mathtext.tt": "soehne", "mathtext.cal": "soehne",
                             "mathtext.default": "regular"})

    @staticmethod
    def style(ax):
        """``fl_style.style`` of 2026-08-02 (moved unchanged)."""
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        ax.tick_params(length=2.6, width=0.6, labelsize=7)
        ax.xaxis.label.set_size(7)
        ax.yaxis.label.set_size(7)
        return ax

    @staticmethod
    def _is_title_text(s):
        """True when a text string is longer than two characters after stripping, i.e. not a panel letter."""
        s = (s or "").strip().rstrip(").:").strip()
        return len(s) > 2

    @classmethod
    def enforce_nbt(cls, fig):
        """``fl_style.enforce_nbt`` of 2026-08-02 (moved unchanged): strip on-canvas titles."""
        stripped = []
        for ax in fig.axes:
            for t in (ax.title, getattr(ax, "_left_title", None), getattr(ax, "_right_title", None)):
                if t is not None and cls._is_title_text(t.get_text()):
                    stripped.append(t.get_text().strip())
                    t.set_text("")
        sup = getattr(fig, "_suptitle", None)
        if sup is not None and sup.get_text().strip():
            stripped.append(sup.get_text().strip())
            sup.set_text("")
        return stripped

    @classmethod
    def save(cls, fig, basepath, dpi=300, formats=("pdf", "png")):
        """``fl_style.save`` plus lab ``save_figure`` (moved unchanged): tight box, white face."""
        cls.enforce_nbt(fig)
        base = os.path.splitext(str(basepath))[0]
        for fmt in formats:
            fig.savefig(f"{base}.{fmt}", dpi=dpi, bbox_inches="tight",
                        facecolor="white", edgecolor="none")
        plt.close(fig)
        return base + ".pdf", base + ".png"


# ---------------------------------------------------------------------------------------
# Panels b-e
# ---------------------------------------------------------------------------------------
def panels_bcde(comparisons, out_stem):
    """Draw panels b, c, d and e on one canvas, as the original producer did.

    Moved from the producer ``build_figure3_cdef.py`` (code sha256 2ba6574e...). b and c
    plot count retention and identity recovery of the 58 target acquisitions with a T1
    (same-design donor) and a T4 (inter-individual-series donor) database; d and e plot
    the 12 technical replicates with their own database and the all-12 and other-11
    unions. Points are acquisitions, bars medians, grey lines join one acquisition.

    Parameters
    ----------
    comparisons : pandas.DataFrame
        ``transfer_pool_comparisons.tsv``.
    out_stem : pathlib.Path
        Output stem; writes ``<stem>.pdf`` and ``<stem>.png``.

    Returns
    -------
    dict
        Plotted medians: ``transfer[(field, arm, design)]`` and ``pool[(field, arm)]``.
    """
    fl = _ProducerStyle
    with mpl.rc_context():
        mpl.rcdefaults()
        fl.init()
        plt.rcParams.update({"axes.labelcolor": fl.INK, "text.color": fl.INK, "xtick.color": fl.INK,
                             "ytick.color": fl.INK, "font.size": 6.5, "axes.labelsize": 6.5,
                             "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6})
        d = comparisons[comparisons.min_length.eq(MIN_LENGTH)]
        assert ((d.overlap + d.gained).eq(d.peptides) & (d.overlap + d.lost).eq(d.own)).all()
        MMW = 1 / 25.4
        W_MM, H_MM = SIZE_BCDE_MM
        fig = plt.figure(figsize=(W_MM * MMW, H_MM * MMW))

        def ax_mm(x, y, w, h):
            """Add an axes placed by x, y (from top), width and height in millimetres."""
            return fig.add_axes([x / W_MM, 1 - (y + h) / H_MM, w / W_MM, h / H_MM])

        def panel(ax, letter):
            """Apply the producer axis style and write the panel letter at the top left."""
            fl.style(ax)
            ax.set_title(letter, loc="left", fontsize=8, color=fl.INK, pad=4)

        groups = list(DESIGNS)
        labels = list(DESIGN_LABELS)
        c = d[d.kind.eq("transfer")]
        axC = ax_mm(12, 5, 76, 36)
        axD = ax_mm(104, 5, 76, 36)
        panel(axC, "b")
        panel(axD, "c")
        med = {}
        for a, field, label in [
            (axC, "count_retention_percent", "Peptides found with the transferred\ndatabase (% of own count)"),
            (axD, "own_set_recovery_percent", "Own peptide identities found again\nwith the transferred database (%)"),
        ]:
            for arm, offset, col, lab in [("T1_within_group", -.16, fl.FL, "Same-design donor (T1)"),
                                          ("T4_different_donor", .16, fl.ORANGE, "Designated donor control (T4)")]:
                for i, g in enumerate(groups):
                    v = c[c.arm.eq(arm) & c.group.eq(g)][field].to_numpy()
                    med[(field, arm, g)] = float(np.median(v))
                    a.scatter(i + offset + np.random.default_rng(JITTER_SEED_BC + i).uniform(
                        -JITTER_HALFWIDTH_BC, JITTER_HALFWIDTH_BC, len(v)), v, s=7, color=col, alpha=.65, lw=0)
                    a.plot([i + offset - .11, i + offset + .11], [np.median(v)] * 2, color=col, lw=1.5,
                           label=lab if i == 0 else None)
            a.set_xticks(range(5), labels)
            a.tick_params(axis="x", length=0)
            a.set_ylabel(label)
            a.set_ylim(0, 112)
            a.axhline(100, color=fl.GREY, lw=.6, ls="--")
        axC.legend(frameon=False, fontsize=5.8, loc="lower left", handlelength=1.2)
        p = d[d.kind.eq("pool")]
        arms = list(POOL_ARMS)
        pl = ["Own\ndatabase", "All 12\npooled", "Other 11\npooled"]
        cols = [fl.GREY, fl.GREEN, fl.FL]
        axE = ax_mm(12, 54, 60, 34)
        axF = ax_mm(104, 54, 60, 34)
        panel(axE, "d")
        panel(axF, "e")
        pm = {}
        for a, field, label in [(axE, "count_retention_percent", "Peptides found\n(% of own-database count)"),
                                (axF, "own_set_recovery_percent", "Own peptide identities\nfound again (%)")]:
            w = p.pivot(index="sample", columns="arm", values=field).reindex(columns=arms)
            assert w.shape == (12, 3) and w.notna().all().all()
            for _, row in w.iterrows():
                a.plot(range(3), row.to_numpy(), color=fl.GRID, lw=.5, zorder=0)
            for i, (arm, col) in enumerate(zip(arms, cols)):
                v = w[arm].to_numpy()
                a.scatter(i + np.random.default_rng(JITTER_SEED_DE + i).uniform(
                    -JITTER_HALFWIDTH_DE, JITTER_HALFWIDTH_DE, 12), v, s=9, color=col, alpha=.8, lw=0, zorder=2)
                a.plot([i - .2, i + .2], [np.median(v)] * 2, color=fl.INK, lw=1.1, zorder=3)
                pm[(field, arm)] = float(np.median(v))
            a.set_xticks(range(3), pl)
            a.tick_params(axis="x", length=0)
            a.set_ylabel(label)
            a.axhline(100, color=fl.GREY, lw=.6, ls="--")
            lo, hi = float(w.min().min()), float(w.max().max())
            a.set_ylim(lo - 2, hi + 2)
            a.yaxis.set_major_locator(plt.MaxNLocator(integer=True))
        plt.rcParams["savefig.bbox"] = None
        Path(out_stem).parent.mkdir(parents=True, exist_ok=True)
        fl.save(fig, out_stem, formats=("pdf", "png"))
    return dict(transfer=med, pool=pm)


# ---------------------------------------------------------------------------------------
# Panel f
# ---------------------------------------------------------------------------------------
def panel_f(pg, ko, design, out_stem):
    """Draw panel f, the pairwise-agreement heatmaps, as the original producer did.

    Moved from the producer ``build_fig3f_cr003.py`` (run as ``verified50``). Upper row:
    protein groups, host pairs above the diagonal and microbial pairs below. Lower row:
    KEGG orthologs, all origins above and human-assigned signal excluded below. Each tile
    is one acquisition pair; labels give the median r, median CV and median number of
    jointly quantified features across the pairs of one design and class.

    Parameters
    ----------
    pg, ko, design : pandas.DataFrame
        ``protein_group_pairwise.tsv``, ``ko_pairwise.tsv`` and ``design.tsv``.
    out_stem : pathlib.Path
        Output stem; writes ``<stem>.pdf`` and ``<stem>.png`` (300 dpi).

    Returns
    -------
    dict
        ``summary_values`` (20 records as the producer's receipt lists them) and
        ``heatmap_matrices`` (10 records with the matrix sha256).
    """
    old, new = pg, ko
    regular, bold = _FONT_DIR / FONT_F_REGULAR, _FONT_DIR / FONT_F_BOLD
    with mpl.rc_context():
        mpl.rcdefaults()
        have_font = regular.is_file() and bold.is_file()
        if have_font:
            font_manager.fontManager.addfont(str(regular))
            font_manager.fontManager.addfont(str(bold))
            family = font_manager.FontProperties(fname=str(regular)).get_name()
            boldprop = font_manager.FontProperties(fname=str(bold))
        else:
            warnings.warn("Söhne absent: panel f uses the default typeface", stacklevel=2)
            family = mpl.rcParams["font.family"]
            boldprop = font_manager.FontProperties(weight="bold")
        plt.rcParams.update({"font.family": family,
                             "font.size": 6.5, "text.color": "#4e4f51", "axes.labelcolor": "#4e4f51",
                             "pdf.fonttype": 42,
                             "axes.labelsize": 7, "xtick.labelsize": 6.5, "ytick.labelsize": 6.5,
                             "savefig.bbox": None})
        tiers = list(DESIGNS)
        labels = list(DESIGN_LABELS)
        cmap = mpl.colors.LinearSegmentedColormap.from_list("OkabeItoBlueSequential", list(HEATMAP_COLOURS))
        norm = mpl.colors.Normalize(*HEATMAP_SCALE)
        fig = plt.figure(figsize=(SIZE_F_PT[0] / 72, SIZE_F_PT[1] / 72))
        checked, matrices, texts = [], [], []
        for row, space in enumerate(["groups", "KO"]):
            for col, tier in enumerate(tiers):
                if space == "groups":
                    b = old[old.space.eq(space) & old.design.eq(tier) & old.cls.isin(list(PG_CLASSES))].copy()
                    classes = list(PG_CLASSES)
                else:
                    b = new[new.design.eq(tier) & new.profile.isin(list(KO_PROFILES))].copy()
                    b["cls"] = b.profile.map({KO_PROFILES[0]: "all", KO_PROFILES[1]: "other"})
                    classes = ["all", "other"]
                samples = design.loc[design.group.eq(tier), "sample_id"].tolist()
                assert set(samples) == set(b.i) | set(b.j)
                idx = {s: i for i, s in enumerate(samples)}
                m = np.full((len(samples), len(samples)), np.nan)
                assert len(b) == len(samples) * (len(samples) - 1)
                for record in b.itertuples():
                    i, j = sorted([idx[record.i], idx[record.j]])
                    assert 0 <= record.r_per_pair <= 1
                    if record.cls == classes[0]:
                        m[i, j] = record.r_per_pair
                    else:
                        m[j, i] = record.r_per_pair
                assert np.isfinite(m).sum() == len(b)
                ax = fig.add_axes([.07 + col * .165, .55 if row == 0 else .12, .145, .30])
                im = ax.imshow(m, cmap=cmap, norm=norm, interpolation="nearest")
                ax.set_axis_off()
                assert np.array_equal(np.ma.filled(im.get_array(), np.nan), m, equal_nan=True)
                mh = hashlib.sha256(m.tobytes()).hexdigest()
                matrices.append(dict(space=space, tier=tier, samples=samples, shape=list(m.shape), matrix_sha256=mh))
                for cls, xy, align in [(classes[0], (.98, .82), "right"), (classes[1], (.02, .18), "left")]:
                    s = b[b.cls.eq(cls)]
                    rmed = s.r_per_pair.median()
                    cv = s.cv_per_pair_pct.median()
                    n = s.n_features_per_pair.median()
                    record = dict(space=space, tier=tier, cls=cls, median_r=rmed, median_cv=cv, median_features=n)
                    texts.append(ax.text(*xy, f"r {rmed:.2f}\nCV {cv:.0f}%", transform=ax.transAxes, ha=align,
                                         va="center", fontsize=6.5,
                                         bbox=dict(facecolor="white", alpha=.8, edgecolor="none", pad=.5)))
                    upper = cls == classes[0]
                    texts.append(ax.text(.5, 1.04 if upper else -.10, f"{n:,.0f} " + cls, transform=ax.transAxes,
                                         ha="center", va="bottom" if upper else "top", fontsize=6.5))
                    checked.append(record)
                if row == 0:
                    texts.append(ax.text(.5, 1.20, labels[col], transform=ax.transAxes, ha="center",
                                         va="bottom", fontsize=7))
        fig.text(.018, .71, "Protein groups", rotation=90, ha="center", va="center", fontsize=7)
        fig.text(.018, .28, "Functions (KO)", rotation=90, ha="center", va="center", fontsize=7)
        fig.text(.01, .965, "f", fontsize=8, fontproperties=boldprop)
        fig.text(.925, .88, "Class labels:\nabove / below", ha="center", fontsize=6.5)
        cax = fig.add_axes([.91, .16, .022, .64])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), cax=cax, ticks=[0, .25, .5, .75, 1],
                          label="Pairwise Pearson r")
        cb.ax.tick_params(labelsize=6.5)
        cb.set_label("Pairwise Pearson r", size=7)
        fig.text(.5, .025, "KO all: all origins. KO other: human-assigned signal excluded; unresolved origin retained.",
                 ha="center", fontsize=6.5)
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()
        for text in texts + fig.texts:
            box = text.get_window_extent(renderer)
            assert box.x0 >= 0 and box.x1 <= fig.bbox.x1 and box.y0 >= 0 and box.y1 <= fig.bbox.y1, (text.get_text(), box)
        Path(out_stem).parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(str(out_stem) + ".pdf")
        fig.savefig(str(out_stem) + ".png", dpi=300)
        plt.close(fig)
    return dict(summary_values=checked, heatmap_matrices=matrices)


def check_panel_f_receipt(result, receipt):
    """Compare ``panel_f`` output with the producer's receipt in ``expected/``.

    Parameters
    ----------
    result : dict
        Return value of ``panel_f``.
    receipt : dict
        ``expected/panel_f_producer_receipt.json`` (summary values and matrix sha256 the
        producer recorded).

    Returns
    -------
    dict
        Counts of equal summary records and equal matrices.
    """
    sv = [r == e for r, e in zip(result["summary_values"], receipt["summary_values"])]
    mx = [r["matrix_sha256"] == e["matrix_sha256"] and r["samples"] == e["samples"]
          for r, e in zip(result["heatmap_matrices"], receipt["heatmap_matrices"])]
    return dict(summary_equal=sum(sv), summary_total=len(receipt["summary_values"]),
                matrices_equal=sum(mx), matrices_total=len(receipt["heatmap_matrices"]))


# ---------------------------------------------------------------------------------------
# Statistics. All are descriptive (medians and counts); Figure 3 has no hypothesis test.
# ---------------------------------------------------------------------------------------
def _q(v):
    """Median (numpy.median), quartiles (numpy.percentile, linear interpolation), min, max and n of v."""
    v = np.asarray(v, dtype=float)
    return dict(median=float(np.median(v)), q1=float(np.percentile(v, 25)), q3=float(np.percentile(v, 75)),
                minimum=float(v.min()), maximum=float(v.max()), n=int(len(v)))


def stat_a_design_counts(comparisons):
    """Acquisitions per sampling series (Fig. 3a; Results line 99; Methods line 739).

    Statistical unit: acquisition. n = 58 target acquisitions in five series.

    Parameters
    ----------
    comparisons : pandas.DataFrame

    Returns
    -------
    dict
        ``per_design`` (design -> acquisitions), ``total``, ``n_designs``.
    """
    t = comparisons[comparisons.kind.eq("transfer") & comparisons.min_length.eq(MIN_LENGTH)
                    & comparisons.arm.eq("own")]
    per = t.groupby("group")["sample"].nunique().reindex(list(DESIGNS))
    return dict(per_design={k: int(v) for k, v in per.items()}, total=int(t["sample"].nunique()),
                n_designs=int(per.notna().sum()))


def stat_bc_transfer(comparisons):
    """Medians and quartiles of transfer count retention and identity recovery (Fig. 3b,c).

    Statistical unit: target acquisition; n per design 12, 15, 10, 11, 10 (58 in all),
    each with one T1 and one T4 donor database. Medians are ``numpy.median``; quartiles
    are ``numpy.percentile`` with linear interpolation.

    Parameters
    ----------
    comparisons : pandas.DataFrame

    Returns
    -------
    pandas.DataFrame
        One row per arm, design and field.
    """
    c = comparisons[comparisons.kind.eq("transfer") & comparisons.min_length.eq(MIN_LENGTH)]
    rows = []
    for arm in TRANSFER_ARMS:
        for g in DESIGNS:
            for field in ("count_retention_percent", "own_set_recovery_percent"):
                rows.append(dict(arm=arm, design=g, field=field, **_q(c[c.arm.eq(arm) & c.group.eq(g)][field])))
    return pd.DataFrame(rows)


def stat_de_pool(comparisons):
    """Pooling medians, quartiles, gains and losses (Fig. 3d,e; Results line 103).

    Statistical unit: acquisition; n = 12 technical replicates of one digest (repeated
    injections, not independent participants). Each acquisition is searched against its
    own database, the union of all 12 and the union of the other 11.

    Parameters
    ----------
    comparisons : pandas.DataFrame

    Returns
    -------
    dict
        Per arm: count retention and identity recovery summaries, median gained and lost
        peptides, and the number of acquisitions with more peptides than with their own
        database.
    """
    p = comparisons[comparisons.kind.eq("pool") & comparisons.min_length.eq(MIN_LENGTH)]
    out = {}
    for arm in POOL_ARMS:
        a = p[p.arm.eq(arm)]
        out[arm] = dict(count=_q(a.count_retention_percent), identity=_q(a.own_set_recovery_percent),
                        median_gained=float(a.gained.median()), median_lost=float(a.lost.median()),
                        acquisitions_with_more_peptides=int((a.peptides > a.own).sum()),
                        n=int(a["sample"].nunique()))
    return out


def stat_methods_donor_assignments(comparisons):
    """T1 and T4 donor agreement for the inter-individual targets (Methods line 739).

    Statistical unit: target acquisition; n = 10 inter-individual targets.

    Parameters
    ----------
    comparisons : pandas.DataFrame

    Returns
    -------
    dict
        ``targets``, ``differ``, ``coincide`` and the coinciding (target, donor) pairs.
    """
    c = comparisons[comparisons.kind.eq("transfer") & comparisons.min_length.eq(MIN_LENGTH)
                    & comparisons.group.eq("Interindividual") & comparisons.arm.isin(list(TRANSFER_ARMS))]
    w = c.pivot(index="sample", columns="arm", values="donor")
    same = w[w[TRANSFER_ARMS[0]] == w[TRANSFER_ARMS[1]]]
    return dict(targets=int(len(w)), differ=int((w[TRANSFER_ARMS[0]] != w[TRANSFER_ARMS[1]]).sum()),
                coincide=int(len(same)),
                coinciding=[(t, d) for t, d in same[TRANSFER_ARMS[0]].items()])


def _pair_summary(s):
    """Medians (pandas) and quartiles (linear interpolation) of r and CV over pairs, median features."""
    return dict(median_r=s.r_per_pair.median(), median_cv=s.cv_per_pair_pct.median(),
                median_features=s.n_features_per_pair.median(),
                r_q1=s.r_per_pair.quantile(.25), r_q3=s.r_per_pair.quantile(.75),
                cv_q1=s.cv_per_pair_pct.quantile(.25), cv_q3=s.cv_per_pair_pct.quantile(.75))


def stat_f_summaries(pg, ko):
    """Median r, CV and joint feature count per design and class (Fig. 3f; Results line 109).

    Statistical unit: acquisition pair (pairs share acquisitions, so they are dependent);
    n pairs per design 66, 105, 45, 55, 45 from 12, 15, 10, 11, 10 acquisitions. Per pair,
    r is ``numpy.corrcoef`` of log2 intensities and CV the median over features of the
    population SD (ddof 0) divided by the mean x 100, both over features positive in both
    acquisitions (computed upstream; see UPSTREAM.json). Medians here are pandas medians.

    Parameters
    ----------
    pg, ko : pandas.DataFrame

    Returns
    -------
    pandas.DataFrame
        One row per space, design and class: n_pairs, median r, CV and features, and
        the quartiles of r and CV across pairs.
    """
    rows = []
    for g in DESIGNS:
        for cls in PG_CLASSES:
            s = pg[pg.space.eq("groups") & pg.design.eq(g) & pg.cls.eq(cls)]
            rows.append(dict(space="protein groups", design=g, cls=cls, n_pairs=len(s), **_pair_summary(s)))
    for g in DESIGNS:
        for prof in KO_PROFILES:
            s = ko[ko.design.eq(g) & ko.profile.eq(prof)]
            rows.append(dict(space="KO", design=g, cls=prof, n_pairs=len(s), **_pair_summary(s)))
    return pd.DataFrame(rows)


def stat_f_ko_vs_microbial_pg(summaries):
    """KO median r against microbial protein-group median r in every design (Results line 109).

    Statistical unit: design (five designs), comparing two medians over the same
    acquisition pairs; no test is made.

    Parameters
    ----------
    summaries : pandas.DataFrame
        Output of ``stat_f_summaries``.

    Returns
    -------
    pandas.DataFrame
        One row per design and KO profile with both medians and whether KO r is larger.
    """
    pgm = summaries[summaries.space.eq("protein groups") & summaries.cls.eq("microbial")].set_index("design").median_r
    k = summaries[summaries.space.eq("KO")].copy()
    k["microbial_pg_median_r"] = k.design.map(pgm)
    k["ko_greater"] = k.median_r > k.microbial_pg_median_r
    return k[["design", "cls", "median_r", "microbial_pg_median_r", "ko_greater"]].reset_index(drop=True)


_WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven", 8: "eight",
          9: "nine", 10: "ten", 11: "eleven", 12: "twelve"}


def as_word(n):
    """Return the English word the manuscript uses for a count from 1 to 12."""
    return _WORDS[int(n)]


def numbers(d):
    """Compute every value in ``NUMBERS.tsv``, in file order.

    Parameters
    ----------
    d : dict
        Output of ``load``.

    Returns
    -------
    list of tuple
        ``(key, full-precision value as str)``; keys match the ``code`` column suffix.
    """
    cmp_ = d["comparisons"]
    a = stat_a_design_counts(cmp_)
    bc = stat_bc_transfer(cmp_).set_index(["arm", "design", "field"])
    de = stat_de_pool(cmp_)
    don = stat_methods_donor_assignments(cmp_)
    f = stat_f_summaries(d["pg"], d["ko"]).set_index(["space", "design", "cls"])
    kv = stat_f_ko_vs_microbial_pg(stat_f_summaries(d["pg"], d["ko"]))
    t1 = "T1_within_group"
    out = [
        ("L99 acquisitions", a["total"]),
        ("L99 series", as_word(a["n_designs"])),
        ("L99 technical count", bc.loc[(t1, "TechnicalReplicates", "count_retention_percent"), "median"]),
        ("L99 technical identity", bc.loc[(t1, "TechnicalReplicates", "own_set_recovery_percent"), "median"]),
        ("L99 interindividual count", bc.loc[(t1, "Interindividual", "count_retention_percent"), "median"]),
        ("L99 interindividual identity", bc.loc[(t1, "Interindividual", "own_set_recovery_percent"), "median"]),
        ("L103 technical replicates", de["pool_other11"]["n"]),
        ("L103 other injections", de["pool_other11"]["n"] - 1),
        ("L103 median gain %", de["pool_other11"]["count"]["median"] - 100),
        ("L103 identity retained", de["pool_other11"]["identity"]["median"]),
        ("L103 median gained", de["pool_other11"]["median_gained"]),
        ("L103 median lost", de["pool_other11"]["median_lost"]),
        ("L107 five series (a)", as_word(a["n_designs"])),
        ("L107 ten-donor series", as_word(a["per_design"]["Interindividual"])),
        ("L107 all 12", de["pool_all12"]["n"]),
        ("L107 other 11", de["pool_other11"]["n"] - 1),
        ("L107 five series (f)", as_word(a["n_designs"])),
        ("L109 five series", as_word(a["n_designs"])),
        ("L109 microbial PG CV technical", f.loc[("protein groups", "TechnicalReplicates", "microbial"), "median_cv"]),
        ("L109 microbial PG CV interindividual", f.loc[("protein groups", "Interindividual", "microbial"), "median_cv"]),
        ("L109 microbial PG r technical", f.loc[("protein groups", "TechnicalReplicates", "microbial"), "median_r"]),
        ("L109 microbial PG r interindividual", f.loc[("protein groups", "Interindividual", "microbial"), "median_r"]),
        ("L109 KO r technical", f.loc[("KO", "TechnicalReplicates", "new_all"), "median_r"]),
        ("L109 KO r interindividual", f.loc[("KO", "Interindividual", "new_all"), "median_r"]),
        ("L739 acquisitions", a["total"]),
        ("L739 inter-individual acquisitions", as_word(a["per_design"]["Interindividual"])),
        ("L739 inter-individual targets", as_word(don["targets"])),
        ("L739 differ", as_word(don["differ"])),
        ("L739 coincide", as_word(don["coincide"])),
        # v30-only values (NUMBERS.tsv rows with NA v28d columns): v30 paragraph 50 "(+10.7%, 97.5% retained)"
        ("v30 P50 all-12 gain %", de["pool_all12"]["count"]["median"] - 100),
        ("v30 P50 all-12 identity retained", de["pool_all12"]["identity"]["median"]),
    ]
    assert bool(kv.ko_greater.all()), "KO r not above microbial protein-group r in every design"
    return [(k, v if isinstance(v, str) else repr(float(v)) if isinstance(v, float) else str(v)) for k, v in out]


# ---------------------------------------------------------------------------------------
# Composite placement checks (visual). The v28f composite was assembled by PDF patching
# (ASSEMBLY.md). These checks place a regenerated panel where the assembly chain put it
# and compare that band of the page with the v28f composite at 300 dpi. Differences are
# expected exactly where the chain edited text; the receipts and diff images record them.
# ---------------------------------------------------------------------------------------
COMPOSITE_PAGE_PT = (518.7401733398438, 570.8800048828125)   # expected/Fig3_v28f.pdf
SCALE_183_180 = 183 / 180          # fix_current_figures.py: 180 mm design exported at 183 mm
A_BAND_PT = 91.0                   # fig3_restore/compose_fig3.py: panel a band kept from JAG19
BCDE_TOP_MM = 33                   # compose_figure3.py: b-e canvas placed 33 mm from the top
BCDE_NATIVE_WIDTH_PT = 514.9598388671875   # width of the producer's tight-cropped b-e PDF
F_TOP_PT = 344.0 + 24.0            # build_fig3f_cr003.py fstart 344 pt, then +24 pt in compose_fig3.py
A_NATIVE_WIDTH_PT = 510.2362060546875      # artwork/panel_a_design.pdf width (180 mm)


def composite_band_check(panel_pdf, composite_pdf, rect, out_stem, work_dir):
    """Place ``panel_pdf`` into ``rect`` on a page the size of the composite and compare that band.

    Parameters
    ----------
    panel_pdf, composite_pdf : pathlib.Path
    rect : tuple of float
        (x0, y0, x1, y1) in PDF points where the assembly chain placed the panel.
    out_stem : pathlib.Path
        Receipt stem passed to ``tools.check_pixels.compare``.
    work_dir : pathlib.Path
        Directory for the two band renders (``build/``, gitignored).

    Returns
    -------
    dict
        The pixel-check receipt.
    """
    import fitz
    from tools.check_pixels import compare

    work_dir = Path(work_dir)
    work_dir.mkdir(parents=True, exist_ok=True)
    comp = fitz.open(composite_pdf)
    w, h = comp[0].rect.width, comp[0].rect.height
    doc = fitz.open()
    page = doc.new_page(width=w, height=h)
    page.show_pdf_page(fitz.Rect(*rect), fitz.open(panel_pdf), 0)
    clip = fitz.Rect(0, rect[1], w, min(rect[3], h))
    name = Path(out_stem).name
    new_png, exp_png = work_dir / f"{name}_regenerated.png", work_dir / f"{name}_v28f.png"
    page.get_pixmap(dpi=300, alpha=False, clip=clip).save(new_png)
    comp[0].get_pixmap(dpi=300, alpha=False, clip=clip).save(exp_png)
    return compare(new_png, exp_png, out_stem, dpi=300)


def composite_rects():
    """Return the placement rectangle of panel a, b-e and f in the v28f composite (points)."""
    w = COMPOSITE_PAGE_PT[0]
    a_h = 85.03936767578125 * w / A_NATIVE_WIDTH_PT
    y0 = BCDE_TOP_MM * 72 / 25.4 * SCALE_183_180
    bcde_h = 266.02435302734375 * w / BCDE_NATIVE_WIDTH_PT
    return dict(a=(0.0, 0.0, w, a_h),
                bcde=(0.0, y0, w, y0 + bcde_h),
                f=(0.0, F_TOP_PT, w, F_TOP_PT + SIZE_F_PT[1]))


def render_artwork(pdf, out_png, dpi=300):
    """Rasterise a schematic PDF from ``artwork/`` for display (no change to the artwork).

    Parameters
    ----------
    pdf, out_png : pathlib.Path
    dpi : int

    Returns
    -------
    pathlib.Path
        ``out_png``.
    """
    import fitz

    Path(out_png).parent.mkdir(parents=True, exist_ok=True)
    fitz.open(pdf)[0].get_pixmap(dpi=dpi, alpha=False).save(out_png)
    return Path(out_png)


#: Width in pixels of the preview images stored in the notebook (the full-resolution
#: renders stay in build/; a preview keeps the notebook small).
PREVIEW_WIDTH_PX = 900


def preview(png, width_px=PREVIEW_WIDTH_PX):
    """Write a downscaled copy ``<png stem>_preview.png`` beside ``png`` and return its path.

    Display only; pixel checks always use the full-resolution renders.
    """
    from PIL import Image

    png = Path(png)
    im = Image.open(png)
    h = round(im.height * width_px / im.width)
    out = png.with_name(png.stem + "_preview.png")
    im.convert("RGB").resize((width_px, h), Image.LANCZOS).save(out)
    return out


# ---------------------------------------------------------------------------------------
# v30 composite. Moved from the v30 producer ``build_fig3_len950.py`` (logic unchanged):
# the b-e band of the relabelled v28f page (expected/Fig3_v28f_relabel_base.pdf) is covered
# with white, the 9-50 b-e canvas is placed at composite_rects()['bcde'], its text is
# re-typeset as the v28f chain did, and the b legend is relabelled "Same series" /
# "Different donor".
# ---------------------------------------------------------------------------------------
V30_BAND_ROWS_300DPI = (333, 1520)                 # build_fig3_len950.py: replaced rows at 300 dpi
V30_LEGEND_OLD = ("Same-design donor (T1)", "Designated donor control (T4)")  # producer legend, not retyped
V30_RELABEL_BOX = (54.0, 186.8, 156.0, 204.7)      # fig3_labels/relabel.py: covered legend box (pt)
V30_RELABEL = ((193.862, "Same series"), (202.560, "Different donor"))  # relabel.py: baseline y (pt), text
V30_RELABEL_X = 54.517                             # relabel.py: legend text x (pt)
V30_RELABEL_SIZE = 6.5                             # relabel.py: Soehne Buch 6.5 pt
V30_LETTER_SIZE, V30_NUMERIC_SIZE, V30_TEXT_SIZE = 8, 6.5, 7   # build_FAR21.py re-typesetting sizes (pt)
V30_INK = (0x4e / 255, 0x4f / 255, 0x51 / 255)     # relabel.py: #4e4f51


def assemble_v30(panel_bcde_pdf, base_pdf, out_stem):
    """Build the v30 Figure 3 composite from the 9-50 b-e canvas, as ``build_fig3_len950.py`` did.

    Parameters
    ----------
    panel_bcde_pdf : pathlib.Path
        Output of ``panels_bcde`` with ``MIN_LENGTH = 9``.
    base_pdf : pathlib.Path
        The relabelled v28f page, ``expected/Fig3_v28f_relabel_base.pdf``.
    out_stem : pathlib.Path
        Writes ``<stem>.pdf`` and ``<stem>.png`` (300 dpi).

    Returns
    -------
    dict
        ``pdf``, ``png`` paths and the number of re-typeset spans.
    """
    import re
    import fitz

    reg, bold = _FONT_DIR / FONT_F_REGULAR, _FONT_DIR / FONT_F_BOLD
    band = fitz.Rect(0, V30_BAND_ROWS_300DPI[0] * 72 / 300, COMPOSITE_PAGE_PT[0], V30_BAND_ROWS_300DPI[1] * 72 / 300)
    src = fitz.open(panel_bcde_pdf)
    pw, ph = src[0].rect.width, src[0].rect.height
    rect = fitz.Rect(*composite_rects()["bcde"])
    s = min(rect.width / pw, rect.height / ph)
    ox = rect.x0 + (rect.width - pw * s) / 2
    oy = rect.y0 + (rect.height - ph * s) / 2

    def tr(b):
        """Map a span bbox from canvas to page coordinates."""
        return fitz.Rect(ox + b[0] * s, oy + b[1] * s, ox + b[2] * s, oy + b[3] * s)

    doc = fitz.open(base_pdf)
    page = doc[0]
    page.draw_rect(band, color=None, fill=(1, 1, 1), overlay=True)
    page.show_pdf_page(rect, src, 0, clip=None)
    page.insert_font(fontname="L950Reg", fontfile=str(reg))
    page.insert_font(fontname="L950Bold", fontfile=str(bold))
    font, bf = fitz.Font(fontfile=str(reg)), fitz.Font(fontfile=str(bold))
    spans = [dict(sp, direction=ln["dir"]) for b in src[0].get_text("dict")["blocks"]
             for ln in b.get("lines", []) for sp in ln["spans"]]
    n = 0
    for sp in spans:
        text = sp["text"]
        box = tr(sp["bbox"])
        assert band.y0 <= box.y0 and box.y1 <= band.y1, (text, box)
        if text in V30_LEGEND_OLD:
            continue
        letter = text in ["b", "c", "d", "e"]
        numeric = bool(re.fullmatch(r"\d+(?:\.\d+)?", text))
        size = V30_LETTER_SIZE if letter else V30_NUMERIC_SIZE if numeric else V30_TEXT_SIZE
        f, fname = (bf, "L950Bold") if letter else (font, "L950Reg")
        length = f.text_length(text, fontsize=size)
        x, y = ox + sp["origin"][0] * s, oy + sp["origin"][1] * s
        rotate = 0
        if sp["direction"] == (0.0, -1.0):
            rotate = 90
            y = (box.y0 + box.y1) / 2 + length / 2
        elif letter:
            x = box.x0
        elif numeric:
            x = box.x1 - length
        else:
            x = (box.x0 + box.x1) / 2 - length / 2
        page.draw_rect(box + (-.2, -.2, .2, .2), color=None, fill=(1, 1, 1), overlay=True)
        color = tuple(((sp["color"] >> v) & 255) / 255 for v in [16, 8, 0])
        page.insert_text((x, y), text, fontname=fname, fontsize=size, color=color, rotate=rotate)
        n += 1
    page.draw_rect(fitz.Rect(*V30_RELABEL_BOX), color=None, fill=(1, 1, 1), overlay=True)
    page.insert_font(fontname="L950Buch", fontfile=str(reg))
    for y, t in V30_RELABEL:
        page.insert_text((V30_RELABEL_X, y), t, fontname="L950Buch", fontsize=V30_RELABEL_SIZE, color=V30_INK)
    out_stem = Path(out_stem)
    out_stem.parent.mkdir(parents=True, exist_ok=True)
    doc.save(str(out_stem) + ".pdf", garbage=4, deflate=True)
    fitz.open(str(out_stem) + ".pdf")[0].get_pixmap(dpi=300, alpha=False).save(str(out_stem) + ".png")
    return dict(pdf=Path(str(out_stem) + ".pdf"), png=Path(str(out_stem) + ".png"), retyped_spans=n)
