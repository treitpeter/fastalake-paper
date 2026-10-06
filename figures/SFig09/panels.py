"""Supplementary Figure 9: predicted metagenome database size and FastaLake peptide gain.

Panels
------
a  Specimen-level FastaLake gain against predicted-protein records in the
   specimen-matched metagenome (MG) database, 250 CapScan specimens.
b  Median acquisition-level gain in thirds of intestinal-fluid MG database size,
   226 acquisitions (76, 75, 75).

The original producer drew both panels in one matplotlib figure and saved one file
(``Supplementary_Figure_09_proposed.{pdf,svg,png}``). ``figure_ab`` reproduces that call;
``panel_a`` and ``panel_b`` hold the drawing code of each axis, moved unchanged.

Source data (``source_data/SFig09/``)
------------------------------------
capscan_assembly_dna_yield_gain_by_acquisition.tsv
    282 matched CapScan acquisitions: region, participant, BioSample, MG database size,
    recorded DNA quantity, FastaLake and MG peptide counts. Input of the aggregation.
capsule_assembly_gain_per_specimen.tsv
    250 specimens (panel a), written by the aggregation from the table above.
capsule_assembly_gain_acquisitions.tsv
    226 intestinal-fluid acquisitions with their tertile (panel b).
capsule_assembly_gain_tertiles.tsv
    Tertile medians and n (panel b bars and labels).
CAPSCAN_sanitize_report.tsv
    Per-assembly record counts after sanitising (``n_kept``); the MG database sizes are
    checked against it.
capscan_roster_347.tsv
    The v28f CapScan roster (347 acquisitions, 15 participants): region and participant of
    each acquisition; used to cross-check the matched table.
a08_benchmark282_old.tsv, a08_benchmark282_new.tsv
    Per-acquisition FastaLake and MG peptide counts recounted from the peptide sets of the
    production public-sequence arm (old) and of its 347-acquisition rerun (new); used to
    cross-check the counts and as a v28f roster sensitivity.
"""

import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import font_manager

# ---------------------------------------------------------------- choices (named constants)
#: Gain definition, Supplementary Fig. 9 legend (v28d line 445):
#: 100 x (FastaLake peptides / metagenome-search peptides - 1).
GAIN_SCALE = 100
#: Specimen summary in panel a: median over repeated acquisitions (v28d line 445).
SPECIMEN_SUMMARY = "median"
#: Panel b bins: equal-count quantile bins (pandas.qcut) of acquisition-level MG database
#: size over intestinal-fluid acquisitions only; "increasing thirds" (v28d lines 129, 445).
N_BINS = 3
BIN_LABELS = ["Smallest", "Middle", "Largest"]
#: Compartment excluded from panel b and pooled as "stool" on v28d line 129.
STOOL = "stool"
#: Symmetric-log threshold of the gain axis in panel a: "linear up to 50%" (v28d line 445).
LINTHRESH = 50
#: Expected roster sizes (v28d lines 129, 445).
N_MATCHED_ACQUISITIONS = 282
N_SPECIMENS = 250
N_FLUID_ACQUISITIONS = 226
EXPECTED_BIN_N = [76, 75, 75]
#: Regions in plotting order with the producer's colours and labels.
REGION_COLORS = {"duodenum": "#2878a3", "jejunum": "#60a5bf", "ileum": "#8ba777",
                 "colon": "#df8f44", "stool": "#777777"}
REGION_LABELS = dict(zip(REGION_COLORS, ["Region 1", "Region 2", "Region 3", "Region 4", "Stool"]))
#: Figure geometry and output resolution of the producer.
FIGSIZE = (7.2, 3.65)
DPI = 300
#: Font file the producer registered (Söhne Buch, TrueType), read from FASTALAKE_FONT_DIR.
FONT_FILE = "soehne-buch.ttf"
#: Counts the manuscript spells out as words (v28d line 129: "fourteen participants").
NUMBER_WORDS = {14: "fourteen", 15: "fifteen"}


def as_word(n):
    """Return the manuscript's spelled-out form of a small count (e.g. 14 -> 'fourteen')."""
    return NUMBER_WORDS.get(int(n), str(n))


# ---------------------------------------------------------------- loading
def load(sd):
    """Read every source_data table of this unit.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig09``.

    Returns
    -------
    dict of pandas.DataFrame
        Keys ``acq282``, ``spec``, ``caps``, ``summary``, ``census``, ``roster``, ``a08_old``,
        ``a08_new``.
    """
    r = lambda n: pd.read_csv(sd / n, sep="\t")
    return dict(
        acq282=r("capscan_assembly_dna_yield_gain_by_acquisition.tsv"),
        spec=r("capsule_assembly_gain_per_specimen.tsv"),
        caps=r("capsule_assembly_gain_acquisitions.tsv"),
        summary=r("capsule_assembly_gain_tertiles.tsv"),
        census=r("CAPSCAN_sanitize_report.tsv"),
        roster=r("capscan_roster_347.tsv"),
        a08_old=r("a08_benchmark282_old.tsv"),
        a08_new=r("a08_benchmark282_new.tsv"),
    )


# ---------------------------------------------------------------- aggregation (ported)
def derive_tables(acq282):
    """Aggregate the 282 matched acquisitions into the three plotted tables.

    Moved unchanged from the original aggregation step (v15 figure builder). Unit:
    acquisition (282 matched; 226 intestinal fluid) and specimen (250).

    Parameters
    ----------
    acq282 : pandas.DataFrame
        ``capscan_assembly_dna_yield_gain_by_acquisition.tsv``.

    Returns
    -------
    tuple of pandas.DataFrame
        ``(spec, caps, summary)`` as written to the three ``capsule_assembly_gain_*`` tables.
    """
    d = acq282.copy()
    d["gain_exact"] = GAIN_SCALE * (d.FL_peptides / d.MG_peptides - 1)
    spec = d.groupby(["biosample", "subject", "region"], as_index=False).agg(
        assembly_proteins=("assembly_proteins", SPECIMEN_SUMMARY),
        gain_percent=("gain_exact", SPECIMEN_SUMMARY), acquisitions=("sample", "size"))
    caps = d[d.region.ne(STOOL)].copy()
    caps["assembly_tertile"] = pd.qcut(caps.assembly_proteins, N_BINS, labels=BIN_LABELS)
    summary = caps.groupby("assembly_tertile", observed=True).agg(
        median_gain=("gain_exact", "median"), n=("sample", "size"),
        min_assembly=("assembly_proteins", "min"), max_assembly=("assembly_proteins", "max")).reset_index()
    return spec, caps, summary


def check_derived_equal_shipped(t):
    """Compare ``derive_tables`` output with the shipped tables.

    Parameters
    ----------
    t : dict
        Output of ``load``.

    Returns
    -------
    dict
        Maximum absolute difference per table and whether row sets agree.
    """
    spec, caps, summary = derive_tables(t["acq282"])
    out = {}
    a = spec.sort_values("biosample").reset_index(drop=True)
    b = t["spec"].sort_values("biosample").reset_index(drop=True)
    out["spec_rows_equal"] = bool((a.biosample.values == b.biosample.values).all() and len(a) == len(b))
    out["spec_max_abs_diff"] = float(np.nanmax(np.abs(a[["assembly_proteins", "gain_percent", "acquisitions"]].to_numpy(float)
                                                      - b[["assembly_proteins", "gain_percent", "acquisitions"]].to_numpy(float))))
    a = caps.sort_values("sample").reset_index(drop=True)
    b = t["caps"].sort_values("sample").reset_index(drop=True)
    out["caps_rows_equal"] = bool((a["sample"].values == b["sample"].values).all() and len(a) == len(b))
    out["caps_tertile_equal"] = bool((a.assembly_tertile.astype(str).values == b.assembly_tertile.values).all())
    out["caps_gain_max_abs_diff"] = float(np.max(np.abs(a.gain_exact.values - b.gain_exact.values)))
    s = t["summary"]
    out["summary_labels_equal"] = bool((summary.assembly_tertile.astype(str).values == s.assembly_tertile.values).all())
    out["summary_median_max_abs_diff"] = float(np.max(np.abs(summary.median_gain.values - s.median_gain.values)))
    out["summary_n_equal"] = bool((summary.n.values == s.n.values).all())
    out["summary_bounds_equal"] = bool((summary.min_assembly.values == s.min_assembly.values).all()
                                       and (summary.max_assembly.values == s.max_assembly.values).all())
    return out


def check_producer_asserts(t):
    """Run the original producer's input assertions (moved unchanged).

    Checks that every plotted database size equals the sanitizer ``n_kept`` of its
    BioSample, that ``gain_exact`` follows the gain formula, and that tertile medians and
    n reproduce from the acquisition table.

    Parameters
    ----------
    t : dict
        Output of ``load``.

    Returns
    -------
    bool
        True if all assertions pass (an AssertionError is raised otherwise).
    """
    spec, caps, summary = t["spec"], t["caps"], t["summary"]
    census = t["census"].set_index("sample")
    assert census.index.is_unique and spec.biosample.is_unique
    assert len(spec) == N_SPECIMENS and len(caps) == N_FLUID_ACQUISITIONS
    for frame in [spec, caps]:
        assert frame.biosample.isin(census.index).all()
        assert np.array_equal(frame.assembly_proteins, census.loc[frame.biosample, "n_kept"].to_numpy())
    assert np.allclose(caps.gain_exact, GAIN_SCALE * (caps.FL_peptides / caps.MG_peptides - 1), rtol=1e-12, atol=1e-10)
    recomputed = caps.groupby("assembly_tertile").gain_exact.agg(["median", "size"]).loc[summary.assembly_tertile]
    assert np.allclose(recomputed["median"], summary.median_gain, rtol=1e-12)
    assert np.array_equal(recomputed["size"], summary.n)
    assert np.array_equal(summary.n, EXPECTED_BIN_N)
    return True


def check_counts_against_recount(t):
    """Compare FastaLake and MG peptide counts with the independent recount.

    ``a08_benchmark282_old.tsv`` recounts both arms from the per-acquisition peptide sets
    of the production public-sequence arm (5-50 residues, peptide q <= 0.01).

    Parameters
    ----------
    t : dict
        Output of ``load``.

    Returns
    -------
    dict
        Rows joined and numbers of disagreeing FastaLake / MG counts and regions.
    """
    m = t["acq282"].merge(t["a08_old"], on="sample", how="outer", indicator=True)
    both = m[m._merge == "both"]
    return dict(rows_joined=int(len(both)), unmatched=int((m._merge != "both").sum()),
                fl_disagree=int((both.FL_peptides != both.n_fl_5).sum()),
                mg_disagree=int((both.MG_peptides != both.n_mg_5).sum()),
                region_disagree=int((both.region != both.group).sum()))


def check_roster(t):
    """Compare region and participant of the 282 matched acquisitions with the v28f roster.

    Parameters
    ----------
    t : dict
        Output of ``load``.

    Returns
    -------
    dict
        Roster size, participants, matched acquisitions found in the roster, disagreeing
        regions and participants, and roster acquisitions without an MG database by
        participant.
    """
    ro = t["roster"]
    m = t["acq282"].merge(ro, left_on="sample", right_on="sample_id", how="left", indicator=True)
    missing = ro[~ro.sample_id.isin(t["acq282"]["sample"])]
    return dict(roster_acquisitions=int(len(ro)), roster_participants=int(ro.subject_id.nunique()),
                matched_in_roster=int((m._merge == "both").sum()),
                region_disagree=int((m.region_x != m.region_y).sum()),
                participant_disagree=int((m.subject != m.subject_id).sum()),
                unmatched_by_participant={int(k): int(v) for k, v in missing.subject_id.value_counts().sort_index().items()})


# ---------------------------------------------------------------- plotting (moved unchanged)
def _font_family(font_dir):
    """Register the producer's font file and return its family name, or None if absent."""
    if not font_dir:
        return None
    font = Path(font_dir) / FONT_FILE
    if not font.is_file():
        return None
    font_manager.fontManager.addfont(str(font))
    return font_manager.FontProperties(fname=str(font)).get_name()


def panel_a(ax, spec):
    """Draw panel a: specimen-level gain against MG database size.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    spec : pandas.DataFrame
        ``capsule_assembly_gain_per_specimen.tsv`` (250 specimens).
    """
    colors, labels = REGION_COLORS, REGION_LABELS
    for reg, c in colors.items():
        q = spec[spec.region.eq(reg)]
        ax.scatter(q.assembly_proteins, q.gain_percent, s=11, c=c, label=labels[reg], alpha=.8, linewidths=0)
    ax.set(xscale='log', xlabel='Predicted protein records in\nspecimen-matched metagenome database', ylabel='FastaLake peptide gain over MG (%)')
    ax.set_title('a  Specimen-level association', loc='left'); ax.legend(frameon=False, fontsize=7, ncol=2)
    ax.set_ylim(0, spec.gain_percent.max()*1.2); ax.set_yscale('symlog', linthresh=LINTHRESH)
    ax.set_yticks([0,50,100,200,500,1000,2000], [0,50,100,200,500,1000,2000])


def panel_b(ax, summary):
    """Draw panel b: median acquisition-level gain per MG database-size third.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    summary : pandas.DataFrame
        ``capsule_assembly_gain_tertiles.tsv``.
    """
    ax.bar(range(3), summary.median_gain, color='#2878a3', width=.65)
    ax.set_xticks(range(3), summary.assembly_tertile)
    ax.set(xlabel='Capsule metagenome database-size tertile', ylabel='Median acquisition-level gain (%)', ylim=(0,100))
    ax.set_title('b  Capsule acquisitions', loc='left')
    for i, r in summary.iterrows():
        ax.text(i, r.median_gain+2, f'{r.median_gain:.1f}%', ha='center')
        ax.text(i, 4, f'n = {r.n}', ha='center', color='white', fontsize=7)


def figure_ab(spec, summary, font_dir=None):
    """Draw the composite exactly as the original producer did.

    matplotlib defaults are restored first so that the repository style module does not
    alter the producer's rendering; then the producer's rcParams are applied.

    Parameters
    ----------
    spec, summary : pandas.DataFrame
        Panel a and panel b tables.
    font_dir : str or None
        Directory holding ``soehne-buch.ttf`` (``FASTALAKE_FONT_DIR``). Without it text
        falls back to DejaVu Sans and text pixels differ.

    Returns
    -------
    matplotlib.figure.Figure
    """
    plt.rcdefaults()
    family = _font_family(font_dir)
    rc = {'font.size': 8, 'axes.titlesize': 9, 'pdf.fonttype': 42, 'svg.fonttype': 'none'}
    if family:
        rc['font.family'] = family
    plt.rcParams.update(rc)
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE, gridspec_kw={'width_ratios': [1.5, 1]})
    panel_a(axes[0], spec)
    panel_b(axes[1], summary)
    for ax in axes: ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    return fig


def save(fig, stem):
    """Save ``fig`` as ``<stem>.pdf`` and ``<stem>.png`` at the producer's dpi."""
    for ext in ["pdf", "png"]:
        fig.savefig(f"{stem}.{ext}", dpi=DPI)


# ---------------------------------------------------------------- statistics
def stat_a_n_specimens(spec):
    """Number of specimens plotted in panel a.

    Unit: specimen; n = 250 (repeated acquisitions summarised by their median).

    Returns
    -------
    int
    """
    return int(spec.biosample.nunique())


def stat_b_tertile_medians(caps):
    """Median acquisition-level gain per third of MG database size (panel b).

    Unit: acquisition; n = 226 intestinal-fluid acquisitions (76, 75, 75). Descriptive;
    no test, no interval. Repeated acquisitions per specimen and participant are kept.

    Returns
    -------
    pandas.DataFrame
        Columns ``bin``, ``median_gain``, ``n``, ``min_size``, ``max_size``.
    """
    c = caps.copy()
    c["gain"] = GAIN_SCALE * (c.FL_peptides / c.MG_peptides - 1)
    c["bin"] = pd.qcut(c.assembly_proteins, N_BINS, labels=BIN_LABELS)
    return c.groupby("bin", observed=True).agg(median_gain=("gain", "median"), n=("sample", "size"),
                                              min_size=("assembly_proteins", "min"),
                                              max_size=("assembly_proteins", "max")).reset_index()


def stat_129_db_size(acq282):
    """Median MG database size (predicted-protein records) by compartment.

    Unit: acquisition; n = 226 intestinal fluid, 56 stool.

    Returns
    -------
    dict
        ``fluid_median``, ``stool_median``, ``n_fluid``, ``n_stool``.
    """
    fl = acq282.region.ne(STOOL)
    return dict(fluid_median=float(acq282[fl].assembly_proteins.median()),
                stool_median=float(acq282[~fl].assembly_proteins.median()),
                n_fluid=int(fl.sum()), n_stool=int((~fl).sum()))


def stat_129_within_participant(acq282):
    """Participants whose median intestinal-fluid MG database is smaller than stool.

    Unit: participant; medians over that participant's acquisitions in each compartment.
    n = 15 participants, 14 with both compartments.

    Returns
    -------
    dict
        ``n_participants``, ``n_both``, ``n_fluid_smaller``, ``stool_only`` (ids).
    """
    g = (acq282.assign(c=np.where(acq282.region.ne(STOOL), "fluid", "stool"))
         .groupby(["subject", "c"]).assembly_proteins.median().unstack())
    both = g.dropna()
    return dict(n_participants=int(len(g)), n_both=int(len(both)),
                n_fluid_smaller=int((both.fluid < both.stool).sum()),
                stool_only=[int(s) for s in g.index[g.fluid.isna()]])


def stat_129_dna(acq282):
    """Median recorded DNA quantity (ng) by region over acquisitions with a value.

    Unit: acquisition; n = 99 of 282 with a recorded value.

    Returns
    -------
    dict
        ``n_with_dna``, ``n_total``, ``fluid_region_min``, ``fluid_region_max``,
        ``stool``, ``by_region`` (dict).
    """
    d = acq282[acq282.dna_ng.notna()]
    med = d.groupby("region").dna_ng.median()
    fluid = med.drop(STOOL)
    return dict(n_with_dna=int(len(d)), n_total=int(len(acq282)), fluid_region_min=float(fluid.min()),
                fluid_region_max=float(fluid.max()), stool=float(med[STOOL]),
                by_region={k: float(v) for k, v in med.items()})


def stat_129_pooled_yield(acq282):
    """Pooled median intestinal-fluid peptide yield relative to stool, per search arm.

    Ratio of the median peptide count over intestinal-fluid acquisitions (n = 226) to the
    median over stool acquisitions (n = 56). Unit: acquisition.

    Returns
    -------
    dict
        ``mg_percent``, ``fl_percent``.
    """
    fl = acq282.region.ne(STOOL)
    r = lambda col: 100 * acq282[fl][col].median() / acq282[~fl][col].median()
    return dict(mg_percent=float(r("MG_peptides")), fl_percent=float(r("FL_peptides")))


def stat_v28f_roster_sensitivity(acq282, a08_new):
    """Panel b tertile medians with peptide counts from the 347-acquisition rerun.

    The rerun of the public-sequence arm on the 347-acquisition roster changed the
    FastaLake count of some matched acquisitions. MG database sizes and bins are unchanged.
    Unit: acquisition; n = 226.

    Returns
    -------
    tuple
        ``(tertile table, n acquisitions whose FastaLake count changed,
        n specimens whose plotted gain changed)``.
    """
    m = acq282.merge(a08_new[["sample", "n_fl_5", "n_mg_5"]], on="sample")
    assert len(m) == N_MATCHED_ACQUISITIONS
    changed = int(((m.n_fl_5 != m.FL_peptides) | (m.n_mg_5 != m.MG_peptides)).sum())
    new = m.assign(FL_peptides=m.n_fl_5, MG_peptides=m.n_mg_5)
    spec_old, caps_old, _ = derive_tables(acq282)
    spec_new, caps_new, _ = derive_tables(new)
    a = spec_old.sort_values("biosample").gain_percent.values
    b = spec_new.sort_values("biosample").gain_percent.values
    return stat_b_tertile_medians(caps_new), changed, int((~np.isclose(a, b, rtol=0, atol=0)).sum())


# =======================================================================================
# v30: 9-50-residue counts (len950 build.py, build_ylim.py) and the v30 house-style
# restyle (build_sfig09_style.py). Logic moved unchanged.
# =======================================================================================
import io  # noqa: E402
import json  # noqa: E402
import re  # noqa: E402

from matplotlib.collections import PathCollection  # noqa: E402
from matplotlib.patches import Rectangle  # noqa: E402

from style import fl_style as fl  # noqa: E402

# ---- INPUTS (parameterised). The CapScan metagenome benchmark is being extended to 295
#      acquisitions (lab folder final/mg295); SF09 will then be redrawn from new tables.
#      Override by environment, or pass paths to ``render_v30``, without editing this file.
#      Relative values resolve against source_data/SFig09.
ACQ_9_50_TSV = os.environ.get("SFIG09_ACQ_TSV", "capscan_acquisitions_9_50.tsv")   # per-acquisition FL/MG peptides, 9-50 residues
TERTILES_9_50_TSV = os.environ.get("SFIG09_TERTILES_TSV", "capsule_assembly_gain_tertiles_9_50.tsv")  # expected panel b tertile medians/n
REFERENCE_PNG = os.environ.get("SFIG09_BASE_PNG", "")   # image the producer reference render must equal ("" = the unit's expected/SFig09_v29_reference.png)
V29_ACQ_9_50_SHA256 = "5bb9ec898849e853e93ef0d01da9b7eaebab64db2b1a075b89f8d9d44fe5ad39"   # the table the v30 figure was drawn from

LOG_EXP_PT = 5.0                    # exponent of the composed 10^n x tick labels (Söhne, >= 5 pt floor)
V30_W_MM, V30_H_MM = 180.0, 91.2    # placed size in the v29/v30 docx
V30_TICK_PT, V30_LAB_PT, V30_LEG_PT, V30_VAL_PT, V30_N_PT = 6.5, 7.0, 6.5, 6.5, 6.0   # type sizes (pt)
REGION_COLOUR = {"duodenum": fl.SKY, "jejunum": fl.BLUE, "ileum": fl.GREEN, "colon": fl.PURPLE, "stool": fl.GREY}  # Okabe-Ito; orange is reserved for MG
MARKER_S = 11                       # panel a scatter size (pt^2), as the producer


def _sd_path(sd, name):
    p = Path(name)
    return p if p.is_absolute() else Path(sd) / p


def acquisitions_9_50(t):
    """Rebuild the 9-50-residue per-acquisition table from the 282-acquisition table.

    As ``len950/SFig09/build/build.py``: the FastaLake and MG peptide counts of each matched
    acquisition are replaced by the 9-50 counts of ``a08_benchmark282_old.tsv`` (whose 5-50
    columns equal the plotted counts in all 282 rows). Unit: acquisition, n = 282.

    Parameters
    ----------
    t : dict
        Output of ``load``.

    Returns
    -------
    pandas.DataFrame
        Same columns as ``capscan_assembly_dna_yield_gain_by_acquisition.tsv``.
    """
    acq = t['acq282']
    chk = check_counts_against_recount(t)
    assert chk == dict(rows_joined=282, unmatched=0, fl_disagree=0, mg_disagree=0, region_disagree=0), chk
    m = acq.merge(t['a08_old'][['sample', 'n_fl_5', 'n_mg_5', 'n_fl_9', 'n_mg_9']], on='sample', validate='1:1')
    assert len(m) == 282 and (m.n_fl_5 == m.FL_peptides).all() and (m.n_mg_5 == m.MG_peptides).all()
    acq9 = m.assign(FL_peptides=m.n_fl_9, MG_peptides=m.n_mg_9).drop(columns=['n_fl_5', 'n_mg_5', 'n_fl_9', 'n_mg_9'])
    acq9['gain_percent'] = (100 * (acq9.FL_peptides / acq9.MG_peptides - 1)).round(3)
    return acq9


def lowest_ok(ax, fig, spec):
    """Smallest |ymin| step (0.01 % units) so the lowest marker's full disc lies above the x axis."""
    r_px = np.sqrt(MARKER_S) / 2 * fig.dpi / 72
    ymin = 0.0
    lo = spec.gain_percent.min()
    while True:
        ax.set_ylim(ymin, spec.gain_percent.max() * 1.2)
        y_pt = ax.transData.transform((1, lo))[1]
        y_ax = ax.transAxes.transform((0, 0))[1]
        if y_pt - r_px >= y_ax - 1e-9:
            return ymin
        ymin -= 0.01


def figure_v29(spec, summary, font_dir):
    """Producer ``figure_ab`` call for call plus the lowered panel-a y minimum (build_ylim.py).

    Returns
    -------
    tuple
        ``(figure, ymin)``.
    """
    plt.rcdefaults()
    family = _font_family(font_dir)
    rc = {'font.size': 8, 'axes.titlesize': 9, 'pdf.fonttype': 42, 'svg.fonttype': 'none'}
    if family:
        rc['font.family'] = family
    plt.rcParams.update(rc)
    fig, axes = plt.subplots(1, 2, figsize=FIGSIZE, gridspec_kw={'width_ratios': [1.5, 1]})
    panel_a(axes[0], spec)
    ymin = lowest_ok(axes[0], fig, spec) if spec.gain_percent.min() < 0 else 0.0
    axes[0].set_ylim(ymin, spec.gain_percent.max() * 1.2)
    panel_b(axes[1], summary)
    for ax in axes:
        ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    return fig, ymin


def v30_collect(axes):
    """Every plotted quantity of each axes, colour-free, for the before/after check."""
    axes[0].figure.canvas.draw()
    out = []
    for ax in axes:
        out.append(dict(
            xlim=list(ax.get_xlim()), ylim=list(ax.get_ylim()), xscale=ax.get_xscale(), yscale=ax.get_yscale(),
            xticks=[float(v) for v in ax.get_xticks()], yticks=[float(v) for v in ax.get_yticks()],
            xticklabels=[t.get_text() for t in ax.get_xticklabels()], yticklabels=[t.get_text() for t in ax.get_yticklabels()],
            scatter=[np.asarray(c.get_offsets(), float) for c in ax.collections if isinstance(c, PathCollection)],
            bars=[(p.get_x(), p.get_y(), p.get_width(), p.get_height()) for p in ax.patches if isinstance(p, Rectangle)],
            texts=[(tuple(map(float, t.get_position())), t.get_text()) for t in ax.texts],
            ylabel=ax.get_ylabel(), xlabel=ax.get_xlabel(), title=ax.get_title(loc="left"),
            legend=[t.get_text() for t in ax.get_legend().get_texts()] if ax.get_legend() else []))
    return out


def v30_compare(before, after):
    """Exact comparison of producer and restyled plotted arrays; one record per check."""
    rec = []
    for L, b, a in zip("ab", before, after):
        def ok(name, cond, detail=""):
            rec.append(dict(panel=L, check=name, passed=bool(cond), detail=detail))
        ok("scatter arrays identical (count, order, x, y)",
           len(b["scatter"]) == len(a["scatter"]) and all(x.shape == y.shape and np.array_equal(x, y) for x, y in zip(b["scatter"], a["scatter"])),
           f"{len(a['scatter'])} collections, {sum(len(x) for x in a['scatter'])} points")
        ok("bars identical", b["bars"] == a["bars"], f"{len(a['bars'])} bars")
        ok("data-label texts identical (value and n labels)", b["texts"] == a["texts"], str([t for _, t in a["texts"]]))
        ok("x and y limits identical", np.allclose(b["xlim"], a["xlim"], rtol=0, atol=1e-12) and np.allclose(b["ylim"], a["ylim"], rtol=0, atol=1e-12),
           f"x {a['xlim']} y {a['ylim']}")
        ok("scales identical", (b["xscale"], b["yscale"]) == (a["xscale"], a["yscale"]), f"{a['xscale']}/{a['yscale']}")
        ok("x ticks and labels identical", b["xticks"] == a["xticks"] and b["xticklabels"] == a["xticklabels"])
        ok("y ticks and labels identical", b["yticks"] == a["yticks"] and b["yticklabels"] == a["yticklabels"], str(a["yticklabels"]))
        ok("axis labels identical", (b["xlabel"], b["ylabel"]) == (a["xlabel"], a["ylabel"]), f"{a['xlabel']!r} / {a['ylabel']!r}")
        ok("legend entries identical", b["legend"] == a["legend"], str(a["legend"]))
        home = {"a": "legend a 'Each point is one of 250 matched specimens'",
                "b": "x label 'Capsule metagenome database-size tertile' and legend b"}[L]
        ok("producer title removed (content already in v29 axis label/legend)", a["title"] in ("", L), f"{b['title']!r} -> {home}")
    return rec


def v30_panel_a(ax, spec, ymin):
    """Panel a (v30): producer ``panel_a`` calls in the same order; house colours and sizes."""
    for reg, c in REGION_COLOUR.items():
        q = spec[spec.region.eq(reg)]
        ax.scatter(q.assembly_proteins, q.gain_percent, s=MARKER_S, c=c, label=REGION_LABELS[reg], alpha=.8, linewidths=0)
    ax.set(xscale='log', xlabel='Predicted protein records in\nspecimen-matched metagenome database', ylabel='FastaLake peptide gain over MG (%)')
    ax.legend(frameon=False, fontsize=V30_LEG_PT, ncol=2, loc="upper right", handletextpad=0.2, columnspacing=1.0, borderaxespad=0.2)
    ax.set_ylim(0, spec.gain_percent.max() * 1.2)
    ax.set_yscale('symlog', linthresh=LINTHRESH)
    ax.set_yticks([0, 50, 100, 200, 500, 1000, 2000], [0, 50, 100, 200, 500, 1000, 2000])
    ax.set_ylim(ymin, spec.gain_percent.max() * 1.2)


def v30_panel_b(ax, summary):
    """Panel b (v30): median acquisition-level gain per database-size third (226 acquisitions)."""
    ax.bar(range(3), summary.median_gain, color=fl.FL, width=.65)
    ax.set_xticks(range(3), summary.assembly_tertile)
    ax.set(xlabel='Capsule metagenome database-size tertile', ylabel='Median acquisition-level gain (%)', ylim=(0, 100))
    for i, r in summary.iterrows():
        ax.text(i, r.median_gain + 2, f'{r.median_gain:.1f}%', ha='center', fontsize=V30_VAL_PT, color=fl.INK)
        ax.text(i, 4, f'n = {r.n}', ha='center', color='white', fontsize=V30_N_PT)


def log_labels(ax, ref, rec):
    """Draw the producer's 10^n x tick labels as Söhne text: '10' at tick size, exponent at ``LOG_EXP_PT``.

    The tick label strings stay on the axis (hidden) and are checked against the producer.
    """
    from matplotlib.transforms import blended_transform_factory, offset_copy
    ax.tick_params(axis="x", which="both", labelbottom=False)
    tr = blended_transform_factory(ax.transData, ax.transAxes)
    lo, hi = ax.get_xlim()
    drawn = []
    for x, lab in zip(ref["xticks"], ref["xticklabels"]):
        m = re.fullmatch(r"\$\\mathdefault\{10\^\{(-?\d+)\}\}\$", lab)
        if not lab or not (lo <= x <= hi):
            continue
        assert m, lab
        pad = ax.xaxis.get_tick_padding() + ax.xaxis.majorTicks[0].tick1line.get_markersize()
        ax.text(x, 0, "10", transform=offset_copy(tr, fig=ax.figure, x=1.6, y=-pad, units="points"),
                ha="right", va="top", fontsize=V30_TICK_PT, color=fl.INK)
        ax.text(x, 0, m.group(1), transform=offset_copy(tr, fig=ax.figure, x=1.8, y=-pad + 2.6, units="points"),
                ha="left", va="top", fontsize=LOG_EXP_PT, color=fl.INK)
        drawn.append(f"10^{m.group(1)}")
    ax.xaxis.labelpad = ax.xaxis.labelpad + 1.3 * V30_TICK_PT
    rec.append(dict(panel="a", check="log x tick labels composed from the producer strings (10^n), Söhne, base 6.5 pt, exponent 5 pt",
                    passed=len(drawn) > 0, detail=str(drawn)))


def render_v30(sd, out, unit_dir, font_dir, acq_tsv=None, tertiles_tsv=None, reference_png=None):
    """Draw the v30 figure from the 9-50 per-acquisition table and check every plotted value.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig09``.
    out : pathlib.Path
        Output directory: ``SupplementaryFigure9.pdf/.png``, ``values_check.json``.
    unit_dir : pathlib.Path
        ``figures/SFig09`` (for the default reference PNG).
    font_dir : str or None
        ``FASTALAKE_FONT_DIR``.
    acq_tsv, tertiles_tsv, reference_png : str or None
        Override ``ACQ_9_50_TSV``, ``TERTILES_9_50_TSV``, ``REFERENCE_PNG``.

    Returns
    -------
    dict
        ``checks``, ``pdf``, ``png``, ``acq``, ``spec``, ``caps``, ``summary``, ``ymin``,
        ``is_v29_input``.
    """
    import hashlib
    from PIL import Image
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    acq_path = _sd_path(sd, acq_tsv or ACQ_9_50_TSV)
    is_v29 = hashlib.sha256(acq_path.read_bytes()).hexdigest() == V29_ACQ_9_50_SHA256
    ref_png = Path(reference_png or REFERENCE_PNG or (Path(unit_dir) / "expected" / "SFig09_v29_reference.png"))
    acq9 = pd.read_csv(acq_path, sep="\t")
    spec, caps, summary = derive_tables(acq9)
    rec = []
    # 1. reference render == v29 PNG, pixel for pixel (only required for the v29 input table)
    fig0, ymin = figure_v29(spec, summary, font_dir)
    buf = io.BytesIO()
    fig0.savefig(buf, format="png", dpi=DPI)
    buf.seek(0)
    a = np.asarray(Image.open(buf).convert("RGB")).astype(int)
    v = np.asarray(Image.open(ref_png).convert("RGB")).astype(int)
    same = a.shape == v.shape and int((np.abs(a - v).max(axis=2) > 0).sum()) == 0
    rec.append(dict(panel="all", check="reference render equals v29 PNG pixel for pixel" if is_v29 else
                    "reference render vs v29 PNG: not required (new input tables)", passed=same or not is_v29,
                    detail=f"{v.shape[1]}x{v.shape[0]} px, panel-a ymin {ymin:.2f}"))
    before = v30_collect(fig0.axes)
    plt.close(fig0)
    # 2. restyled
    fl.init()
    fig, axes = plt.subplots(1, 2, figsize=fl.mm(V30_W_MM, V30_H_MM), gridspec_kw={'width_ratios': [1.5, 1]})
    fig.subplots_adjust(left=0.075, right=0.99, bottom=0.16, top=0.95, wspace=0.24)
    v30_panel_a(axes[0], spec, ymin)
    v30_panel_b(axes[1], summary)
    for ax, L in zip(axes, "ab"):
        fl.style(ax)
        ax.tick_params(labelsize=V30_TICK_PT)
        ax.xaxis.label.set_size(V30_LAB_PT)
        ax.yaxis.label.set_size(V30_LAB_PT)
        fl.panel(ax, L)
    after = v30_collect(axes)
    log_labels(axes[0], before[0], rec)
    rec += v30_compare(before, after)
    exp = pd.read_csv(_sd_path(sd, tertiles_tsv or TERTILES_9_50_TSV), sep="\t")
    rec.append(dict(panel="b", check="tertile medians and n equal the expected tertile table (9-50)",
                    passed=[round(x, 9) for x in exp.median_gain] == [round(x, 9) for x in summary.median_gain]
                    and list(exp.n) == list(summary.n), detail=str([(round(x, 1), int(n)) for x, n in zip(summary.median_gain, summary.n)])))
    rec.append(dict(panel="a", check="one point per specimen (250 for the v29 inputs)",
                    passed=sum(len(x) for x in after[0]["scatter"]) == len(spec) and (len(spec) == N_SPECIMENS or not is_v29),
                    detail=f"{len(spec)} specimens"))
    (out / "values_check.json").write_text(json.dumps(rec, indent=1) + "\n")
    bad = [r for r in rec if not r["passed"]]
    if bad:
        raise SystemExit("VALUE CHECK FAILED: " + json.dumps(bad, indent=1))
    pdf, png = fl.save(fig, out / "SupplementaryFigure9", formats=("pdf", "png"))
    plt.close(fig)
    return dict(checks=rec, pdf=Path(pdf), png=Path(png), acq=acq9, spec=spec, caps=caps, summary=summary,
                ymin=ymin, is_v29_input=is_v29)
