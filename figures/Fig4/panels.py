"""Figure 4: public-sequence benchmarking and PREDICT molecular profiles.

Panels and their source data (all under ``source_data/Fig4/``):

* a  schematic, shipped as ``artwork/panel_a.pdf`` (see ARTWORK_NUMBERS.md)
* b,c ``crossfit_per_acquisition.tsv``           10-fold cross-fitted peptides and recovery, 9-50 residues
* d  ``spectrum_audit_per_acquisition.tsv``       MG-search outcome of FL-only spectra, PREDICT
* e  controlled ``paired_antibiotic_profiles.tsv`` + ``clinical_paired_statistics.tsv``
* f  ``clinical_composition_means.tsv``           mean origin composition by antibiotic state (aggregate)
* g  controlled ``paired_antibiotic_profiles.tsv`` + ``clinical_paired_statistics.tsv``
* h  ``heatmap_features.tsv``, ``heatmap_z.tsv``, ``heatmap_tracks.tsv`` + controlled
     ``heatmap_85_patient_order.tsv`` (antibiotic track only)
* i  ``crossomics_native_correlations.tsv``, ``crossomics_functional_correlations.tsv``;
     the correlations are recomputed from ``crossomics_distances_predict77.tsv``

Text statistics also use ``same_spectrum_benchmark.tsv``, ``predict_acquisition_roster.tsv``
and ``capscan_dna_yield_by_acquisition.tsv``.

Every plotting function is the original producer's code moved here with paths made
relative and data read from ``source_data``; the producer is named in each docstring
(LAB_ROOT locations in PROVENANCE.md). Each function starts from matplotlib defaults and
applies the producer's own settings, because the producers ran as separate processes.
"""

import itertools
import os
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import cbook, font_manager
from matplotlib.colors import LinearSegmentedColormap, ListedColormap, Normalize
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle

# ----------------------------------------------------------------------------------
# Choices (README "Choices" table). Methods line numbers refer to v28d_accepted.txt.
# ----------------------------------------------------------------------------------
# Peptide length windows (min, max residues). The source tables hold counts for minimum
# lengths 5 and 9 (maximum 50 in both), in columns suffixed _5 and _9.
LENGTH_WINDOW_FIG4BC = (9, 50)      # cross-fitted b,c; Results v28d line 123 (legend line 127 says 5-50; changes in v28f)
LENGTH_WINDOW_SAME_SPECTRUM = (5, 50)  # same-spectrum benchmark; Results v28d line 123, Methods v28d line 267
L_BC = str(LENGTH_WINDOW_FIG4BC[0])            # column suffix for the b,c window
L_SAME = str(LENGTH_WINDOW_SAME_SPECTRUM[0])   # column suffix for the 5-50 window
LENGTHS = (L_SAME, L_BC)            # order in which the bootstrap evaluates the two windows
BOOT_DRAWS_CROSSFIT = 10_000        # participant bootstrap resamples; Methods v28d lines 267 and 777
BOOT_SEED_CROSSFIT = 20261004       # seed of that bootstrap; Methods v28d line 777
CROSSFIT_STRATA = ["PREDICT", "CAPSCAN_all", "CAPSCAN_stool", "CAPSCAN_intestinal"]  # stream order; not stated in the text (v28d line 777 names a separate region stream)
CAPSCAN_REGIONS = ["duodenum", "jejunum", "ileum", "colon", "stool"]  # region stream order; v28d line 777 ("separate bootstrap stream with the same seed")
REGION_LABEL = {"duodenum": "R1", "jejunum": "R2", "ileum": "R3", "colon": "R4", "stool": "stool"}  # regions 1-4, Results v28d line 129
CI_PERCENTILES = (2.5, 97.5)        # percentile bootstrap 95% interval; Methods v28d lines 267, 275, 777
BOOT_SENSITIVITY_SEEDS = tuple(range(20261005, 20261025))  # 20 further seeds for the stool-CI Monte Carlo sensitivity (repository addition, not in the text)
CLIN_BOOT_DRAWS = 50_000            # patient bootstrap resamples; Methods v28d line 275
CLIN_SIGNFLIP_DRAWS = 200_000       # sign assignments (100 x 2,000); Methods v28d line 275
CLIN_SEED = 20260907                # bootstrap seed; sign-flip seed is CLIN_SEED + endpoint index; not stated in the text (OPEN_ITEMS TEXT-GAP)
MANTEL_PERMUTATIONS = 1_999         # patient-label permutations; Methods v28d line 815
MANTEL_SEED = 20260907              # seed of those permutations; Methods v28d line 815
Q_STAR = 0.05                       # single star in Fig. 4i; legend v28d line 127 defines three tiers (OPEN_ITEMS AUTHOR)
BOX_WHIS = 1.5                      # whisker length in IQR, Fig. 4e,g; legend v28d line 127
HEATMAP_CLIP = 2.0                  # z scores clipped to +-2 for display, Fig. 4h; not stated (legend v28d line 127 says z scores)
MIXED_MODEL_LENGTH = L_BC           # mixed model on 9-50-residue cross-fitted counts; v30 Results paragraph 69 (v28d line 129 printed the 5-50 fit), Methods v28d line 229
SAME_VS_CROSSFIT_LENGTH = L_BC      # same-spectrum vs cross-fitted gain comparison at 9-50 residues; v30 Results paragraph 70 (v28d line 133 used 5-50)
TEXT_LENGTH = int(L_BC)             # same-spectrum text values at 9-50 residues; v30 Results paragraphs 64, 65 and 70 (v28d lines 123 and 129 used 5-50)
H_EDS_TRACK_LABEL = 'Enterotype dysbiosis score (0–1)'  # Fig. 4h top track label in v30 (fig4_eds_label/relabel.py; producer fig_v27 printed 'EDS (0–1)')

UNIT = Path(__file__).resolve().parent
REPO = UNIT.parents[1]
SD = REPO / "source_data" / UNIT.name

#: Page geometry of the v28f composite (points); used by b,c and the panel-i colour bar.
PAGE_W = 518.7401733398438
PAGE_H_BC = 688.0                     # fig4bc_crossfit build.py
PAGE_H_V27 = 687.9999996757813        # fig_v27 page height after crop (COMPOSE_fig4.json)
I_DXI = 66.89008666992186             # panel-i body shift measured by fig_v27 (ko_block_right - 382)


# ----------------------------------------------------------------------------------
# fonts and controlled inputs
# ----------------------------------------------------------------------------------
def font_dir():
    """Return the Söhne font directory.

    Parameters
    ----------
    None. Reads the environment variable FASTALAKE_FONT_DIR (alias FL_STYLE_DIR).

    Returns
    -------
    pathlib.Path or None
        The directory, or None when it is not set or does not exist.
    """
    d = os.environ.get("FASTALAKE_FONT_DIR") or os.environ.get("FL_STYLE_DIR")
    return Path(d) if d and Path(d).is_dir() else None


def fonts_present():
    """Report whether the Söhne fonts needed for pixel-exact text are available.

    Parameters
    ----------
    None.

    Returns
    -------
    bool
        True when soehne-buch.ttf and soehne-buch.otf are both in :func:`font_dir`.
    """
    d = font_dir()
    return d is not None and (d / "soehne-buch.ttf").is_file() and (d / "soehne-buch.otf").is_file()


def _buch_family():
    """Register soehne-buch.ttf as the producers did and return its family name.

    Returns
    -------
    str
        The family name of Söhne Buch, or ``"DejaVu Sans"`` when the font is absent.
    """
    d = font_dir()
    if d is None or not (d / "soehne-buch.ttf").is_file():
        return "DejaVu Sans"
    font = d / "soehne-buch.ttf"
    font_manager.fontManager.addfont(str(font))
    return font_manager.FontProperties(fname=str(font)).get_name()


def controlled_path(name):
    """Return the path of a controlled input, or None (printing why) when it is absent.

    Parameters
    ----------
    name : str
        File name listed in controlled/CONTROLLED_INPUTS.md.

    Returns
    -------
    pathlib.Path or None
        Path under FASTALAKE_CONTROLLED_DIR, or None when the variable or file is absent.
    """
    root = os.environ.get("FASTALAKE_CONTROLLED_DIR")
    p = Path(root) / name if root else None
    if p is None or not p.is_file():
        print(f"SKIPPED: needs controlled input {name} (see controlled/CONTROLLED_INPUTS.md)")
        return None
    return p


def _rd(name, **kw):
    """Read a source_data table with round-trip float parsing.

    Floats were written with 17 significant digits, so round-trip parsing returns exactly
    the values the original producers held.

    Parameters
    ----------
    name : str
        File name in ``source_data/Fig4/``.
    **kw
        Passed to :func:`pandas.read_csv`.

    Returns
    -------
    pandas.DataFrame
    """
    return pd.read_csv(SD / name, sep="\t", float_precision="round_trip", **kw)


# ----------------------------------------------------------------------------------
# loaders
# ----------------------------------------------------------------------------------
#: Tables copied byte for byte from the upstream analysis; read with pandas' default
#: float parser, as the original producers read them. All other tables were rewritten
#: with 17 significant digits and are read with round-trip parsing.
COPIED_VERBATIM = {"clinical_paired_statistics.tsv", "crossomics_native_correlations.tsv",
                   "crossomics_functional_correlations.tsv", "capscan_dna_yield_by_acquisition.tsv"}


def load_table(name):
    """Read ``source_data/Fig4/<name>`` with the parser that reproduces the producer's values.

    Parameters
    ----------
    name : str
        File name in ``source_data/Fig4/``.

    Returns
    -------
    pandas.DataFrame
    """
    if name in COPIED_VERBATIM:
        return pd.read_csv(SD / name, sep="\t")
    return _rd(name)


def load_crossfit():
    """Cross-fit table with the stratum used by the producer.

    Parameters
    ----------
    None. Reads ``source_data/Fig4/crossfit_per_acquisition.tsv``.

    Returns
    -------
    pandas.DataFrame
        578 acquisitions (296 PREDICT, 282 CapScan) in the original row order.
    """
    x = _rd("crossfit_per_acquisition.tsv", dtype={"participant_id": str}, keep_default_na=False, na_values=[""])
    x["lda_fallback"] = x["lda_fallback"].fillna("")
    assert len(x) == 591 and list(x.crossfit_index) == list(range(591))
    x["stratum"] = np.where(x.cohort.eq("PREDICT"), "PREDICT",
                            np.where(x.sampling_site.eq("stool"), "CAPSCAN_stool", "CAPSCAN_intestinal"))
    return x


def load_paired_controlled():
    """Read the controlled paired antibiotic profiles.

    Parameters
    ----------
    None. Reads ``paired_antibiotic_profiles.tsv`` from FASTALAKE_CONTROLLED_DIR.

    Returns
    -------
    pandas.DataFrame or None
        84 rows (42 patients x off/on), or None when the controlled input is absent.
    """
    p = controlled_path("paired_antibiotic_profiles.tsv")
    return None if p is None else pd.read_csv(p, sep="\t")


def load_heatmap():
    """Read the panel-h inputs.

    Parameters
    ----------
    None. Reads ``heatmap_features.tsv``, ``heatmap_z.tsv`` and ``heatmap_tracks.tsv``.

    Returns
    -------
    tuple
        ``(Z, feat, tracks)``: z array (20 features x 85 patients, display order), the
        feature table and the per-patient track table (85 patients).
    """
    feat = _rd("heatmap_features.tsv")
    zl = _rd("heatmap_z.tsv")
    Z = zl.pivot(index="feature_id", columns="display_position", values="z").loc[feat.feature_id].to_numpy()
    tracks = _rd("heatmap_tracks.tsv")
    assert Z.shape == (20, 85) and list(tracks.display_position) == list(range(1, 86))
    return Z, feat, tracks


def load_antibiotic_track():
    """Read the antibiotic track of panel h from the controlled file.

    Parameters
    ----------
    None. Reads ``heatmap_85_patient_order.tsv`` from FASTALAKE_CONTROLLED_DIR and aligns
    it to ``heatmap_tracks.tsv`` by EDS score.

    Returns
    -------
    numpy.ndarray or None
        Antibiotic state per heatmap column (85 patients, display order), or None when
        the controlled input is absent.
    """
    p = controlled_path("heatmap_85_patient_order.tsv")
    if p is None:
        return None
    o = pd.read_csv(p, sep="\t")
    tr = _rd("heatmap_tracks.tsv")
    assert len(o) == 85 and np.array_equal(o.dysbiosis_score.to_numpy(), tr.eds_score.to_numpy())
    return o.antibiotic_state.to_numpy()


# ----------------------------------------------------------------------------------
# Panels b,c
# ----------------------------------------------------------------------------------
def panel_bc(x, out_pdf):
    """Fig. 4b,c: cross-fitted peptides (b) and MG-identity recovery (c), 9-50 residues.

    Producer: fig4bc_crossfit/build/build.py, section "Fig. 4b,c (new)". Draws both panels
    on a transparent page of the composite's size (518.74 x 688 pt) at their page
    positions, as the producer did.

    Parameters
    ----------
    x : pandas.DataFrame
        Output of :func:`load_crossfit` (acquisitions; 296 PREDICT, 56 CapScan stool,
        226 CapScan intestinal fluid).
    out_pdf : pathlib.Path
        PDF to write.

    Returns
    -------
    pathlib.Path
        ``out_pdf``.
    """
    from style import fl_style as fl
    BUCH = _buch_family()
    W, H = PAGE_W, PAGE_H_BC
    LG = '#cfcfcf'
    STRATA = [('PREDICT', 'PREDICT'), ('CAPSCAN_stool', 'CapScan stool'), ('CAPSCAN_intestinal', 'Intestinal fluid')]
    x = x.rename(columns={f"recovery_cf_{L_BC}_pct": "rec_cf_9"})
    with mpl.rc_context():
        mpl.rcdefaults()
        plt.rcParams.update({'font.family': BUCH, 'font.size': 6.0, 'axes.labelsize': 6.0, 'xtick.labelsize': 6.0, 'ytick.labelsize': 6.0,
                             'axes.linewidth': .46, 'xtick.major.width': .46, 'ytick.major.width': .46, 'xtick.major.size': 3.2, 'ytick.major.size': 3.2,
                             'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'savefig.bbox': None,
                             'text.color': fl.INK, 'axes.labelcolor': fl.INK, 'axes.edgecolor': fl.INK, 'xtick.color': fl.INK, 'ytick.color': fl.INK,
                             'axes.unicode_minus': True})
        fig = plt.figure(figsize=(W / 72, H / 72))

        def ax_pt(x0, y0, x1, y1):          # page points, y from top
            """Add axes spanning page points (x0, y0)-(x1, y1), y measured from the top; return them."""
            return fig.add_axes([x0 / W, 1 - y1 / H, (x1 - x0) / W, (y1 - y0) / H])
        Y0, Y1 = 88.1, 171.0
        axB = ax_pt(64.8, Y0, 196.0, Y1); axC = ax_pt(232.8, Y0, 318.0, Y1)
        rng = np.random.default_rng(4)
        POS = {'PREDICT': (0, 1), 'CAPSCAN_stool': (2.3, 3.3), 'CAPSCAN_intestinal': (4.6, 5.6)}
        YMAX = 55
        for s, lab in STRATA:
            g = x[x.stratum.eq(s)]; x0, x1 = POS[s]
            j0 = rng.uniform(-.14, .14, len(g)); j1 = rng.uniform(-.14, .14, len(g))
            mg, flv = g[f'MG_pooled_{L_BC}'].to_numpy() / 1e3, g[f'FL_pooled_{L_BC}'].to_numpy() / 1e3
            assert max(mg.max(), flv.max()) < 50
            for a, b, ja, jb in zip(mg, flv, j0, j1):
                axB.plot([x0 + ja, x1 + jb], [a, b], color=LG, lw=.28, alpha=.6, zorder=1)
            axB.scatter(x0 + j0, mg, s=2.5, color=fl.METAG, alpha=.7, lw=0, zorder=2)
            axB.scatter(x1 + j1, flv, s=2.5, color=fl.FL, alpha=.7, lw=0, zorder=2)
            for xx, v in ((x0, np.median(mg)), (x1, np.median(flv))):
                axB.plot([xx - .3, xx + .3], [v, v], color='k', lw=1.0, zorder=3)
            # label: median paired change 100 x (FL/MG - 1), the statistic the text reports
            # (per-arm medians can point the other way, e.g. PREDICT 18.8k vs 19.0k at -3.5%)
            assert (mg > 0).all()
            d = float(np.median(100 * (flv / mg - 1)))
            axB.text((x0 + x1) / 2, 53.0, f'{d:+.1f}%'.replace('-', '−'), ha='center', va='center', fontsize=5.3, zorder=4)
        axB.set_xticks([p for s, _ in STRATA for p in POS[s]], ['MG', 'FL'] * 3)
        axB.set_ylabel('Peptides per acquisition\n(thousands; cross-fitted,\n9–50 residues)', labelpad=2)
        axB.set_ylim(0, YMAX); axB.set_yticks([0, 10, 20, 30, 40, 50]); axB.set_xlim(-.55, 6.15)
        GL = {'PREDICT': 'PREDICT\n(n = 296)', 'CAPSCAN_stool': 'CapScan stool\n(n = 56)', 'CAPSCAN_intestinal': 'Intestinal fluid\n(n = 239)'}
        for s, _ in STRATA:
            axB.text(np.mean(POS[s]), -.155, GL[s], ha='center', va='top', transform=axB.get_xaxis_transform(), fontsize=6.0, linespacing=1.1)
        for i, (s, lab) in enumerate(STRATA):
            g = x[x.stratum.eq(s)]; v = g.rec_cf_9.dropna().to_numpy()
            j = rng.uniform(-.15, .15, len(v))
            axC.scatter(i + j, v, s=2.5, color=fl.FL, alpha=.7, lw=0)
            m = float(np.median(v)); axC.plot([i - .27, i + .27], [m, m], color='k', lw=1.0, zorder=3)
            axC.text(i, 100.0, f'{m:.1f}%', ha='center', va='center', fontsize=5.3, zorder=4)
            assert v.max() < 95
            assert v.min() > 0 and v.max() <= 100
        axC.set_xticks([0, 1, 2], ['PREDICT', 'CapScan\nstool', 'Intestinal\nfluid'])
        axC.set_ylabel('MG peptides found by FL (%)\ncross-fitted, 9–50 residues', labelpad=2)
        axC.set_ylim(0, 105); axC.set_yticks([0, 20, 40, 60, 80, 100]); axC.set_xlim(-.45, 2.45)
        for ax in (axB, axC):
            ax.tick_params(axis='both', pad=2)
        fig.savefig(out_pdf, transparent=True)
        plt.close(fig)
    return out_pdf


# ----------------------------------------------------------------------------------
# Panel d
# ----------------------------------------------------------------------------------
AUDIT_CATEGORIES = [
    ('D_displaced', 'Different MG peptide\nnot accepted'),
    ('A_absent', 'FL sequence absent\nfrom MG database'),
    ('C_top_fails_q', 'Same MG peptide\nnot accepted'),
    ('B_not_cleavable', 'FL sequence present,\noutside audited digest'),
    ('E_no_target', 'No MG target\ncandidate reported')]


def stat_d_audit(df):
    """Partition of PREDICT spectra accepted only by FastaLake into five MG-search outcomes.

    Unit: spectrum (pooled over 296 acquisitions). No test; percentages of the total.

    Parameters
    ----------
    df : pandas.DataFrame
        ``spectrum_audit_per_acquisition.tsv``.

    Returns
    -------
    dict
        ``total`` (spectra), ``counts`` and ``percent`` per category.
    """
    assert len(df) == 296
    assert (df[[k for k, _ in AUDIT_CATEGORIES]].sum(axis=1) == df.total).all()
    assert (df.n_fl_spectra_q01 - df.n_shared_spectra_q01 == df.total).all()
    counts = {k: int(df[k].sum()) for k, _ in AUDIT_CATEGORIES}
    tot = sum(counts.values())
    return dict(total=tot, counts=counts, percent={k: 100 * v / tot for k, v in counts.items()})


def _legacy_fl_style_init():
    """fl.init() of the 2026-08 style module the panel-d producer imported.

    That module (LAB_ROOT/figure_style/fl_style.py, sha256 58b64c34...) differs from the
    current style/fl_style.py: 0.6-pt lines from the lab styleguide, no 1-pt floor.
    Reproduced here so panel d is drawn exactly as it was.
    """
    d = font_dir()
    font_path = str(d / 'soehne-buch.otf') if d else ''
    mpl.rcParams['figure.dpi'] = 300
    mpl.rcParams['pdf.fonttype'] = 42
    mpl.rcParams['ps.fonttype'] = 42
    mpl.rcParams['svg.fonttype'] = 'none'
    if font_path and os.path.isfile(font_path):
        fe = font_manager.FontEntry(fname=font_path, name='soehne')
        font_manager.fontManager.ttflist.insert(0, fe)
        mpl.rcParams['font.family'] = fe.name
    else:
        mpl.rcParams['font.family'] = 'sans-serif'
    for k in ('text.color', 'axes.labelcolor', 'axes.edgecolor', 'xtick.color', 'ytick.color'):
        mpl.rcParams[k] = '#4e4f51'
    for k in ('lines.linewidth', 'axes.linewidth', 'xtick.major.width', 'ytick.major.width',
              'xtick.minor.width', 'ytick.minor.width'):
        mpl.rcParams[k] = 0.6
    mpl.rcParams.update({'font.size': 7, 'axes.titleweight': 'normal',
                         'axes.titlesize': 7, 'legend.fontsize': 6.5,
                         'xtick.major.size': 2.6, 'ytick.major.size': 2.6,
                         'figure.dpi': 300, 'savefig.dpi': 300,
                         'savefig.bbox': 'tight', 'savefig.pad_inches': 0.02,
                         'svg.fonttype': 'none'})
    mpl.rcParams.update({'mathtext.fontset': 'custom',
                         'mathtext.rm': 'soehne', 'mathtext.it': 'soehne',
                         'mathtext.bf': 'soehne', 'mathtext.sf': 'soehne',
                         'mathtext.tt': 'soehne', 'mathtext.cal': 'soehne',
                         'mathtext.default': 'regular'})


def _legacy_style(ax):
    """fl.style(ax) of the same 2026-08 module."""
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.tick_params(length=2.6, width=0.6, labelsize=7)
    ax.xaxis.label.set_size(7)
    ax.yaxis.label.set_size(7)
    return ax


def panel_d(df, out_pdf):
    """Fig. 4d: MG-search outcome of 6.69 million PREDICT spectra accepted only by FL.

    Producer: cleanup_2026-09-22/figure_updates_2026-09-22/build_saved_figures.py (lines
    21-52), output Figure4d.pdf. Two strings it drew ("MG search outcome", "Accepted:
    peptide q at most 0.01") were later clipped out of the composite (ASSEMBLY.md).

    Parameters
    ----------
    df : pandas.DataFrame
        ``spectrum_audit_per_acquisition.tsv`` (296 acquisitions).
    out_pdf : pathlib.Path
        PDF to write.

    Returns
    -------
    pathlib.Path
        ``out_pdf``.
    """
    INK, FL = '#4e4f51', '#0072B2'
    a = stat_d_audit(df)
    counts = a['counts']
    assert a['total'] == 6689188
    assert counts == {'D_displaced': 2925662, 'A_absent': 1939387, 'C_top_fails_q': 1663345,
                      'B_not_cleavable': 90012, 'E_no_target': 70782}
    percent = [100 * counts[k] / 6689188 for k, _ in AUDIT_CATEGORIES]
    with mpl.rc_context():
        mpl.rcdefaults()
        _legacy_fl_style_init(); plt.rcParams['savefig.bbox'] = None
        fw, fh = 69, 45
        fig = plt.figure(figsize=(fw / 25.4, fh / 25.4))
        ax = fig.add_axes([25 / fw, 12 / fh, 40 / fw, 27 / fh])
        _legacy_style(ax)
        y = np.arange(5)[::-1]
        ax.barh(y, percent, color=FL, height=.64)
        ax.set_yticks(y, [label for _, label in AUDIT_CATEGORIES], fontsize=6.1)
        ax.tick_params(axis='y', length=0, pad=3)
        ax.set_xlim(0, 55); ax.set_xticks([0, 20, 40], ['0', '20', '40'], fontsize=6)
        ax.set_ylim(-.6, 4.6)
        for yy, pct in zip(y, percent):
            ax.text(pct + 1.2, yy, f'{pct:.1f}%', va='center', fontsize=6, color=INK)
        fig.text(25 / fw, 1 - 2.7 / fh, 'd', fontsize=7.2, va='center', color=INK)
        fig.text(35 / fw, 1 - 2.7 / fh, 'MG search outcome', fontsize=6.5, va='center', color=INK)
        fig.text(45 / fw, 1 - 39.2 / fh, 'FL-only accepted spectra (%)', fontsize=6.1, ha='center', va='center', color=INK)
        fig.text(.52, 1 - 43.2 / fh, 'Accepted: peptide q at most 0.01', fontsize=5.8, ha='center', va='center', color=INK)
        fig.savefig(out_pdf)
        plt.close(fig)
    return out_pdf


# ----------------------------------------------------------------------------------
# Panels e, f, g (antibiotic exposure, 42 patients)
# ----------------------------------------------------------------------------------
def _d12_rc(BUCH):
    """Apply the rcParams of the D12 producer (build_figures.py) with font family ``BUCH``."""
    plt.rcParams.update({'font.family': BUCH, 'font.size': 6.4, 'axes.labelsize': 6.4, 'xtick.labelsize': 6, 'ytick.labelsize': 6,
                         'axes.linewidth': .5, 'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'savefig.bbox': None})


def panel_e(ab, st, out_pdf):
    """Fig. 4e: human-assigned LFQ share off and on antibiotics, 42 patients.

    Producer: iteration4/figures/shadow_d12/B/build_figures.py (D12), first block.

    Parameters
    ----------
    ab : pandas.DataFrame
        Controlled ``paired_antibiotic_profiles.tsv`` (84 rows = 42 patients x 2 states).
    st : pandas.DataFrame
        ``clinical_paired_statistics.tsv`` indexed by metric.
    out_pdf : pathlib.Path
        PDF to write.

    Returns
    -------
    pathlib.Path
        ``out_pdf``.
    """
    blue = '#0072B2'; grey = '#cccccc'
    with mpl.rc_context():
        mpl.rcdefaults(); _d12_rc(_buch_family())
        w = ab.pivot(index='patient_id', columns='antibiotic_state', values='Human')[['Off', 'On']] * 100; s = st.loc['Human']
        assert len(w) == 42 and np.isclose((w.On - w.Off).mean(), s.mean_on_minus_off)
        fig = plt.figure(figsize=(136 / 72, 120 / 72)); ax = fig.add_axes([31.85 / 136, 17.03 / 120, 98.54 / 136, 88.19 / 120])
        for a, b in w.itertuples(index=False, name=None):
            ax.plot([0, 1], [a, b], color=grey, lw=.4, zorder=1)
        for x, c in [(0, '#888888'), (1, '#D55E00')]:
            ax.scatter([x] * 42, w.iloc[:, x], s=4.2, color=c, lw=0, zorder=3)
        ax.boxplot([w.Off, w.On], positions=[0, 1], widths=.46, whis=BOX_WHIS, showfliers=False, boxprops={'linewidth': .65},
                   whiskerprops={'linewidth': .65}, capprops={'linewidth': .65}, medianprops={'color': 'black', 'linewidth': 1.1}, zorder=4)
        ax.set(xlim=(-.5, 1.5), ylim=(0, 110), xticks=[0, 1], xticklabels=['Antibiotics\noff', 'Antibiotics\non'], yticks=[0, 25, 50, 75, 100])
        ax.set_ylabel('Human-assigned LFQ signal (%)', labelpad=2); ax.set_title('e', loc='left', fontsize=8); ax.tick_params(axis='x', length=0, pad=1.5)
        ax.plot([0, 0, 1, 1], np.array([.93, .96, .96, .93]) * 110, color='#333333', lw=.6); ax.text(.5, 106.7, '***', ha='center', va='bottom')
        ax.text(.38, 94.6, f'+{s.mean_on_minus_off:.1f} points\n95% CI {s.CI95_low:.1f}–{s.CI95_high:.1f}\nn = 42 patients', ha='center', va='top', fontsize=5.3, linespacing=1.05)
        fig.savefig(out_pdf); plt.close(fig)
    del blue
    return out_pdf


def panel_f(mean_table, out_pdf, with_labels=True):
    """Fig. 4f: mean origin composition of LFQ signal by antibiotic state (aggregate of 42 patients).

    Producer: final/fig_v27/build_fig4_v27.py, ``build_f`` (D12 panel f plus in-bar
    percentages for segments of at least 10%).

    Parameters
    ----------
    mean_table : pandas.DataFrame
        ``clinical_composition_means.tsv`` (rows Off, On; column ``condition``).
    out_pdf : pathlib.Path
    with_labels : bool
        False reproduces the unlabelled D12 panel. ``out_pdf`` is the PDF to write.

    Returns
    -------
    list of dict
        Printed in-bar labels.
    """
    with mpl.rc_context():
        mpl.rcdefaults()
        plt.rcParams.update({'font.family': _buch_family(), 'font.size': 6.4,
                             'axes.labelsize': 6.4, 'xtick.labelsize': 6, 'ytick.labelsize': 6, 'axes.linewidth': .5,
                             'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'savefig.bbox': None})
        pink = '#CC79A7'; blue = '#0072B2'; orange = '#E69F00'; grey = '#cccccc'
        classes = ['Human', 'Microbial_species', 'Microbial_above_species', 'Other_annotated_origins', 'Unresolved_origin']
        labels = ['Human assigned', 'Microbial: species', 'Microbial: above species', 'Other annotated origins', 'Unresolved origin']
        colors = [pink, blue, '#56B4E9', orange, grey]
        mean = mean_table.set_index('condition').loc[['Off', 'On']]
        assert np.allclose(mean[classes].sum(1), 100)
        fig = plt.figure(figsize=(150 / 72, 128 / 72)); ax = fig.add_axes([.2, .13, .33, .68]); bottom = np.zeros(2)
        printed = []
        for k, label, color in zip(classes, labels, colors):
            ax.bar(range(2), mean[k], bottom=bottom, color=color, width=.6, label=label)
            if with_labels:
                for x in range(2):
                    v = mean[k].iloc[x]
                    if v >= 10:
                        ax.text(x, bottom[x] + v / 2, f'{v:.0f}%', ha='center', va='center', fontsize=5.0, color='black')
                        printed.append(dict(condition=mean.index[x], origin_class=k, mean_LFQ_signal_pct=v, printed_label=f'{v:.0f}%'))
            bottom += mean[k].to_numpy()
        ax.set(ylim=(0, 100), yticks=[0, 50, 100], xticks=[0, 1], xticklabels=['Off', 'On']); ax.set_title('f', loc='left', fontsize=8)
        ax.set_ylabel('LFQ signal (%)', labelpad=1); ax.tick_params(axis='x', length=0)
        ax.legend(loc='upper left', bbox_to_anchor=(1.01, 1.02), fontsize=5.0, handlelength=.8, handletextpad=.4, borderaxespad=0, frameon=False)
        fig.savefig(out_pdf); plt.close(fig)
    return printed


def panel_g(ab, st, out_pdf):
    """Fig. 4g: quantified protein groups off and on antibiotics, 42 patients.

    Producer: four_author_integration_2026-10-02/actual_FAR_261002/build_FAR23.py (loop
    entry ``g``). The producer read the A33 copies of the paired table and statistics;
    their ``pg_parsimony`` values and the pg_parsimony statistics row are identical to
    the clinical_refit files used here (PROVENANCE.md).

    Parameters
    ----------
    ab : pandas.DataFrame
        Controlled ``paired_antibiotic_profiles.tsv`` (84 rows = 42 patients x 2 states).
    st : pandas.DataFrame
        ``clinical_paired_statistics.tsv`` indexed by metric.
    out_pdf : pathlib.Path
        PDF to write.

    Returns
    -------
    list of dict
        matplotlib box statistics for off and on (checked against numpy percentiles).
    """
    with mpl.rc_context():
        mpl.rcdefaults()
        plt.rcParams.update({'font.family': _buch_family(), 'font.size': 6.4, 'axes.labelsize': 6.4, 'xtick.labelsize': 5.75, 'ytick.labelsize': 6.4,
                             'axes.linewidth': .5, 'xtick.major.width': .5, 'ytick.major.width': .5, 'ytick.major.size': 2,
                             'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'savefig.bbox': None,
                             'text.color': '#333333', 'axes.labelcolor': '#333333', 'xtick.color': '#333333', 'ytick.color': '#333333'})
        letter, metric, scale, ylim = 'g', 'pg_parsimony', 1, 19500
        yticks, ylabels, ylabel, star = [0, 4000, 8000, 12000, 16000], ['0', '4k', '8k', '12k', '16k'], 'Quantified protein groups', '**'
        w = ab.pivot(index='patient_id', columns='antibiotic_state', values=metric)[['Off', 'On']].dropna() * scale
        assert len(w) == 42 and np.isfinite(w.to_numpy()).all()
        s = st.loc[metric]
        assert np.isclose((w.On - w.Off).mean(), s.mean_on_minus_off, rtol=0, atol=1e-10)
        assert np.isclose(w.Off.mean(), s.mean_off, rtol=0, atol=1e-10)
        fig = plt.figure(figsize=(136 / 72, 120 / 72)); ax = fig.add_axes([31.85 / 136, 17.03 / 120, 98.54 / 136, 88.19 / 120])
        for a, b in w.itertuples(index=False, name=None):
            ax.plot([0, 1], [a, b], color='#cfcfcf', lw=.4, zorder=1)
        for x, c in [(0, '#888888'), (1, '#D55E00')]:
            ax.scatter([x] * 42, w.iloc[:, x], s=4.2, color=c, lw=0, zorder=3)
        box = ax.boxplot([w.Off, w.On], positions=[0, 1], widths=.46, whis=BOX_WHIS, showfliers=False, patch_artist=False,
                         boxprops={'color': '#333333', 'linewidth': .65}, whiskerprops={'color': '#333333', 'linewidth': .65},
                         capprops={'color': '#333333', 'linewidth': .65}, medianprops={'color': 'black', 'linewidth': 1.1}, zorder=4)
        stats = cbook.boxplot_stats([w.Off, w.On], whis=BOX_WHIS)
        for k, sb in enumerate(stats):
            q = np.percentile(w.iloc[:, k], [25, 50, 75]); assert np.allclose(q, [sb['q1'], sb['med'], sb['q3']], rtol=0, atol=1e-12)
            assert np.allclose(box['medians'][k].get_ydata(), [q[1], q[1]])
        ax.set_xlim(-.5, 1.5); ax.set_ylim(0, ylim); ax.set_xticks([0, 1], ['Antibiotics\noff', 'Antibiotics\non']); ax.tick_params(axis='x', length=0, pad=1.5)
        ax.set_yticks(yticks, ylabels); ax.tick_params(axis='y', pad=2)
        ax.set_ylabel(ylabel, labelpad=2); ax.set_title(letter, loc='left', fontsize=8.13, pad=2)
        ax.plot([0, 0, 1, 1], np.array([.93, .96, .96, .93]) * ylim, color='#333333', lw=.6)
        ax.text(.5, .97 * ylim, star, ha='center', va='bottom', fontsize=6.4)
        annotation = f'{s.mean_on_minus_off:,.0f} groups\n95% CI {s.CI95_low:,.0f} to {s.CI95_high:,.0f}'.replace('-', '−')
        ax.text(.5, .86 * ylim, annotation, ha='center', va='top', fontsize=5.3, linespacing=1.05)
        fig.savefig(out_pdf); plt.close(fig)
    return stats


# ----------------------------------------------------------------------------------
# Panel h
# ----------------------------------------------------------------------------------
def _panel_h_figure(Z20, feat, tracks, abx):
    """Build the panel-h figure exactly as final/fig_v27/build_fig4_v27.py (panel h).

    Returns the open figure and the colour-bar x position and width in points (used by
    the panel-i colour bar, which fig_v27 aligned with this one).
    """
    from style import fl_style as fl
    fl.init()
    W = PAGE_W
    keep = np.flatnonzero(feat.n_positive_85 > 0); assert len(keep) == 19
    labels = [s.replace('Dimethylsuccinic acid', 'Dimethylsuccinic / methylglutaric / adipic acid') for s in feat.label.iloc[keep]]
    Z = Z20[keep].astype(float)
    blocks = feat.block.iloc[keep].tolist()
    assert labels[0] == 'HSPB6' and labels[4] == 'SODM' and labels[5] == 'Nitrogen metabolism' \
        and labels[10] == 'Alanine, aspartate and glutamate metabolism' and labels[11] == 'Hydroxyindoleacetic acid' and labels[18] == 'Dihydrothymine'
    assert len(set(blocks[0:5])) == len(set(blocks[5:11])) == len(set(blocks[11:19])) == 1
    GROUPS = [('Host\nproteins', 0, 4), ('Microbial\npathways', 5, 10), ('Metabolites', 11, 18)]
    HH = 218.0; FS = 6.4
    cm = LinearSegmentedColormap.from_list('fl', [fl.BLUE, 'white', fl.ORANGE]); cm.set_bad('#dedede')
    fig = plt.figure(figsize=(W / 72, HH / 72))
    pt = lambda x, y, w, h: [x / W, 1 - (y + h) / HH, w / W, h / HH]  # noqa: E731
    hx, hw = 0.0, 300.0
    row = 8.5; top = 5.0; trk_h = 7.0; trk_gap = 2.6
    trks = [(tracks.eds_score.to_numpy(), H_EDS_TRACK_LABEL, plt.get_cmap('Greys'), 1),
            (100 * tracks.human_fraction.to_numpy(), 'Human-assigned LFQ share (%)', LinearSegmentedColormap.from_list('hum', ['white', fl.PURPLE]), 100),
            (pd.Series(abx).map({'Off': 0, 'On': 1}).to_numpy(), 'Antibiotics (dark: on)', ListedColormap(['#ededed', '#333333']), 1)]
    axes = []
    for i, (v, lab, c, vmax) in enumerate(trks):
        ax = fig.add_axes(pt(hx, top + i * (trk_h + trk_gap), hw, trk_h))
        ax.imshow(v[None, :], cmap=c, vmin=0, vmax=vmax, aspect='auto', interpolation='nearest')
        ax.set(xticks=[], yticks=[0]); ax.set_yticklabels([lab], fontsize=FS); ax.tick_params(length=0, pad=2)
        for sp in ax.spines.values():
            sp.set_visible(False)
        axes.append(ax)
    y0 = top + 3 * (trk_h + trk_gap) + 1.0
    ax = fig.add_axes(pt(hx, y0, hw, 19 * row)); axes.append(ax); hm = ax
    im = ax.imshow(np.clip(Z, -HEATMAP_CLIP, HEATMAP_CLIP), cmap=cm, vmin=-2, vmax=2, aspect='auto', interpolation='nearest')
    ax.set_yticks(range(19)); ax.set_yticklabels(labels, fontsize=FS); ax.tick_params(axis='y', length=0, pad=2)
    ax.set_xticks([0, 21, 42, 63, 84]); ax.set_xticklabels(['1', '22', '43', '64', '85'], fontsize=FS)
    ax.tick_params(axis='x', length=2, pad=1.5)
    ax.set_xlabel('85 patients ordered by metagenomic Enterotype Dysbiosis Score', fontsize=FS, labelpad=2)
    for s in ['left', 'top', 'right']:
        ax.spines[s].set_visible(False)
    for sep in [4.5, 10.5]:
        ax.axhline(sep, color=fl.INK, lw=fl.MIN_LINE_PT)
    cax = fig.add_axes(pt(hx + hw + 9, y0 + 19 * row * .25, 5.5, 19 * row * .5)); axes.append(cax)
    cb = fig.colorbar(im, cax=cax, ticks=[-2, 0, 2]); cb.outline.set_linewidth(fl.MIN_LINE_PT)
    cb.set_label('Within-feature z score', fontsize=FS, labelpad=2); cb.ax.tick_params(length=2, pad=1, labelsize=FS)
    fig.canvas.draw(); r = fig.canvas.get_renderer()
    lab_x0 = min(t.get_window_extent(r).x0 for t in hm.get_yticklabels())
    wlab = (hm.bbox.x0 - lab_x0) / fig.dpi * 72
    BR = wlab + 4.0
    TICK = 2.2; TXT_GAP = 2.0
    probe = fig.text(0, 0, 'Microbial\npathways', fontsize=FS, rotation=90, linespacing=1.05); fig.canvas.draw()
    txt_w = probe.get_window_extent(r).width / fig.dpi * 72; probe.remove()
    bbs = [a.get_tightbbox(r) for a in axes]
    x0 = min(min(b.x0 for b in bbs) / fig.dpi * 72, hm.bbox.x0 / fig.dpi * 72 - BR - TXT_GAP - txt_w)
    x1 = max(b.x1 for b in bbs) / fig.dpi * 72
    dx = (W - (x1 - x0)) / 2 - x0
    for a in axes:
        p = a.get_position(); a.set_position([p.x0 + dx / W, p.y0, p.width, p.height])
    fig.canvas.draw()
    hm_left = hm.get_position().x0 * W

    def ypt(row_edge):
        """Page y (points from the top) of a heatmap row edge."""
        return y0 + (row_edge + .5) * row
    for name, a, b in GROUPS:
        ya, yb = ypt(a - .5) + 1.0, ypt(b + .5) - 1.0
        xb = hm_left - BR
        for xs, ys in [([xb, xb], [ya, yb]), ([xb, xb + TICK], [ya, ya]), ([xb, xb + TICK], [yb, yb])]:
            fig.add_artist(Line2D([x / W for x in xs], [1 - y / HH for y in ys], transform=fig.transFigure,
                                  color=fl.INK, lw=fl.MIN_LINE_PT, solid_capstyle='butt'))
        fig.text((xb - TXT_GAP) / W, 1 - ((ya + yb) / 2) / HH, name, rotation=90, ha='right', va='center',
                 fontsize=FS, fontweight='normal', color='black', multialignment='center', linespacing=1.05)
    fig.text(3 / W, 1 - 3 / HH, 'h', fontsize=9, va='top', ha='left', fontweight='normal', color='black')
    fig.canvas.draw()
    hcb = cax.get_position()
    return fig, hcb.x0 * W, hcb.width * W


def panel_h(Z20, feat, tracks, abx, out_pdf):
    """Fig. 4h: 19 host-protein, microbial-pathway and metabolite features across 85 patients.

    Producer: final/fig_v27/build_fig4_v27.py, panel h. Needs the controlled antibiotic
    track; without it nothing is drawn. The top track label is H_EDS_TRACK_LABEL, the v30
    wording; v30 changed it on the composite (final/fig4_eds_label/relabel.py: white box
    over 'EDS (0–1)', new text right-aligned at the same baseline), here it is drawn by
    the producer code directly.

    Parameters
    ----------
    Z20 : numpy.ndarray
        Within-feature z scores, 20 features x 85 patients (display order).
    feat, tracks : pandas.DataFrame
        ``heatmap_features.tsv``, ``heatmap_tracks.tsv``.
    abx : numpy.ndarray or None
        Antibiotic state per column (controlled).
    out_pdf : pathlib.Path
        PDF to write.

    Returns
    -------
    pathlib.Path or None
        ``out_pdf``, or None when the antibiotic track is absent and nothing is drawn.
    """
    if abx is None:
        print("panel h SKIPPED: antibiotic track is a controlled input")
        return None
    with mpl.rc_context():
        mpl.rcdefaults()
        fig, _, _ = _panel_h_figure(Z20, feat, tracks, abx)
        fig.savefig(out_pdf); plt.close(fig)
    return out_pdf


# ----------------------------------------------------------------------------------
# Panel i
# ----------------------------------------------------------------------------------
def panel_i(native, functional, out_pdf):
    """Fig. 4i: patient-profile distance correlations (77 patients), native, eggNOG OG, KEGG KO.

    Producer: iteration4/figures/shadow_d12/B/build_figures.py, ``pyramid('PREDICT',
    'panel_i', 518.74, 160, ['native','OG','KO'], 'i')`` (big=False path). The title
    "PREDICT · 77 patients" and the two footnotes it draws were clipped out of the
    composite by fig_v27 (ASSEMBLY.md).

    Parameters
    ----------
    native : pandas.DataFrame
        ``crossomics_native_correlations.tsv``.
    functional : pandas.DataFrame
        ``crossomics_functional_correlations.tsv``.
    out_pdf : pathlib.Path
        PDF to write.

    Returns
    -------
    list of dict
        One record per drawn off-diagonal cell (rho, q, printed label).
    """
    blue = '#0072B2'; orange = '#E69F00'
    cm = LinearSegmentedColormap.from_list('fl', [blue, 'white', orange]); cm.set_bad('#dedede')
    records = []
    with mpl.rc_context():
        mpl.rcdefaults(); _d12_rc(_buch_family())
        co, name, width, height, levels, letter, big = 'PREDICT', 'panel_i', PAGE_W, 160, ['native', 'OG', 'KO'], 'i', False
        fs = (lambda small, large: large if big else small)  # noqa: E731
        layers = ['DNA', 'RNA', 'Protein', 'Human protein', 'Metabolites'] if co == 'PREDICT' else ['DNA', 'Protein', 'Human protein', 'Metabolites']; n = len(layers)
        km = {'DNA': 'metaG', 'RNA': 'metaT', 'Protein': 'microbial_metaP', 'Human protein': 'host_metaP', 'Metabolites': 'metabolites'}
        kn = {'DNA': 'DNA', 'RNA': 'RNA', 'Protein': 'Non-host protein', 'Human protein': 'Host protein', 'Metabolites': 'Metabolites'}
        gx = 76; gy = 32; cell = min(fs(18, 27), (width - gx - 12) / (len(levels) * n + len(levels) - 1)); gap = cell
        fig = plt.figure(figsize=(width / 72, height / 72)); ax = fig.add_axes([0, 0, 1, 1]); ax.set(xlim=(0, width), ylim=(height, 0)); ax.axis('off')
        ax.text(3, 10, letter, weight='bold', fontsize=8)
        ax.text(17, 10, 'PREDICT · 77 patients' if co == 'PREDICT' else 'CapScan · provisional frozen 12-participant roster', fontsize=fs(7, 8.5))
        for col, lev in enumerate(levels):
            x = gx + col * (n * cell + gap); ax.text(x + n * cell / 2, 25, {'native': 'Native', 'OG': 'eggNOG OG', 'KO': 'KEGG KO', 'Pathway': 'KEGG pathways'}[lev], ha='center', fontsize=fs(7, 8.5))
            t = native[(native.cohort == co) & (native.analysis == 'sensitivity_all_features')] if lev == 'native' else \
                functional[(functional.cohort == co) & (functional.level == lev) & (functional.prevalence == 0) & (functional.scope == ('displayed77' if co == 'PREDICT' else 'matched12'))]
            maps = kn if lev == 'native' else km; lookup = {frozenset((r.layer_a, r.layer_b)): r for r in t.itertuples()}
            for i, la in enumerate(layers):
                if col == 0:
                    ax.text(gx - 5, gy + (i + .5) * cell, f'{i + 1} ' + la, ha='right', va='center', fontsize=fs(6, 8))
                ax.text(x + (i + .5) * cell, gy + n * cell + 5, str(i + 1), ha='center', va='top', fontsize=fs(6, 8))
                for j, lb in enumerate(layers[:i + 1]):
                    r = lookup.get(frozenset((maps[la], maps[lb])))
                    if i == j:
                        rho = 1.; label = '1.00'
                    elif r is None:
                        rho = 0.; label = '—'; assert lev != 'native' and {la, lb} == {'Human protein', 'Metabolites'}
                    else:
                        rho = r.rho; label = f'{rho:.2f}' + ('*' if r.q < Q_STAR else '')
                        records.append(dict(cohort=co, panel=name, level=lev, layer_a=la, layer_b=lb, rho=rho, q=r.q, label=label))
                    ax.add_patch(Rectangle((x + j * cell, gy + i * cell), cell, cell, facecolor=cm((rho + 1) / 2) if label != '—' else '#f5f5f5', edgecolor='white', lw=.5))
                    ax.text(x + (j + .5) * cell, gy + (i + .5) * cell, label, ha='center', va='center', fontsize=fs(5.9, 8))
        y = gy + n * cell + fs(19, 21)
        ax.text(3, y, 'Protein excludes human-assigned signal; unresolved origin is retained.', fontsize=fs(5.6, 8), va='top')
        ax.text(3, y + fs(10, 12), 'All observed features. * q < 0.05. Complete families and prevalence sensitivities in Source Data.', fontsize=fs(5.4, 8), va='top')
        fig.savefig(out_pdf); plt.close(fig)
    return records


def panel_i_colourbar(feat, tracks, out_pdf):
    """Fig. 4i colour scale (-1 to 1, 'Spearman ρ'), placed on the full composite page.

    Producer: final/fig_v27/build_fig4_v27.py, "panel i colour scale". Its x position and
    width are those of the panel-h colour bar, which depend only on the panel-h layout
    (labels and fonts), not on any data value; the layout is rebuilt here with the
    public tracks and an all-off antibiotic track that is never drawn or saved.

    Parameters
    ----------
    feat, tracks : pandas.DataFrame
        ``heatmap_features.tsv``, ``heatmap_tracks.tsv`` (layout only).
    out_pdf : pathlib.Path
        PDF to write (full composite page, 518.74 x 688 pt).

    Returns
    -------
    dict
        Colour-bar position and size in points (x, y, w, h) and the KO block's right edge.
    """
    with mpl.rc_context():
        mpl.rcdefaults()
        fig_h, H_CBAR_X0, H_CBAR_W = _panel_h_figure(load_heatmap()[0], feat, tracks, np.array(['Off'] * 85))
        plt.close(fig_h)
        from style import fl_style as fl
        W, NH, FS = PAGE_W, PAGE_H_V27, 6.4
        I_TOP, grow = 549.0, 218.0 - (546.0 - 331.0)
        i_top = I_TOP + grow
        gx, gy, cell = 76, 32, 18.0
        mat_top = i_top + gy
        ko_right = gx + 2 * (5 * cell + cell) + 5 * cell + I_DXI
        CB_X = max(H_CBAR_X0, ko_right + 12); CB_W = H_CBAR_W
        CB_Y0 = mat_top + 5 * cell * .25; CB_H = 5 * cell * .5
        fig = plt.figure(figsize=(W / 72, NH / 72))
        cax = fig.add_axes([CB_X / W, 1 - (CB_Y0 + CB_H) / NH, CB_W / W, CB_H / NH])
        cmi = LinearSegmentedColormap.from_list('fl', ['#0072B2', 'white', '#E69F00'])
        sm = plt.cm.ScalarMappable(norm=Normalize(-1, 1), cmap=cmi)
        cbi = fig.colorbar(sm, cax=cax, ticks=[-1, 0, 1]); cbi.outline.set_linewidth(fl.MIN_LINE_PT)
        cbi.ax.set_yticklabels(['−1', '0', '1'])
        cbi.ax.tick_params(length=2, pad=1, labelsize=FS)
        cbi.set_label('Spearman ρ', fontsize=FS, labelpad=2, fontfamily=['soehne', 'DejaVu Sans'])
        fig.savefig(out_pdf, transparent=True); plt.close(fig)
    return dict(x=CB_X, y=CB_Y0, w=CB_W, h=CB_H, ko_block_right=ko_right)


# ----------------------------------------------------------------------------------
# Statistics
# ----------------------------------------------------------------------------------
def _gain(fl, mg):
    """Per-acquisition gain 100 x (FL / MG - 1), NaN where MG is 0. Unit: acquisition."""
    fl, mg = np.asarray(fl, float), np.asarray(mg, float)
    with np.errstate(divide='ignore', invalid='ignore'):
        return np.where(mg > 0, 100 * (fl / mg - 1), np.nan)


def stat_bc_crossfit_bootstrap(x, strata=None, seed=BOOT_SEED_CROSSFIT, draws=BOOT_DRAWS_CROSSFIT):
    """Median per-acquisition gain and participant cluster bootstrap CI (crossfit analyze_final.py).

    Unit: acquisition (median over acquisitions); resampling unit: participant (whole
    participants drawn with replacement; acquisitions without a participant ID are their
    own cluster). One RNG stream, strata in fixed order, the same resampled index sets for
    every column within a stratum, exactly as the original ``summarize``.
    n: PREDICT 296 acquisitions / 97 clusters; CapScan 282 / 15; stool 56 / 15;
    intestinal fluid 226 / 14; regions 51-65 / 13-14.

    Parameters
    ----------
    x : pandas.DataFrame
        :func:`load_crossfit` output.
    strata : list of (name, mask) or None
        Default: the strata of CROSSFIT_STRATA, in that order.
    seed : int
        Seed of the single random stream (default BOOT_SEED_CROSSFIT).
    draws : int
        Bootstrap resamples (default BOOT_DRAWS_CROSSFIT).

    Returns
    -------
    pandas.DataFrame
        One row per stratum x length x comparison (gain_cf_pooled, gain_frozen).
    """
    cap = x.cohort.eq('CAPSCAN')
    if strata is None:
        masks = {'PREDICT': x.cohort.eq('PREDICT'), 'CAPSCAN_all': cap, 'CAPSCAN_stool': cap & x.sampling_site.eq('stool'),
                 'CAPSCAN_intestinal': cap & x.sampling_site.ne('stool')}
        strata = [(name, masks[name]) for name in CROSSFIT_STRATA]
    rng = np.random.default_rng(seed)
    out = []
    for name, m in strata:
        g = x[m].reset_index(drop=True)
        subj = g.participant_sort_rank.to_numpy(); us = np.unique(subj)
        members = [np.flatnonzero(subj == s) for s in us]
        ix_draws = [np.concatenate([members[k] for k in rng.integers(0, len(us), len(us))]) for _ in range(draws)]
        for L in LENGTHS:
            cols = {'gain_cf_pooled': _gain(g[f'FL_pooled_{L}'], g[f'MG_pooled_{L}']),
                    'gain_frozen': _gain(g[f'frozen_fl_{L}'], g[f'frozen_mg_{L}'])}
            for col, v in cols.items():
                bs = np.array([np.nanmedian(v[ix]) for ix in ix_draws])
                out.append(dict(stratum=name, min_length=L, comparison=col, n_acquisitions=len(g), n_participants=len(us),
                                median=float(np.nanmedian(v)), ci95_low=float(np.percentile(bs, CI_PERCENTILES[0])),
                                ci95_high=float(np.percentile(bs, CI_PERCENTILES[1])), n_pos=int((v > 0).sum()), n_neg=int((v < 0).sum())))
    return pd.DataFrame(out)


def stat_bc_region_bootstrap(x):
    """Participant bootstrap of median gains per CapScan sampling site, on its own stream.

    Reproduces SUMMARY_CAPSCAN_BY_REGION.tsv of the cross-fit analysis: a fresh stream with
    seed BOOT_SEED_CROSSFIT over CAPSCAN_REGIONS in that order (v28d line 777). Unit:
    acquisition; resampling unit: participant. n: duodenum 56 acquisitions / 13
    participants, jejunum 65 / 14, ileum 54 / 14, colon 51 / 13, stool 56 / 15.

    Parameters
    ----------
    x : pandas.DataFrame
        :func:`load_crossfit` output.

    Returns
    -------
    pandas.DataFrame
        As :func:`stat_bc_crossfit_bootstrap`, one stratum per sampling site.
    """
    cap = x.cohort.eq('CAPSCAN')
    return stat_bc_crossfit_bootstrap(x, [(f'CAPSCAN_{g}', cap & x.sampling_site.eq(g)) for g in CAPSCAN_REGIONS])


def stat_bc_recovery(x):
    """Median cross-fitted MG-identity recovery (%) at 9-50 residues per stratum.

    Unit: acquisition. n: 296 PREDICT, 56 CapScan stool, 226 intestinal-fluid
    acquisitions. No test.

    Parameters
    ----------
    x : pandas.DataFrame
        :func:`load_crossfit` output.

    Returns
    -------
    dict
        stratum -> median recovery (%).
    """
    return x.groupby('stratum')[f'recovery_cf_{L_BC}_pct'].median().to_dict()


def stat_bc_bootstrap_seed_sensitivity(x, seeds=BOOT_SENSITIVITY_SEEDS):
    """Monte Carlo sensitivity of the CapScan stool 9-50 interval to the random stream.

    Repository addition, not in the manuscript. The printed stool interval (v28d line 123)
    comes from one stream shared with three other strata; this reruns the same
    participant bootstrap (BOOT_DRAWS_CROSSFIT draws) for the stool stratum alone, with
    BOOT_SEED_CROSSFIT and each seed in ``seeds``. Unit: acquisition (56); resampling unit:
    participant (15).

    Parameters
    ----------
    x : pandas.DataFrame
        :func:`load_crossfit` output.
    seeds : iterable of int
        Further seeds (default BOOT_SENSITIVITY_SEEDS).

    Returns
    -------
    pandas.DataFrame
        One row per seed: seed, ci95_low, ci95_high (gain_cf_pooled, 9-50 residues).
    """
    m = x.cohort.eq('CAPSCAN') & x.sampling_site.eq('stool')
    out = []
    for s in (BOOT_SEED_CROSSFIT, *seeds):
        b = stat_bc_crossfit_bootstrap(x, [('CAPSCAN_stool', m)], seed=s)
        r = b[(b.min_length == L_BC) & (b.comparison == 'gain_cf_pooled')].iloc[0]
        out.append(dict(seed=s, ci95_low=r.ci95_low, ci95_high=r.ci95_high))
    return pd.DataFrame(out)


def stat_bc_same_vs_crossfit(x, L=SAME_VS_CROSSFIT_LENGTH):
    """Acquisitions whose same-spectrum gain exceeds the cross-fitted gain.

    Unit: acquisition. n: 296 PREDICT, 282 CapScan. No test (count). The v30 text
    quotes 9-50 residues (296 of 296, 282 of 282); at 5-50 the counts are 295 and 281
    (v28d line 133). Rule as in LAB_ROOT/.../final/len950/Fig4/build/same_vs_crossfit.py.

    Parameters
    ----------
    x : pandas.DataFrame
        :func:`load_crossfit` output.
    L : str
        Minimum length suffix ('5' or '9'; default SAME_VS_CROSSFIT_LENGTH).

    Returns
    -------
    dict
        cohort -> (n exceeding, n acquisitions).
    """
    out = {}
    for coh, g in x.groupby('cohort'):
        same = _gain(g[f'frozen_fl_{L}'], g[f'frozen_mg_{L}']); cf = _gain(g[f'FL_pooled_{L}'], g[f'MG_pooled_{L}'])
        out[coh] = (int((same > cf).sum()), len(g))
    return out


def _mixed_fit(f, d):
    """Best ML fit of four optimisers (C3_region_test/region_test.py ``fit``, unchanged)."""
    import warnings
    import statsmodels.formula.api as smf
    warnings.filterwarnings('ignore')
    best, last = None, None
    for meth in ['lbfgs', 'bfgs', 'powell', 'nm']:
        try:
            m = smf.mixedlm(f, d, groups=d.subject_id).fit(reml=False, method=meth, maxiter=5000)
            cp = m.cov_params().loc[m.fe_params.index, m.fe_params.index].values
            if m.converged and np.isfinite(m.llf) and np.all(np.isfinite(cp)) and np.all(np.linalg.eigvalsh(cp) > 0):
                if best is None or m.llf > best.llf + 1e-9:
                    m.method_used = meth; best = m
        except Exception as e:  # noqa: BLE001 - as in the original
            last = e
    if best is None:
        raise RuntimeError(f'no method converged for {f}: {last}')
    return best


def stat_mixed_model(x, view='crossfit', L=MIXED_MODEL_LENGTH):
    """Linear mixed model of ln(FL/MG) on sampling site with participant random intercept.

    Port of C3_region_test/region_test.py for one view and length (log scale). Unit:
    acquisition (282 CapScan; 15 participants, participant 15 stool only); repeated:
    acquisitions within participant. ML fit; Wald z for the pooled fluid-vs-stool
    contrast (mean of the four region coefficients), 3-df Wald for the region omnibus,
    4-df LRT against the intercept-only model. No multiplicity adjustment is applied to
    any of these (the upstream script applies Holm only to the four region-vs-stool and
    six pairwise region contrasts, which the manuscript does not print). Methods v28d
    line 229 describes an LRT omnibus with Holm-adjusted Wald contrasts (OPEN_ITEMS
    MISMATCH).

    Parameters
    ----------
    x : pandas.DataFrame
        :func:`load_crossfit` output.
    view : str
        'crossfit' (pooled cross-fitted counts) or any other value for same-spectrum counts.
    L : str
        Minimum length suffix ('5' or '9'; default MIXED_MODEL_LENGTH).

    Returns
    -------
    dict
        estimates (log scale and fold), Wald CI, z, P, chi-square statistics.
    """
    import scipy.stats as st
    a = x[x.cohort.eq('CAPSCAN')].copy(); assert len(a) == 295
    a['g'] = a.sampling_site.map(REGION_LABEL); a['subject_id'] = a.participant_id.astype(str)
    fl, mg = ('FL_pooled_{L}', 'MG_pooled_{L}') if view == 'crossfit' else ('frozen_fl_{L}', 'frozen_mg_{L}')
    d = a[['subject_id', 'g']].copy(); d['y'] = np.log(a[fl.format(L=L)].astype(float) / a[mg.format(L=L)].astype(float))
    d = d[np.isfinite(d.y)]
    sd = d.y.std(); d['ys'] = d.y / sd
    m1 = _mixed_fit("ys ~ C(g, Treatment('stool'))", d); m0 = _mixed_fit("ys ~ 1", d)
    lr = 2 * (m1.llf - m0.llf)
    REG = ['R1', 'R2', 'R3', 'R4']
    names = [f"C(g, Treatment('stool'))[T.{r}]" for r in REG]
    b = m1.fe_params[names].values * sd; V = m1.cov_params().loc[names, names].values * sd ** 2
    C3 = np.array([[-1, 1, 0, 0], [-1, 0, 1, 0], [-1, 0, 0, 1]], float)
    cb = C3 @ b; w3 = float(cb @ np.linalg.solve(C3 @ V @ C3.T, cb))
    c = np.array([.25] * 4); est = float(c @ b); se = float(np.sqrt(c @ V @ c)); z = est / se
    return dict(view=view, length=f'{L}-50', n_acq=len(d), n_part=d.subject_id.nunique(), fit_method=m1.method_used,
                fluid_minus_stool_log=est, se=se, z=z, p=2 * st.norm.sf(abs(z)),
                fold=float(np.exp(est)), fold_ci_low=float(np.exp(est - 1.96 * se)), fold_ci_high=float(np.exp(est + 1.96 * se)),
                region_wald_chi2=w3, region_wald_df=3, region_p=float(st.chi2.sf(w3, 3)),
                lrt_chi2=float(lr), lrt_df=4, lrt_p=float(st.chi2.sf(lr, 4)),
                participant_var_frac=float(m1.cov_re.iloc[0, 0] / (m1.cov_re.iloc[0, 0] + m1.scale)))


def _bh(p):
    """Benjamini-Hochberg (clinical_refit.py ``bh``)."""
    p = np.asarray(p, dtype=float)
    order = np.argsort(p); ranked = p[order] * len(p) / np.arange(1, len(p) + 1)
    out = np.empty_like(p); out[order] = np.minimum(1, np.minimum.accumulate(ranked[::-1])[::-1])
    return out


def stat_eg_paired(ab):
    """Paired on-minus-off contrast for human-assigned share and protein groups (clinical_refit.py ``derive``).

    Unit: patient (n = 42, each measured off and on antibiotics; per-state values are
    means (human share) or medians (protein groups) over that patient's acquisitions).
    Estimate: mean paired difference. CI: percentile bootstrap, 50,000 draws resampling
    patients, seed 20260907. Test: two-sided Monte Carlo sign-flip permutation of the
    mean difference, 100 x 2,000 flips, seed 20260907 + endpoint index,
    P = (count + 1) / 200,001. Correction: Benjamini-Hochberg over the two endpoints.

    Parameters
    ----------
    ab : pandas.DataFrame
        Controlled ``paired_antibiotic_profiles.tsv`` (84 rows = 42 patients x 2 states).

    Returns
    -------
    pandas.DataFrame
        One row per endpoint.
    """
    stats = []
    for idx, (metric, scale) in enumerate([('Human', 100), ('pg_parsimony', 1)]):
        wide = ab.pivot(index='patient_id', columns='antibiotic_state', values=metric)[['Off', 'On']] * scale
        delta = (wide.On - wide.Off).to_numpy()
        rng = np.random.default_rng(CLIN_SEED)
        boots = delta[rng.integers(0, len(delta), (CLIN_BOOT_DRAWS, len(delta)))].mean(axis=1)
        lo, hi = np.quantile(boots, [CI_PERCENTILES[0] / 100, CI_PERCENTILES[1] / 100])
        rng = np.random.default_rng(CLIN_SEED + idx); count = 0
        for _ in range(CLIN_SIGNFLIP_DRAWS // 2000):
            null = ((rng.integers(0, 2, (2000, len(delta))) * 2 - 1) * delta).mean(axis=1)
            count += int((abs(null) >= abs(delta.mean()) - 1e-12).sum())
        d_sd = delta.std(ddof=1)
        stats.append(dict(metric=metric, n_patients=len(delta), mean_off=wide.Off.mean(), mean_on=wide.On.mean(),
                          mean_on_minus_off=delta.mean(), CI95_low=lo, CI95_high=hi, p_signflip_MC=(count + 1) / (CLIN_SIGNFLIP_DRAWS + 1),
                          n_increase=int((delta > 0).sum()), sd_difference=d_sd, cohens_dz=delta.mean() / d_sd))
    stats = pd.DataFrame(stats); stats['q_two_clinical_endpoints'] = _bh(stats.p_signflip_MC)
    return stats


def _rank_norm(dist):
    """Centre and unit-normalise the ranks of a condensed distance vector (dot product = Spearman rho)."""
    from scipy.stats import rankdata
    r = rankdata(dist); r = r - r.mean(); return r / np.linalg.norm(r)


def stat_i_native(dist):
    """Mantel-type Spearman correlations between native profile distances (native_pyramids.py ``correlations``).

    Unit: patient (77; distance pairs 2,926). Statistic: Spearman rho of the upper
    triangles. Test: two-sided permutation of patient labels of one matrix, 1,999
    permutations, seed 20260907, P = (1 + extreme) / 2000. Correction: Benjamini-
    Hochberg over the 10 layer pairs.

    Parameters
    ----------
    dist : pandas.DataFrame
        ``crossomics_distances_predict77.tsv`` (2,926 patient pairs; columns ``native__<layer>``).

    Returns
    -------
    pandas.DataFrame
        One row per layer pair: rho, permutation P, n, BH q.
    """
    from scipy.spatial.distance import squareform
    n = 77; ix = np.triu_indices(n, 1)
    layers = ['DNA', 'RNA', 'Non-host protein', 'Host protein', 'Metabolites']
    rng = np.random.default_rng(MANTEL_SEED)
    perms = np.array([rng.permutation(n) for _ in range(MANTEL_PERMUTATIONS)])
    norm = {k: _rank_norm(dist[f'native__{k}'].to_numpy()) for k in layers}
    rows = []
    for a, b in itertools.combinations(layers, 2):
        rho = float(norm[a] @ norm[b]); bmat = squareform(norm[b]); extreme = 0
        for chunk in np.array_split(perms, 20):
            null = bmat[chunk[:, ix[0]], chunk[:, ix[1]]] @ norm[a]
            extreme += int((abs(null) >= abs(rho) - 1e-12).sum())
        rows.append(dict(level='native', layer_a=a, layer_b=b, rho=rho, p=(1 + extreme) / 2000, n=n))
    t = pd.DataFrame(rows); t['q'] = _bh(t.p); return t


def stat_i_functional(dist):
    """Functional-level distance correlations, PREDICT, all observed features (check_prevalence.py ``evaluate``).

    Unit: patient (77). Same statistic and permutation scheme as :func:`stat_i_native`
    (a fresh RNG with seed 20260907 per level, P = (1 + extreme) / 2000). The host-protein
    vs metabolite pair is not tested at KO and pathway level (identical profiles).
    Correction: Benjamini-Hochberg over the 28 tests of the three levels.

    Parameters
    ----------
    dist : pandas.DataFrame
        ``crossomics_distances_predict77.tsv`` (columns ``<level>__<layer>``).

    Returns
    -------
    pandas.DataFrame
        One row per level and layer pair (28): rho, permutation P, n, BH q.
    """
    from scipy.spatial.distance import squareform
    from statsmodels.stats.multitest import multipletests
    n = 77; ix = np.triu_indices(n, 1)
    layers = ['metaG', 'metaT', 'microbial_metaP', 'host_metaP', 'metabolites']
    rows = []
    for lev in ['OG', 'KO', 'Pathway']:
        rng = np.random.default_rng(MANTEL_SEED)
        perms = np.array([rng.permutation(n) for _ in range(MANTEL_PERMUTATIONS)])
        ranks = {k: _rank_norm(dist[f'{lev}__{k}'].to_numpy()) for k in layers}
        for a, b in itertools.combinations(layers, 2):
            if lev != 'OG' and {a, b} == {'host_metaP', 'metabolites'}:
                continue
            rho = float(ranks[a] @ ranks[b]); bm = squareform(ranks[b]); ext = 0
            for start in range(0, MANTEL_PERMUTATIONS, 200):
                pp = perms[start:start + 200]; null = bm[pp[:, ix[0]], pp[:, ix[1]]] @ ranks[a]; ext += np.count_nonzero(np.abs(null) >= abs(rho))
            rows.append(dict(level=lev, layer_a=a, layer_b=b, rho=rho, p=(1 + ext) / (1 + MANTEL_PERMUTATIONS), n=n))
    t = pd.DataFrame(rows); t['q'] = multipletests(t.p, method='fdr_bh')[1]; return t


def stat_text_capscan(x, dna):
    """Text line 129 descriptive values (same-spectrum gains by site, DNA amounts, database sizes).

    Units: acquisition (282 CapScan; 56 stool, 51-65 per region; DNA recorded for 99);
    participant (14 with databases in both compartments) for the database comparison.
    No test.

    Parameters
    ----------
    x : pandas.DataFrame
        :func:`load_crossfit` output.
    dna : pandas.DataFrame
        ``capscan_dna_yield_by_acquisition.tsv`` (282 acquisitions).

    Returns
    -------
    dict
        Named values (medians, counts, participant lists).
    """
    cap = x[x.cohort.eq('CAPSCAN')]
    out = {}
    for L in LENGTHS:
        same = pd.Series(_gain(cap[f'frozen_fl_{L}'], cap[f'frozen_mg_{L}']), index=cap.index)
        for site, g in cap.groupby('sampling_site'):
            out[f'same_gain_median_{L}_{site}'] = float(same[g.index].median())
            out[f'same_n_{site}'] = len(g)
            out[f'same_n_more_{L}_{site}'] = int((g[f'frozen_fl_{L}'] > g[f'frozen_mg_{L}']).sum())
        fl_ = cap.sampling_site.ne('stool')
        out[f'same_n_more_{L}_fluid'] = int((cap.loc[fl_, f'frozen_fl_{L}'] > cap.loc[fl_, f'frozen_mg_{L}']).sum())
        out['n_fluid'] = int(fl_.sum())
    d = dna.copy()
    d['fluid'] = d.region.ne('stool')
    rec = d[d.dna_ng.notna()]
    out['n_dna_recorded'] = len(rec)
    out['n_dna_table'] = len(d)
    for site, g in rec.groupby('region'):
        out[f'dna_ng_median_{site}'] = float(g.dna_ng.median())
    out['db_median_fluid'] = float(d.loc[d.fluid, 'assembly_proteins'].median())
    out['db_median_stool'] = float(d.loc[~d.fluid, 'assembly_proteins'].median())
    per = d.groupby(['subject', 'fluid']).assembly_proteins.median().unstack()
    both = per.dropna()
    out['participants_both'] = len(both)
    out['participants_fluid_db_smaller'] = int((both[True] < both[False]).sum())
    out['participants_stool_only'] = sorted(per.index[per[True].isna()].astype(str).tolist())
    for col, key in (('MG_peptides', 'MG'), ('FL_peptides', 'FL')):
        out[f'yield_fluid_over_stool_{key}_pct'] = float(100 * d.loc[d.fluid, col].median() / d.loc[~d.fluid, col].median())
    return out


def stat_text_same_spectrum(bench, L=TEXT_LENGTH):
    """Same-spectrum benchmark values the v30 text prints at 9-50 residues.

    FastaLake databases selected from the spectra being searched (Supplementary Fig. 6d,e
    in v30). Unit: acquisition; n = 296 PREDICT, 282 CapScan (56 stool, 226 intestinal
    fluid). No test. Gain = median of FL/MG peptide-count ratios minus one; recovery =
    median fraction of MG peptide identities also identified by FL; yield = median
    intestinal-fluid peptides over median stool peptides (as
    ``figures/SFig09/panels.py`` and LAB_ROOT/.../final/len950/SFig09, 34.4% and 52.7%).

    Parameters
    ----------
    bench : pandas.DataFrame
        ``same_spectrum_benchmark.tsv``.
    L : int
        Minimum peptide length (default TEXT_LENGTH = 9).

    Returns
    -------
    dict
        ``<cohort>_gain_pct``, ``<cohort>_n_more``, ``<cohort>_n``, ``<cohort>_recovery_pct``,
        ``yield_fluid_over_stool_MG_pct``, ``yield_fluid_over_stool_FL_pct``.
    """
    out = {}
    for coh in ('PREDICT', 'CAPSCAN'):
        d = bench[(bench.cohort == coh) & (bench.min_length == L)]
        out[f'{coh}_n'] = len(d)
        out[f'{coh}_gain_pct'] = float(100 * (d.count_ratio.median() - 1))
        out[f'{coh}_n_more'] = int((d.n_fl > d.n_mg).sum())
        out[f'{coh}_recovery_pct'] = float(100 * d.mg_recovery_fraction.median())
    c = bench[(bench.cohort == 'CAPSCAN') & (bench.min_length == L)]
    assert len(c) == 295
    fl_ = c.sampling_site.ne('stool')
    for col, key in (('n_mg', 'MG'), ('n_fl', 'FL')):
        out[f'yield_fluid_over_stool_{key}_pct'] = float(100 * c.loc[fl_, col].median() / c.loc[~fl_, col].median())
    return out


# ----------------------------------------------------------------------------------
# Numbers, in NUMBERS.tsv order
# ----------------------------------------------------------------------------------
def sci_superscript(p, sig=1):
    """Format a P value the way the manuscript prints it, e.g. 2.9e-07 -> '3 × 10⁻⁷'.

    Parameters
    ----------
    p : float
        Value to format.
    sig : int
        Significant figures of the mantissa.

    Returns
    -------
    str
    """
    e = int(np.floor(np.log10(p))); m = p / 10 ** e
    m = round(m, sig - 1)
    if m >= 10:
        m, e = m / 10, e + 1
    sup = str.maketrans('-0123456789', '⁻⁰¹²³⁴⁵⁶⁷⁸⁹')
    ms = f'{m:.{sig - 1}f}' if sig > 1 else f'{int(m)}'
    return f'{ms} × 10{str(e).translate(sup)}'


def collect_numbers(x, boot, region, mm, paired, native, functional, audit, text, roster, bench, n_heatmap):
    """Return every recomputed number in NUMBERS.tsv row order as (key, value).

    Keys name the v28d line (``L<n>``) and/or the v30 paragraph (``P<n>``, index in the
    v30 visible text) where the number is printed; rows with only ``L`` are no longer
    printed in v30, rows with only ``P`` are new in v30. Same-spectrum text values are at
    9-50 residues (TEXT_LENGTH), as printed in v30.

    ``paired`` may be None (controlled input absent); its rows then carry the values of
    the public ``clinical_paired_statistics.tsv`` and are labelled as such.

    Parameters
    ----------
    x, boot, region, mm, paired, native, functional, audit, text, roster, bench, n_heatmap
        Outputs of the loaders and ``stat_*`` functions as called in the notebook.

    Returns
    -------
    list of tuple
        ``(key, value)``, full precision, in NUMBERS.tsv row order.
    """
    def bt(stratum, L, comp, col):
        """Look up one column of one row of the cross-fit bootstrap table."""
        r = boot[(boot.stratum == stratum) & (boot.min_length == L) & (boot.comparison == comp)].iloc[0]
        return r[col]
    reg9 = region[(region.min_length == L_BC) & (region.comparison == 'gain_cf_pooled') & (region.stratum != 'CAPSCAN_stool')]
    bp5 = bench[(bench.cohort == 'PREDICT') & (bench.min_length == 5)]
    bc5 = bench[(bench.cohort == 'CAPSCAN') & (bench.min_length == 5)]
    rec = stat_bc_recovery(x)
    svc = stat_bc_same_vs_crossfit(x)
    fluid_sites = ['duodenum', 'jejunum', 'ileum', 'colon']
    pr = paired.set_index('metric') if paired is not None else load_table('clinical_paired_statistics.tsv').set_index('metric')
    nat = native.set_index(['layer_a', 'layer_b'])
    fun = functional.set_index(['level', 'layer_a', 'layer_b'])
    ss = stat_text_same_spectrum(bench)
    rows = [
        ('L121/P63 PREDICT acquisitions', len(roster)),
        ('L123/P64 matched PREDICT acquisitions', int(roster.matched_metagenome.sum())),
        ('L123/P64 same-spectrum PREDICT median gain 9-50 (%)', ss['PREDICT_gain_pct']),
        ('L123/P64 PREDICT acquisitions with more FL peptides 9-50', ss['PREDICT_n_more']),
        ('L123/P64 PREDICT median MG-identity recovery 9-50 (%)', ss['PREDICT_recovery_pct']),
        ('L123/P65 cross-fitting folds (design constant, not computed)', ''),
        ('L123 cross-fit PREDICT median gain 9-50, signed (%)', bt('PREDICT', L_BC, 'gain_cf_pooled', 'median')),
        ('L123/P65 cross-fit PREDICT CI low', bt('PREDICT', L_BC, 'gain_cf_pooled', 'ci95_low')),
        ('L123/P65 cross-fit PREDICT CI high', bt('PREDICT', L_BC, 'gain_cf_pooled', 'ci95_high')),
        ('L123 cross-fit CapScan stool median gain 9-50, signed (%)', bt('CAPSCAN_stool', L_BC, 'gain_cf_pooled', 'median')),
        ('L123/P65 cross-fit stool CI low', bt('CAPSCAN_stool', L_BC, 'gain_cf_pooled', 'ci95_low')),
        ('L123/P65 cross-fit stool CI high', bt('CAPSCAN_stool', L_BC, 'gain_cf_pooled', 'ci95_high')),
        ('L123/P65 cross-fit intestinal fluid median gain 9-50 (%)', bt('CAPSCAN_intestinal', L_BC, 'gain_cf_pooled', 'median')),
        ('L123/P65 cross-fit fluid CI low', bt('CAPSCAN_intestinal', L_BC, 'gain_cf_pooled', 'ci95_low')),
        ('L123/P65 cross-fit fluid CI high', bt('CAPSCAN_intestinal', L_BC, 'gain_cf_pooled', 'ci95_high')),
        ('L123/P65 cross-fit lowest region median (%)', reg9['median'].min()),
        ('L123/P65 cross-fit highest region median (%)', reg9['median'].max()),
        ('P65 cross-fit recovery PREDICT (%)', rec['PREDICT']),
        ('P65 cross-fit recovery CapScan stool (%)', rec['CAPSCAN_stool']),
        ('P65 cross-fit recovery intestinal fluid (%)', rec['CAPSCAN_intestinal']),
        ('L123/P65 matched CapScan acquisitions', len(bc5)),
        ('L123/P65 same-spectrum CapScan median gain 9-50 (%)', ss['CAPSCAN_gain_pct']),
        ('P65 CapScan acquisitions with more FL peptides 9-50', ss['CAPSCAN_n_more']),
        ('L123/P65 CapScan median MG-identity recovery 9-50 (%)', ss['CAPSCAN_recovery_pct']),
        ('L127/P68 legend PREDICT acquisitions', int((x.stratum == 'PREDICT').sum())),
        ('L127 legend CapScan acquisitions', len(bc5)),
        ('P68 legend CapScan stool acquisitions', int((x.stratum == 'CAPSCAN_stool').sum())),
        ('P68 legend CapScan intestinal-fluid acquisitions', int((x.stratum == 'CAPSCAN_intestinal').sum())),
        ('L127/P68 legend FL-only spectra (millions)', audit['total'] / 1e6),
        ('L127/P68 legend patients off and on antibiotics', int(pr.loc['Human', 'n_patients'])),
        ('L127/P68 legend whisker length (IQR; design constant BOX_WHIS, not computed)', BOX_WHIS),
        ('L127/P68 legend heatmap patients', n_heatmap),
        ('L127/P68 legend cross-omics patients', int(native.n.iloc[0])),
        ('L129/P69 same-spectrum stool gain 9-50 (%)', text['same_gain_median_9_stool']),
        ('L129 stool acquisitions', text['same_n_stool']),
        ('L129/P69 lowest region same-spectrum gain 9-50 (%)', min(text[f'same_gain_median_9_{s}'] for s in fluid_sites)),
        ('L129/P69 highest region same-spectrum gain 9-50 (%)', max(text[f'same_gain_median_9_{s}'] for s in fluid_sites)),
        ('L129 fewest acquisitions in a region', min(text[f'same_n_{s}'] for s in fluid_sites)),
        ('L129 most acquisitions in a region', max(text[f'same_n_{s}'] for s in fluid_sites)),
        ('L129/P69 stool acquisitions with more FL peptides 9-50', text['same_n_more_9_stool']),
        ('L129/P69 stool acquisitions (denominator)', text['same_n_stool']),
        ('L129/P69 fluid acquisitions with more FL peptides 9-50', text['same_n_more_9_fluid']),
        ('L129/P69 fluid acquisitions', text['n_fluid']),
        ('L129/P69 mixed model fluid/stool fold 9-50', mm['fold']),
        ('L129/P69 mixed model CI low', mm['fold_ci_low']),
        ('L129/P69 mixed model CI high', mm['fold_ci_high']),
        ('L129/P69 mixed model P (unadjusted Wald z)', mm['p']),
        ('L129/P69 region omnibus P (3-df Wald)', mm['region_p']),
        ('L129/P70 acquisitions with recorded DNA', text['n_dna_recorded']),
        ('L129 matched acquisitions (DNA denominator)', text['n_dna_table']),
        ('L129/P70 lowest regional median DNA (ng)', min(text[f'dna_ng_median_{s}'] for s in fluid_sites)),
        ('L129/P70 highest regional median DNA (ng)', max(text[f'dna_ng_median_{s}'] for s in fluid_sites)),
        ('L129/P70 stool median DNA (ng)', text['dna_ng_median_stool']),
        ('L129/P70 fluid median MG database (thousand proteins)', text['db_median_fluid'] / 1e3),
        ('L129/P70 stool median MG database (thousand proteins)', text['db_median_stool'] / 1e3),
        ('L129/P70 participants with smaller fluid database (of those with both)', text['participants_fluid_db_smaller']),
        ('L129 stool-only participant', ','.join(text['participants_stool_only'])),
        ('L129/P70 fluid/stool median yield 9-50, MG (%)', ss['yield_fluid_over_stool_MG_pct']),
        ('L129/P70 fluid/stool median yield 9-50, FL (%)', ss['yield_fluid_over_stool_FL_pct']),
        ('L131/P71 FL-only spectra (millions)', audit['total'] / 1e6),
        ('L131/P71 sequence absent from MG database (%)', audit['percent']['A_absent']),
        ('L131/P71 present, top MG target failed (%)', audit['percent']['D_displaced'] + audit['percent']['C_top_fails_q']),
        ('P71 different MG peptide not accepted (%)', audit['percent']['D_displaced']),
        ('P71 same MG peptide not accepted (%)', audit['percent']['C_top_fails_q']),
        ('L131 reported candidates per spectrum (Sage setting, not computed)', ''),
        ('L133/P70 PREDICT same-spectrum > cross-fitted 9-50', svc['PREDICT'][0]),
        ('L133/P70 PREDICT acquisitions', svc['PREDICT'][1]),
        ('L133/P70 CapScan same-spectrum > cross-fitted 9-50', svc['CAPSCAN'][0]),
        ('L133/P70 CapScan acquisitions', svc['CAPSCAN'][1]),
        ('L133/P65 PREDICT fewer 9-50 peptides (%)', -bt('PREDICT', L_BC, 'gain_cf_pooled', 'median')),
        ('L133/P65 CapScan stool fewer 9-50 peptides (%)', -bt('CAPSCAN_stool', L_BC, 'gain_cf_pooled', 'median')),
        ('L133 intestinal fluid more peptides (%)', bt('CAPSCAN_intestinal', L_BC, 'gain_cf_pooled', 'median')),
        ('L137/P73 paired patients', int(pr.loc['Human', 'n_patients'])),
        ('L137/P73 human share on-off (pp)', pr.loc['Human', 'mean_on_minus_off']),
        ('L137/P73 human share CI low', pr.loc['Human', 'CI95_low']),
        ('L137/P73 human share CI high', pr.loc['Human', 'CI95_high']),
        ('L137/P73 human share BH q', pr.loc['Human', 'q_two_clinical_endpoints']),
        ('L137/P73 fewer protein groups', -pr.loc['pg_parsimony', 'mean_on_minus_off']),
        ('L137/P73 protein groups CI (smaller decrease)', -pr.loc['pg_parsimony', 'CI95_high']),
        ('L137/P73 protein groups CI (larger decrease)', -pr.loc['pg_parsimony', 'CI95_low']),
        ('L137/P73 protein groups BH q', pr.loc['pg_parsimony', 'q_two_clinical_endpoints']),
        ('L139/P74 cross-omics patients', int(native.n.iloc[0])),
        ('L139/P74 rho protein-DNA', nat.loc[('DNA', 'Non-host protein'), 'rho']),
        ('L139/P74 rho protein-RNA', nat.loc[('RNA', 'Non-host protein'), 'rho']),
        ('L139/P74 rho protein-metabolites', nat.loc[('Non-host protein', 'Metabolites'), 'rho']),
        ('L139/P74 largest q of the three', max(nat.loc[k, 'q'] for k in [('DNA', 'Non-host protein'), ('RNA', 'Non-host protein'), ('Non-host protein', 'Metabolites')])),
        ('L139/P74 pathway RNA-protein rho', fun.loc[('Pathway', 'metaT', 'microbial_metaP'), 'rho']),
        ('L139/P74 pathway RNA-protein q', fun.loc[('Pathway', 'metaT', 'microbial_metaP'), 'q']),
    ]
    return [(k, v.item() if hasattr(v, 'item') else v) for k, v in rows]


# ----------------------------------------------------------------------------------
# Display and checks
# ----------------------------------------------------------------------------------
#: Display crop of the b,c page (points): the region compose4.py replaced on the composite.
BC_CROP_PT = (20.0, 70.0, 330.0, 215.0)
#: Height of panel a on the composite (points); artwork/panel_a.pdf is this band.
PANEL_A_HEIGHT_PT = 77.0


def render_png(pdf, png, dpi=300, clip=None):
    """Rasterise page 1 of a PDF.

    Parameters
    ----------
    pdf, png : pathlib.Path
        Input PDF and output PNG.
    dpi : int
        Resolution.
    clip : tuple of float or None
        Optional clip rectangle (x0, y0, x1, y1) in points.

    Returns
    -------
    pathlib.Path
        ``png``.
    """
    import fitz
    with fitz.open(pdf) as d:
        r = fitz.Rect(*clip) if clip else None
        d[0].get_pixmap(dpi=dpi, clip=r, alpha=False).save(png)
    return png


def held_file(unit_dir, sub, name):
    """Path of a file that contains the controlled panels e, g, h (published composites,
    panel-a artwork). Read from ``unit_dir/sub/name`` if present, else from
    ``$FASTALAKE_CONTROLLED_DIR/Fig4_expected/name`` (composites) or ``.../Fig4_artwork/name``;
    None when neither exists. These files are not in the public repository because they
    show the per-patient antibiotic panels."""
    p = Path(unit_dir) / sub / name
    if p.is_file():
        return p
    root = os.environ.get("FASTALAKE_CONTROLLED_DIR")
    if root:
        for d in ("Fig4_expected", "Fig4_artwork"):
            q = Path(root) / d / name
            if q.is_file():
                return q
    return None



def artwork_preview(composite_png, png):
    """Crop panel a from the published composite PNG (300 dpi) for display.

    The vector artwork renders slowly (it carries the whole page's clipped content), so the
    notebook shows the published raster instead; artwork/panel_a.pdf is the deliverable.

    Parameters
    ----------
    composite_png, png : pathlib.Path
        Published composite raster and output PNG.

    Returns
    -------
    pathlib.Path
        ``png``.
    """
    from PIL import Image as _Im
    im = _Im.open(composite_png)
    im.crop((0, 0, im.width, int(round(PANEL_A_HEIGHT_PT * 300 / 72)))).save(png)
    return png


#: Pixel checks: (name, regenerated file in OUT, expected file, needs controlled input).
PIXEL_CHECKS = [("panel_bc", "panel_bc.pdf", "panel_bc_v42.pdf", False),
                ("panel_d", "panel_d.pdf", "panel_d.pdf", False),
                ("panel_e", "panel_e.pdf", "panel_e.pdf", True),
                ("panel_f", "panel_f.pdf", "panel_f.pdf", False),
                ("panel_g", "panel_g.pdf", "panel_g.pdf", True),
                ("panel_h", "panel_h.pdf", "panel_h.pdf", True),
                ("panel_i", "panel_i.pdf", "panel_i.pdf", False),
                ("panel_i_colourbar", "panel_i_colourbar.pdf", "panel_i_colourbar.pdf", False)]


def run_pixel_checks(out_dir, unit_dir, dpi=300):
    """Compare every regenerated panel with the original producer's output at ``dpi``.

    Public panels: reference ``expected/<file>``, receipts ``checks/<name>.json`` and
    ``checks/<name>.diff.png``. Panels drawn from controlled inputs (e, g, h): their
    producer renders are not shipped in the repository (OPEN_ITEMS AUTHOR); the reference
    is read from ``$FASTALAKE_CONTROLLED_DIR/Fig4_expected/<file>`` and the receipts are
    written to ``<out_dir>/checks_controlled/`` (build/, not committed). Skips all checks
    when the Söhne fonts are absent. Ends with a negative control (the regenerated panel-i
    colour-bar page against the expected b,c page, same canvas size), which must not be
    IDENTICAL.

    Parameters
    ----------
    out_dir : pathlib.Path
        Directory with the regenerated panels (build/Fig4).
    unit_dir : pathlib.Path
        This unit's directory.
    dpi : int
        Comparison resolution.

    Returns
    -------
    list of dict
        One summary per check.
    """
    from tools.check_pixels import compare
    if not fonts_present():
        print("pixel checks SKIPPED: fonts absent")
        return []
    res = []
    for name, new, exp, ctl in PIXEL_CHECKS:
        if not (out_dir / new).is_file():
            print(f"{name}: SKIPPED (not rendered{'; controlled input absent' if ctl else ''})")
            continue
        if ctl:
            root = os.environ.get("FASTALAKE_CONTROLLED_DIR")
            ref = Path(root) / "Fig4_expected" / exp if root else None
            if ref is None or not ref.is_file():
                print(f"{name}: SKIPPED (producer render is held with the controlled inputs, not in the repository)")
                continue
            rec_dir = out_dir / "checks_controlled"; rec_dir.mkdir(exist_ok=True)
            r = compare(out_dir / new, ref, rec_dir / name, dpi=dpi)
        else:
            r = compare(out_dir / new, unit_dir / "expected" / exp, unit_dir / "checks" / name, dpi=dpi)
        res.append(dict(check=name, verdict=r["verdict"], differing_pixels=r.get("differing_pixels"),
                        size_px=r["new_size_px"]))
    if (out_dir / "panel_i_colourbar.pdf").is_file():
        import tempfile
        with tempfile.TemporaryDirectory() as td:
            r = compare(out_dir / "panel_i_colourbar.pdf", unit_dir / "expected" / "panel_bc.pdf", Path(td) / "neg", dpi=dpi)
        assert r["verdict"] != "IDENTICAL", "negative control failed: the check cannot see a different panel"
        res.append(dict(check="negative control (i colour bar vs expected b,c)", verdict=r["verdict"],
                        differing_pixels=r.get("differing_pixels"), size_px=r["new_size_px"]))
    return res


#: Box of the composite replaced by the cross-fitted b,c page (points; fig4bc_crossfit compose4.py).
BC_V30_BOX_PT = (20.0, 86.5, 321.0, 202.5)
#: Band of the composite where fig_v27 placed panel h at scale 1 (points; build_fig4_v27.py H_TOP, H_TOP + HH).
H_V30_BAND_PT = (0.0, 331.0, PAGE_W, 549.0)


def _crop_px(arr, box_pt, dpi):
    """Crop an image array to a box given in points (floor/ceil to whole pixels)."""
    s = dpi / 72
    x0, y0, x1, y1 = int(np.floor(box_pt[0] * s)), int(np.floor(box_pt[1] * s)), int(np.ceil(box_pt[2] * s)), int(np.ceil(box_pt[3] * s))
    return arr[y0:y1, x0:x1]


def _page_rgb(pdf, dpi=300, place=None):
    """Render a PDF page as an RGB array; ``place`` puts it on a blank composite-size page first."""
    import fitz
    if place is not None:
        src = fitz.open(pdf); doc = fitz.open(); pg = doc.new_page(width=PAGE_W, height=PAGE_H_BC)
        pg.show_pdf_page(fitz.Rect(*place), src, 0); page = pg
    else:
        doc = fitz.open(pdf); page = doc[0]
    pix = page.get_pixmap(dpi=dpi, alpha=False)
    return np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, 3).copy()


def check_against_v30(out_dir, unit_dir, dpi=300):
    """Pixel-check regenerated panels against the v30 composite ``expected/Fig4_v30.png``.

    * b,c: the regenerated b,c page (composite size) rendered at ``dpi`` and cropped to
      BC_V30_BOX_PT, after placing it on a blank page as compose4.py did (show_pdf_page),
      against the same crop of the v30 composite. Receipt
      ``checks/v30_panel_bc``.
    * h (controlled): the regenerated panel h placed on a composite-size page at
      H_V30_BAND_PT, cropped to that band, against the same band of the v30 composite.
      Receipt in ``<out_dir>/checks_controlled/v30_panel_h`` (not committed).
    * composite: ``expected/Fig4_v30.png`` against ``expected/Fig4_v28f.png``; documents
      that v30 differs from v28f only in the panel-h track label. Receipt
      ``checks/v30_vs_v28f``.

    Parameters
    ----------
    out_dir, unit_dir : pathlib.Path
        build/Fig4 and this unit's directory.
    dpi : int
        Resolution of the v30 composite PNG (300).

    Returns
    -------
    list of dict
        One summary per check.
    """
    from PIL import Image as _Im
    from tools.check_pixels import compare
    if not fonts_present():
        print("v30 pixel checks SKIPPED: fonts absent")
        return []
    v30p, v28p = held_file(unit_dir, "expected", "Fig4_v30.png"), held_file(unit_dir, "expected", "Fig4_v28f.png")
    # newest held composite is preferred: v36 (b labels give the median paired change), then v33
    # (b,c at 295 CapScan acquisitions), then v30; all have the same layout and size
    refp, tag = v30p, "v30"
    for t in ("v33", "v36", "v42"):
        q = held_file(unit_dir, "expected", f"Fig4_{t}.png")
        if q is not None:
            refp, tag = q, t
    if refp is None:
        print("composite checks SKIPPED: the published composite shows controlled panels e, g, h and is held with the controlled inputs")
        return []
    v30 = np.asarray(_Im.open(refp).convert("RGB"))
    res = []
    jobs = [(f"{tag}_panel_bc", out_dir / "panel_bc.pdf", (0.0, 0.0, PAGE_W, PAGE_H_BC), BC_V30_BOX_PT, unit_dir / "checks")]
    if (out_dir / "panel_h.pdf").is_file():
        cdir = out_dir / "checks_controlled"; cdir.mkdir(exist_ok=True)
        jobs.append((f"{tag}_panel_h", out_dir / "panel_h.pdf", H_V30_BAND_PT, H_V30_BAND_PT, cdir))
    else:
        print(f"{tag}_panel_h: SKIPPED (controlled input absent)")
    for name, pdf, place, box, rdir in jobs:
        new = _crop_px(_page_rgb(pdf, dpi, place), box, dpi)
        exp = _crop_px(v30, box, dpi)
        _Im.fromarray(new).save(out_dir / f"{name}_new.png"); _Im.fromarray(exp).save(out_dir / f"{name}_expected.png")
        r = compare(out_dir / f"{name}_new.png", out_dir / f"{name}_expected.png", rdir / name, dpi=dpi)
        res.append(dict(check=name, verdict=r["verdict"], differing_pixels=r.get("differing_pixels"),
                        diff_bbox_px=r.get("diff_bbox_px"), size_px=r["new_size_px"]))
    if v30p is None or v28p is None:
        return res
    r = compare(v30p, v28p, unit_dir / "checks" / "v30_vs_v28f", dpi=dpi)
    res.append(dict(check="v30_vs_v28f (composites)", verdict=r["verdict"], differing_pixels=r.get("differing_pixels"),
                    diff_bbox_px=r.get("diff_bbox_px"), size_px=r["new_size_px"]))
    return res


def extract_panel_a_artwork(composite_pdf, out_pdf):
    """Write the panel-a band (top PANEL_A_HEIGHT_PT points) of the composite as a vector PDF.

    This is how ``artwork/panel_a.pdf`` was made from ``expected/Fig4_v28f.pdf``. Not run by
    the notebook: rendering the composite is slow (several minutes) because the page
    carries the clipped content of earlier assembly steps.

    Parameters
    ----------
    composite_pdf, out_pdf : pathlib.Path
        Composite PDF and output PDF.

    Returns
    -------
    pathlib.Path
        ``out_pdf``.
    """
    import fitz
    src = fitz.open(composite_pdf)
    W = src[0].rect.width
    clip = fitz.Rect(0, 0, W, PANEL_A_HEIGHT_PT)
    out = fitz.open(); pg = out.new_page(width=W, height=clip.height)
    pg.show_pdf_page(pg.rect, src, 0, clip=clip)
    out.save(out_pdf, garbage=4, deflate=True)
    return out_pdf
