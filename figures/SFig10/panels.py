"""Supplementary Figure 10: PREDICT metabolite-enzyme pairs ordered by Enterotype Dysbiosis Score.

Panel
-----
(single panel, no letter) Three metabolite-enzyme (KO) pairs selected post hoc, 85
patients in EDS order, with EDS, human-signal and antibiotic tracks above the heatmap.

The v30 figure keeps only the former panel a; the former panel b (the 19-feature overview,
identical to Fig. 4h) was removed. ``render`` ports the v30 build script
(build_sfig10.py, see PROVENANCE.md): the drawing functions ``tracks``, ``matrix`` and
``colorbar`` are the earlier producer's (build_suppfig10_merged.py), placed on a
180 x 62 mm canvas in the repository fl_style at 6 pt; layout and colours unchanged.

Statistics (``stat_*``) port the computing code of the upstream analysis
(build_approved_eds_cr003.py and the screen definitions ``corr`` and ``bh`` of
analyze.py), unchanged in logic. The 20-feature overview tables are kept because the
Methods and Supplementary Table 4 still print their counts and ``stat_clip`` reads the
producer's plotted values of both former panels.

Source data (``source_data/SFig10/``): patients.tsv, focused_values.tsv,
focused_features.tsv, legacy_85_values.tsv, legacy_features.tsv, legacy_259_values.tsv,
augmented_pair_screen.tsv, all_EDS_associations.tsv, selected_pair_summary.tsv,
separate_adjustments_all_models.tsv, producer_plotted_ab.tsv. Antibiotic state and
clinical group are controlled inputs (``FASTALAKE_CONTROLLED_DIR``, file
``CONTROLLED_FILE``). The rendered figure and the v30 image it is checked against show
the per-patient antibiotic track, so they are also kept in the controlled directory
(``EXPECTED_PNG``, ``EXPECTED_PDF``) and never written into the repository.
"""

import hashlib
import io
import itertools
import os
from pathlib import Path

import matplotlib
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from scipy.stats import rankdata, spearmanr, t

# ---- choices (named constants) ------------------------------------------------------
N_PATIENTS = 85                 # v30 paragraph 409: earliest matched observation of 85 patients
CLIP_Z = 2                      # v30 paragraphs 145, 239: standardized within feature, clipped at ±2
MIN_POSITIVE_N = 43             # v30 paragraph 409: at least 43 patients for association tests
Q_CUTOFF = 0.05                 # v30 paragraph 409: candidate pairs pass q<0.05
PAIRS = [                       # v30 paragraph 409: the three post hoc pairs (metabolite, KO, labels, reaction)
    ("HILIC_FT23003", "K01812", "Glucuronic acid", "Glucuronate isomerase", "R01482"),
    ("HILIC_FT24070", "K01667", "Tryptophan", "Tryptophanase", "R00673"),
    ("C18_F2423", "K01442", "Glycochenodeoxycholic acid", "Bile salt hydrolase", "R03975"),
]
COVARIATES = ["host", "depth", "antibiotics", "clinical"]   # v30 paragraph 409: 16 combinations
NORMALIZATIONS = ["native", "annotated_relative"]          # v30 paragraph 409
SCOPES = ["pair_complete", "common_84"]                    # v30 paragraph 409: pair-specific and common subsets
JOINT_MODEL = "host_depth_antibiotic_clinical"             # v30 paragraph 409: joint adjustment
SEPARATORS_A = [1.5, 3.5]       # producer: lines between the three pairs (panel a)
ORGANIC_ACID_LABEL = ("Dimethylsuccinic acid", "Dimethylsuccinic / methylglutaric / adipic acid")  # former panel b label, used by stat_clip
FIG_MM = (180.0, 62.0)          # v30 build: authored at the placed size, 180 x 62 mm
FONT_PT = 6.0                   # v30 build: 6 pt text
HEATMAP_LR = (0.30, 0.905)      # v30 build: left and right edge of tracks and heatmap (figure fraction)
TRACK_TOP_HEIGHT = (0.975, 0.055)          # v30 build: tracks top and height per track
HEATMAP_RECT_Y = (0.19, 0.60)              # v30 build: heatmap bottom and height
COLORBAR_RECT = [0.925, 0.25, 0.010, 0.48] # v30 build: colour bar axes
PNG_DPI = 300                   # v30 build: PNG resolution
CONTROLLED_FILE = "patient_order_and_covariates.tsv"   # controlled/CONTROLLED_INPUTS.md
EXPECTED_PNG = "SFig10_v33.png"   # v33 figure PNG, controlled (shows the antibiotic track)
EXPECTED_PDF = "SFig10_v33.pdf"   # v33 figure PDF, controlled (shows the antibiotic track)
WORDS = {3: "three", 6: "six", 19: "nineteen"}   # numbers printed as words in the manuscript
V30_PLOTTED_SHA256 = "de6c8bb6eef41fec138f0253497f0175c4bc785621c54033b73a7c50a9800d62"  # v30 build plotted_values.tsv

CMAP = LinearSegmentedColormap.from_list("fl", ["#0072B2", "white", "#E69F00"])
CMAP.set_bad("#dedede")


# ---- loading ------------------------------------------------------------------------
def _read(path):
    """Read one source table with exact (round-trip) float parsing.

    Parameters
    ----------
    path : pathlib.Path
        Tab-separated file.

    Returns
    -------
    pandas.DataFrame
        The table.
    """
    return pd.read_csv(path, sep="\t", float_precision="round_trip")


def load(sd):
    """Read every SFig10 source table.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/SFig10``.

    Returns
    -------
    dict
        ``patients`` (85 rows, EDS order), ``focused_features``, ``focused_raw`` and
        ``focused_z`` (6 x 85), ``legacy_features``, ``legacy_raw`` and ``legacy_z``
        (20 x 85), ``legacy_259`` (20 x 259 raw), ``screen``, ``eds_assoc``,
        ``summary``, ``adjustments``, ``producer_plotted`` (the producer's plotted
        values for the heatmap cells of panels a and b).
    """
    d = {"patients": _read(sd / "patients.tsv")}
    assert list(d["patients"].patient_position) == list(range(1, N_PATIENTS + 1))
    for name, feats in [("focused", "focused_features.tsv"), ("legacy", "legacy_features.tsv")]:
        f = _read(sd / feats)
        v = _read(sd / f"{'focused_values' if name == 'focused' else 'legacy_85_values'}.tsv")
        col = "feature" if name == "focused" else "feature_id"
        order = list(f[col])
        assert list(dict.fromkeys(v.feature)) == order
        for k in ["raw", "z"]:
            d[f"{name}_{k}"] = v.pivot(index="feature", columns="patient_position", values=k).loc[order].to_numpy()
        d[f"{name}_features"] = f
    v = _read(sd / "legacy_259_values.tsv")
    d["legacy_259"] = v.pivot(index="feature", columns="observation_index", values="raw").loc[list(d["legacy_features"].feature_id)].to_numpy()
    d["screen"] = _read(sd / "augmented_pair_screen.tsv")
    d["eds_assoc"] = _read(sd / "all_EDS_associations.tsv")
    d["summary"] = _read(sd / "selected_pair_summary.tsv")
    d["adjustments"] = _read(sd / "separate_adjustments_all_models.tsv")
    d["producer_plotted"] = _read(sd / "producer_plotted_ab.tsv")
    return d


def load_controlled(patients):
    """Read antibiotic state and clinical group from the controlled input, or return None.

    The file is the producer's per-patient table (one row per patient, EDS order). Rows
    are matched by position and checked against the public EDS, human share and depth.

    Parameters
    ----------
    patients : pandas.DataFrame
        Public per-patient table.

    Returns
    -------
    pandas.DataFrame or None
        Columns ``antibiotic_state`` and ``clinical_group`` in display order; None (with a
        printed SKIPPED message) when the controlled input is absent.
    """
    root = os.environ.get("FASTALAKE_CONTROLLED_DIR")
    path = Path(root) / CONTROLLED_FILE if root else None
    if path is None or not path.is_file():
        print(f"SKIPPED: needs controlled input {CONTROLLED_FILE} (see controlled/CONTROLLED_INPUTS.md)")
        return None
    c = pd.read_csv(path, sep="\t", float_precision="round_trip")
    assert len(c) == N_PATIENTS
    for col in ["dysbiosis_score", "Human", "pg_parsimony"]:
        assert np.array_equal(c[col].to_numpy(), patients[col].to_numpy()), col
    return c[["antibiotic_state", "clinical_group"]].reset_index(drop=True)


# ---- producer plot code (moved as is) ---------------------------------------------
def tracks(fig, order, abx, left, right, top, height, fontsize):
    """Draw the EDS, human-signal and antibiotic tracks (producer ``tracks``).

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Canvas.
    order : pandas.DataFrame
        Public per-patient table in display order (EDS and human share).
    abx : pandas.DataFrame
        Output of ``load_controlled`` (antibiotic state).
    left, right, top, height : float
        Track geometry in figure fractions.
    fontsize : float
        Label size in points.

    Returns
    -------
    list of tuple
        (track label, plotted values) for each track.
    """
    tt = [(order.dysbiosis_score.to_numpy(), 'EDS (0–1)', 'Greys', 1),
          (100 * order.Human.to_numpy(), 'Human-assigned LFQ share (%)', LinearSegmentedColormap.from_list('human', ['white', '#CC79A7']), 100),
          (abx.antibiotic_state.map({'Off': 0, 'On': 1}).to_numpy(), 'Antibiotics (dark: on)', ListedColormap(['#ededed', '#333333']), 1)]
    out = []
    for i, (v, l, c, m) in enumerate(tt):
        ax = fig.add_axes([left, top - (i + 1) * height, right - left, height * .72])
        ax.imshow(v[None, :], cmap=c, vmin=0, vmax=m, aspect='auto', interpolation='nearest')
        ax.set(xticks=[], yticks=[0], yticklabels=[l]); ax.tick_params(length=0, labelsize=fontsize)
        for sp in ax.spines.values(): sp.set_visible(False)
        out.append((l, v))
    return out


def matrix(fig, values, labs, rect, seps, fontsize, xticks=True):
    """Draw one clipped z-score heatmap (producer ``matrix``).

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Canvas.
    values : numpy.ndarray
        Features x patients z scores; clipped at +-``CLIP_Z`` for display.
    labs : list of str
        Row labels.
    rect : list of float
        Axes rectangle in figure fractions.
    seps : list of float
        Row positions of separator lines.
    fontsize : float
        Text size in points.
    xticks : bool
        Draw patient-position ticks and the x label.

    Returns
    -------
    matplotlib.image.AxesImage
        The heatmap image.
    """
    ax = fig.add_axes(rect)
    im = ax.imshow(np.clip(values, -CLIP_Z, CLIP_Z), cmap=CMAP, vmin=-CLIP_Z, vmax=CLIP_Z, aspect='auto', interpolation='nearest')
    ax.set_yticks(range(len(labs)), labs); ax.tick_params(axis='y', length=0, pad=3, labelsize=fontsize)
    if xticks:
        ax.set_xticks([0, 21, 42, 63, 84], ['1', '22', '43', '64', '85'])
        ax.set_xlabel('85 patients ordered by metagenomic EDS', labelpad=3, fontsize=fontsize)
        ax.tick_params(axis='x', labelsize=fontsize)
    else: ax.set_xticks([])
    for sep in seps: ax.axhline(sep, color='#666666', lw=.6)
    for sp in ax.spines.values(): sp.set_visible(False)
    return im


def colorbar(fig, im, rect, fontsize):
    """Draw the shared z-score colour bar (producer ``colorbar``).

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        Canvas.
    im : matplotlib.image.AxesImage
        Image whose colour map and limits the bar shows.
    rect : list of float
        Axes rectangle in figure fractions.
    fontsize : float
        Text size in points.

    Returns
    -------
    None
    """
    cb = fig.colorbar(im, cax=fig.add_axes(rect), ticks=[-2, 0, 2])
    cb.set_label('Within-feature z score', labelpad=2, fontsize=fontsize)
    cb.ax.tick_params(length=2, pad=1, labelsize=fontsize)


def legacy_rows(d):
    """Return the 19 plotted overview rows (features with any positive value) and labels.

    Parameters
    ----------
    d : dict
        Output of ``load``.

    Returns
    -------
    tuple
        (row indices into the 20 overview features, display labels).
    """
    keep = np.flatnonzero(d["legacy_features"].n_positive_85 > 0)
    labels = [s.replace(*ORGANIC_ACID_LABEL) for s in d["legacy_features"].iloc[keep].label]
    return keep, labels


def render(d, abx, out_stem):
    """Render Supplementary Fig. 10 (v30) as the v30 build script did.

    Ported from build_sfig10.py: one 180 x 62 mm canvas, the three tracks, the six-row
    focused heatmap with patient-position ticks, and the colour bar. Writes
    ``<out_stem>.png`` (``PNG_DPI``) and ``<out_stem>.pdf`` with ``fl.save``.

    Parameters
    ----------
    d : dict
        Output of ``load``.
    abx : pandas.DataFrame
        Output of ``load_controlled`` (antibiotic track).
    out_stem : pathlib.Path
        Output path without extension.

    Returns
    -------
    dict
        ``cells_identical`` (the 510 drawn heatmap cells equal the producer's panel-a
        plotted values exactly, as asserted by the v30 build), ``plotted`` (the v30 plotted
        values table) and ``plotted_sha256`` (its sha256 as the v30 build wrote it).
    """
    import matplotlib.pyplot as plt
    from style import fl_style as fl

    fl.init()
    fig = plt.figure(figsize=fl.mm(*FIG_MM))
    L, R = HEATMAP_LR
    tr = tracks(fig, d["patients"], abx, L, R, TRACK_TOP_HEIGHT[0], TRACK_TOP_HEIGHT[1], FONT_PT)
    im = matrix(fig, d["focused_z"], list(d["focused_features"].label),
                [L, HEATMAP_RECT_Y[0], R - L, HEATMAP_RECT_Y[1]], SEPARATORS_A, FONT_PT, True)
    colorbar(fig, im, COLORBAR_RECT, FONT_PT)
    pa = d["producer_plotted"]
    pa = pa[pa.panel.eq("a")]
    mine = np.asarray(im.get_array(), dtype=float).ravel(order="C")
    exp = pa.sort_values(["row", "patient_position"]).plotted.to_numpy(float)
    same = mine.shape == exp.shape and np.array_equal(np.nan_to_num(mine, nan=-99), np.nan_to_num(exp, nan=-99))
    rows = []
    for lab, v in tr:
        for c, x in enumerate(v):
            rows.append(dict(panel="track", label=lab, patient_position=c + 1, plotted=x))
    plotted = pd.concat([pa.drop(columns="panel"), pd.DataFrame(rows)])
    buf = io.StringIO()
    plotted.to_csv(buf, sep="\t", index=False)
    fl.save(fig, out_stem, dpi=PNG_DPI, formats=("pdf", "png"))
    plt.close(fig)
    return dict(cells_identical=bool(same), plotted=plotted,
                plotted_sha256=hashlib.sha256(buf.getvalue().encode()).hexdigest())


def plotted_heatmap_rows(d):
    """Rebuild the heatmap rows of the producer's plotted-values table (panels a and b).

    Parameters
    ----------
    d : dict
        Output of ``load``.

    Returns
    -------
    list of dict
        One dict per heatmap cell: panel, row, label, patient_position, z, plotted (z
        clipped at +-2; NaN for non-detections), in the producer's row order.
    """
    zz, labels = d["focused_z"], list(d["focused_features"].label)
    keep, oldlabels = legacy_rows(d)
    rows = []
    for r, lab in enumerate(labels):
        for c in range(zz.shape[1]):
            rows.append(dict(panel='a', row=r + 1, label=lab, patient_position=c + 1, z=zz[r, c], plotted=float(np.clip(zz[r, c], -2, 2)) if np.isfinite(zz[r, c]) else np.nan))
    for r, lab in enumerate(oldlabels):
        for c in range(d["legacy_z"].shape[1]):
            v = d["legacy_z"][keep][r, c]
            rows.append(dict(panel='b', row=r + 1, label=lab, patient_position=c + 1, z=v, plotted=float(np.clip(v, -2, 2)) if np.isfinite(v) else np.nan))
    return rows


def plotted_values(d, tr):
    """Rebuild the producer's plotted-values table (panel, row, label, position, z, plotted).

    Parameters
    ----------
    d : dict
        Output of ``load``.
    tr : list of tuple
        Output of ``tracks`` (needs the controlled input).

    Returns
    -------
    pandas.DataFrame
        Heatmap rows of panels a and b followed by the track rows, as the producer wrote it.
    """
    rows = plotted_heatmap_rows(d)
    for lab, v in tr:
        for c, x in enumerate(v):
            rows.append(dict(panel='a_track', row='', label=lab, patient_position=c + 1, z='', plotted=x))
    return pd.DataFrame(rows)


# ---- statistics -----------------------------------------------------------------------
def bh(p):
    """Benjamini-Hochberg q values (screen definition ``bh`` of analyze.py); NaN stays NaN.

    Parameters
    ----------
    p : array_like
        P values of one family; non-finite entries are not tested.

    Returns
    -------
    numpy.ndarray
        q values, same length as ``p``; the family size is the number of finite P values.
    """
    p = np.asarray(p, float); v = np.full(len(p), np.nan); idx = np.where(np.isfinite(p))[0]; s = idx[np.argsort(p[idx])]
    v[s] = np.minimum.accumulate((p[s] * len(s) / np.arange(1, len(s) + 1))[::-1])[::-1].clip(0, 1); return v


def corr(x, y, c=None):
    """Spearman or partial rank correlation with a two-sided t test (screen ``corr``).

    Ranks x and y on the complete subset, regresses both on [1, covariates] when given,
    correlates the residuals and tests with df = n - 2 - k (k = covariate rank).
    Statistical unit: patient.

    Parameters
    ----------
    x, y : numpy.ndarray
        Values per patient; non-finite entries are excluded.
    c : numpy.ndarray or None
        Patients x covariates matrix, or None for the plain Spearman correlation.

    Returns
    -------
    tuple
        (n, rho, P, df); rho and P are NaN when n < ``MIN_POSITIVE_N``.
    """
    ok = np.isfinite(x) & np.isfinite(y)
    if c is not None: ok &= np.isfinite(c).all(axis=1)
    n = int(ok.sum())
    if n < MIN_POSITIVE_N or np.unique(x[ok]).size < 3 or np.unique(y[ok]).size < 3: return n, np.nan, np.nan, np.nan
    a = rankdata(x[ok]); b = rankdata(y[ok]); k = 0
    if c is not None:
        C = np.column_stack([np.ones(n), c[ok]]); k = np.linalg.matrix_rank(C) - 1
        a = a - C @ np.linalg.lstsq(C, a, rcond=None)[0]; b = b - C @ np.linalg.lstsq(C, b, rcond=None)[0]
    r = np.corrcoef(a, b)[0, 1]; df = n - k - 2
    return n, float(r), float(2 * t.sf(abs(r) * np.sqrt(df / max(1e-16, 1 - r * r)), df)), df


def _pair_arrays(d, m, k):
    """Return the raw KO and metabolite values of one pair.

    Parameters
    ----------
    d : dict
        Output of ``load``.
    m, k : str
        Metabolite feature id and KO id.

    Returns
    -------
    tuple of numpy.ndarray
        (KO values, metabolite values), 85 patients each, in display order.
    """
    feats = list(d["focused_features"].feature)
    return d["focused_raw"][feats.index(k)], d["focused_raw"][feats.index(m)]


def _positive(v):
    """Keep positive finite values and set every other value to NaN.

    Parameters
    ----------
    v : numpy.ndarray
        Values.

    Returns
    -------
    numpy.ndarray
        Same shape; non-detections and non-positive values are NaN.
    """
    return np.where(np.isfinite(v) & (v > 0), v, np.nan)


def stat_a_unadjusted(d):
    """Spearman correlation of each displayed metabolite with its KO signal.

    Statistical unit: patient (earliest matched observation); n = 85, 84, 85 patients
    with positive finite values for both features. Two-sided t test with df = n - 2.

    Parameters
    ----------
    d : dict
        Output of ``load``.

    Returns
    -------
    pandas.DataFrame
        One row per pair: n, rho, t statistic, df, P (screen ``corr``), scipy rho and P.
    """
    rows = []
    for m, k, ml, kl, _ in PAIRS:
        x, y = _pair_arrays(d, m, k)
        n, r, p, df = corr(_positive(x), _positive(y))
        ok = np.isfinite(_positive(x)) & np.isfinite(_positive(y))
        s = spearmanr(x[ok], y[ok])
        rows.append(dict(metabolite=m, KO=k, label=f"{ml} - {kl}", n=n, rho=r, t=r * np.sqrt(df / (1 - r * r)), df=df,
                         P=p, scipy_rho=float(s.statistic), scipy_P=float(s.pvalue)))
    return pd.DataFrame(rows)


def _blocks(order, abx, idx):
    """Covariate blocks ranked within the analysis subset (build_approved ``blocks``).

    Parameters
    ----------
    order : pandas.DataFrame
        Public per-patient table (human share and protein-group count).
    abx : pandas.DataFrame or None
        Output of ``load_controlled``; without it only host and depth blocks exist.
    idx : numpy.ndarray
        Patient indices of the analysis subset.

    Returns
    -------
    dict
        Covariate name to a patients x columns array.
    """
    b = {'host': rankdata(order.Human.to_numpy()[idx])[:, None],
         'depth': rankdata(np.log2(order.pg_parsimony.to_numpy()[idx]))[:, None]}
    if abx is not None:
        b['antibiotics'] = pd.get_dummies(abx.antibiotic_state, drop_first=True).to_numpy(float)[idx]
        b['clinical'] = pd.get_dummies(abx.clinical_group, drop_first=True).to_numpy(float)[idx]
    return b


def _estimate(x, y, block, subset):
    """Partial rank correlation after regressing ranks on the subset (build_approved ``estimate``).

    Parameters
    ----------
    x, y : numpy.ndarray
        KO and metabolite values of the analysis subset.
    block : dict
        Output of ``_blocks``.
    subset : tuple of str
        Covariate names in the model (empty for the unadjusted rank correlation).

    Returns
    -------
    float
        Pearson correlation of the rank residuals.
    """
    Z = np.column_stack([np.ones(len(x))] + [block[k] for k in subset])
    ranks = np.column_stack([rankdata(x), rankdata(y)])
    residuals = ranks - Z @ np.linalg.lstsq(Z, ranks, rcond=None)[0]
    return float(np.corrcoef(residuals.T)[0, 1])


def stat_a_adjustment_grid(d, abx):
    """All 16 covariate combinations x 2 normalizations x 2 subsets x 3 pairs (192 estimates).

    Statistical unit: patient. n = 85, 84, 85 (pair-complete) and 84 (common to all
    three pairs). Point estimates only: no P values, no intervals (v28d line 799).
    Estimates that need antibiotic state or clinical group are NaN when the controlled
    input is absent.

    Parameters
    ----------
    d : dict
        Output of ``load``.
    abx : pandas.DataFrame or None
        Output of ``load_controlled``.

    Returns
    -------
    pandas.DataFrame
        Columns as the producer's separate_adjustments_all_models.tsv.
    """
    order = d["patients"]
    feats = list(d["focused_features"].feature)
    ko_total = order.ko_total_signal.to_numpy(float)  # annotated-relative normalization
    subsets = [c for n in range(5) for c in itertools.combinations(COVARIATES, n)]
    masks = {}
    for m, k, *_ in PAIRS:
        x, y = _pair_arrays(d, m, k)
        masks[m] = np.isfinite(_positive(x)) & np.isfinite(_positive(y))
    common = np.logical_and.reduce(list(masks.values()))
    rows = []
    for m, k, *_ in PAIRS:
        for scope, mask in [('pair_complete', masks[m]), ('common_84', common)]:
            idx = np.flatnonzero(mask)
            for norm in NORMALIZATIONS:
                x = d["focused_raw"][feats.index(k)].astype(float)
                if norm == 'annotated_relative':
                    x = x / ko_total
                xx, yy = x[idx], d["focused_raw"][feats.index(m)][idx]
                block = _blocks(order, abx, idx)
                est = [(_estimate(xx, yy, block, s) if all(c in block for c in s) else np.nan) for s in subsets]
                for s, r in zip(subsets, est):
                    rows.append(dict(metabolite=m, KO=k, scope=scope, normalization=norm,
                                     model='+'.join(s) if s else 'none', n=len(idx), rho=r,
                                     delta_from_unadjusted=r - est[0]))
    return pd.DataFrame(rows), int(common.sum())


def stat_a_joint_adjusted(d, abx):
    """Partial rank correlation of each pair adjusted for host share, depth, antibiotics and group.

    Uses the screen covariate matrix (ranks over all 85 patients, dummy-coded antibiotic
    state and clinical group) and the screen ``corr``. Statistical unit: patient; n = 85,
    84, 85. Two-sided t test, df = n - 2 - 6. Needs the controlled input.

    Parameters
    ----------
    d : dict
        Output of ``load``.
    abx : pandas.DataFrame or None
        Output of ``load_controlled``.

    Returns
    -------
    pandas.DataFrame or None
        One row per pair (n, rho, df, P); None without the controlled input.
    """
    if abx is None:
        return None
    order = d["patients"]
    C = np.column_stack([rankdata(order.Human), rankdata(np.log2(order.pg_parsimony)),
                         pd.get_dummies(abx[['antibiotic_state', 'clinical_group']], drop_first=True).to_numpy(float)])
    rows = []
    for m, k, *_ in PAIRS:
        x, y = _pair_arrays(d, m, k)
        n, r, p, df = corr(_positive(x), _positive(y), C)
        rows.append(dict(metabolite=m, KO=k, n=n, rho=r, df=df, P=p))
    return pd.DataFrame(rows)


def stat_family_q(d):
    """Recompute Benjamini-Hochberg q in every saved screen family and read the three pairs.

    Family = one (layer, normalization, model) of the augmented linked-pair screen, all
    tested pairs (finite P). Statistical unit of each test: patient; n per pair as saved
    (43 to 85 patients).

    Parameters
    ----------
    d : dict
        Output of ``load``.

    Returns
    -------
    tuple
        (per-pair table with family size and q for the unadjusted and joint models,
        maximum absolute difference to the saved q_augmented).
    """
    s = d["screen"].copy()
    s["q_recomputed"] = s.groupby(["layer", "normalization", "model"]).p.transform(lambda v: bh(v))
    diff = float(np.nanmax(np.abs(s.q_recomputed - s.q_augmented)))
    rows = []
    for m, k, *_ in PAIRS:
        for model in ["unadjusted", JOINT_MODEL]:
            fam = s[(s.layer == "KO") & (s.normalization == "native") & (s.model == model)]
            r = fam[(fam.metabolite == m) & (fam.protein_feature == k)].iloc[0]
            rows.append(dict(metabolite=m, KO=k, model=model, family_size=int(fam.p.notna().sum()),
                             n=int(r.n), rho=float(r.rho), P=float(r.p), q=float(r.q_recomputed)))
    return pd.DataFrame(rows), diff


def stat_candidates(d):
    """Count native-signal pairs that are inverse, opposite in EDS direction and all q<0.05.

    Criteria (v28d line 799): pair rho < 0 and q < 0.05 in the complete unadjusted native
    KO family; metabolite-EDS rho > 0 and KO-EDS rho < 0, each q < 0.05 in its unadjusted
    native EDS family. Statistical unit of each test: patient; n per test as saved (43 to
    85 patients).

    Parameters
    ----------
    d : dict
        Output of ``load``.

    Returns
    -------
    pandas.DataFrame
        The candidate pairs.
    """
    s, e = d["screen"], d["eds_assoc"]
    s = s[(s.layer == "KO") & (s.normalization == "native") & (s.model == "unadjusted")]
    e = e[(e.normalization == "native") & (e.model == "unadjusted")]
    me = e[e.layer == "Metabolite"].set_index("feature")
    ke = e[e.layer == "KO"].set_index("feature")
    j = s.assign(m_rho=me.rho.reindex(s.metabolite).to_numpy(), m_q=me.q.reindex(s.metabolite).to_numpy(),
                 k_rho=ke.rho.reindex(s.protein_feature).to_numpy(), k_q=ke.q.reindex(s.protein_feature).to_numpy())
    sel = (j.rho < 0) & (j.q_augmented < Q_CUTOFF) & (j.m_rho > 0) & (j.m_q < Q_CUTOFF) & (j.k_rho < 0) & (j.k_q < Q_CUTOFF)
    return j[sel]


def stat_family_sizes(d):
    """Number of tests in each Benjamini-Hochberg family used by the candidate rule.

    Statistical unit of each test: patient.

    Parameters
    ----------
    d : dict
        Output of ``load``.

    Returns
    -------
    dict
        Tested pairs in the native KO unadjusted and joint-adjusted pair families, and
        tested features in the native unadjusted metabolite-EDS and KO-EDS families.
    """
    s, e = d["screen"], d["eds_assoc"]
    out = {}
    for model in ["unadjusted", JOINT_MODEL]:
        f = s[(s.layer == "KO") & (s.normalization == "native") & (s.model == model)]
        out[f"pairs_KO_native_{model}"] = int(f.p.notna().sum())
    for layer in ["Metabolite", "KO"]:
        f = e[(e.layer == layer) & (e.normalization == "native") & (e.model == "unadjusted")]
        out[f"EDS_{layer}_native_unadjusted"] = int(f.p.notna().sum())
    return out


def stat_counts(d):
    """Counts printed in the legend and Methods, and the z-score recomputation check.

    Unit: patient (85) or matched observation (259, original selection).

    Parameters
    ----------
    d : dict
        Output of ``load``.

    Returns
    -------
    dict
        patients, focused rows, overview rows plotted, original features, original
        observations, smallest n among tested screen pairs, and the largest absolute
        difference between stored and recomputed z scores.
    """
    import warnings

    def z(raw):
        """Standardize positive finite values within each feature (row).

        Parameters
        ----------
        raw : numpy.ndarray
            Features x patients raw values.

        Returns
        -------
        numpy.ndarray
            z scores; non-detections stay NaN.
        """
        pos = np.where(np.isfinite(raw) & (raw > 0), raw, np.nan)
        with np.errstate(invalid='ignore', divide='ignore'), warnings.catch_warnings():
            warnings.simplefilter('ignore', RuntimeWarning)  # HV366 has no positive value
            return (pos - np.nanmean(pos, axis=1)[:, None]) / np.nanstd(pos, axis=1)[:, None]
    keep, _ = legacy_rows(d)
    zf, zl = z(d["focused_raw"]), z(d["legacy_raw"])
    dz = max(float(np.nanmax(np.abs(zf - d["focused_z"]))), float(np.nanmax(np.abs(zl[keep] - d["legacy_z"][keep]))))
    same_nan = bool(np.array_equal(np.isnan(zf), np.isnan(d["focused_z"])) and np.array_equal(np.isnan(zl), np.isnan(d["legacy_z"])))
    s = d["screen"]
    return dict(patients=len(d["patients"]), focused_rows=d["focused_z"].shape[0],
                overview_rows=len(keep), original_features=d["legacy_259"].shape[0],
                original_observations=d["legacy_259"].shape[1],
                min_n_tested=int(s.loc[s.p.notna(), "n"].min()), max_n_untested=int(s.loc[s.p.isna(), "n"].max()),
                max_abs_z_difference=dz, z_missing_pattern_identical=same_nan)


def stat_clip(d):
    """Evidence for the display clip at +-2, read from the producer's plotted values.

    The producer wrote every plotted heatmap cell with its z score and the value actually
    drawn (``producer_plotted_ab.tsv``, the heatmap rows of the Source Data plotted-values
    sheet). The clip bound is the largest absolute drawn value; it equals 2 only if the
    producer clipped at 2 and some z score exceeded 2. Unit: heatmap cell (feature x
    patient), 6 x 85 cells in panel a and 19 x 85 in panel b.

    Parameters
    ----------
    d : dict
        Output of ``load``.

    Returns
    -------
    dict
        ``bound`` (max absolute drawn value), ``max_abs_z`` (before clipping),
        ``cells`` (finite drawn cells), ``cells_clipped`` (|z| > bound), ``share_clipped``,
        and ``rebuilt_identical`` (True if ``plotted_heatmap_rows`` rebuilds the producer's
        cells exactly from source_data).
    """
    pp = d["producer_plotted"]
    z, drawn = pp.z.to_numpy(float), pp.plotted.to_numpy(float)
    ok = np.isfinite(drawn)
    bound = float(np.max(np.abs(drawn[ok])))
    clipped = int(np.sum(np.abs(z[ok]) > bound))
    rb = pd.DataFrame(plotted_heatmap_rows(d))
    same = (len(rb) == len(pp)
            and all((rb[c].astype(str).to_numpy() == pp[c].astype(str).to_numpy()).all()
                    for c in ["panel", "row", "label", "patient_position"])
            and all(np.array_equal(rb[c].to_numpy(float), pp[c].to_numpy(float), equal_nan=True)
                    for c in ["z", "plotted"]))
    return dict(bound=bound, max_abs_z=float(np.nanmax(np.abs(z))), cells=int(ok.sum()),
                cells_clipped=clipped, share_clipped=clipped / int(ok.sum()), rebuilt_identical=bool(same))


def numbers(d):
    """Every recomputed number in NUMBERS.tsv, in NUMBERS.tsv row order (OPEN rows excluded).

    The ``v28d_line`` column of this unit holds the paragraph index of the v30 text
    (visible_v30f.txt), see README.md.

    None of these needs the controlled input. Numbers printed as words are returned as
    the word (``WORDS``), capitalised where the manuscript starts a sentence with it.

    Parameters
    ----------
    d : dict
        Output of ``load``.

    Returns
    -------
    list of tuple
        (row number in NUMBERS.tsv, recomputed value, function that computed it).
    """
    c = stat_counts(d)
    u = stat_a_unadjusted(d)
    g, common = stat_a_adjustment_grid(d, None)
    dep = g[(g.scope == 'pair_complete') & (g.normalization == 'native') & (g.model == 'depth')]
    cand = stat_candidates(d)
    clip = stat_clip(d)["bound"]
    max_q = float(np.nanmax(cand[["q_augmented", "m_q", "k_q"]].to_numpy()))
    n_models = g.model.nunique()
    n_pairs = len(PAIRS)
    fam = stat_family_sizes(d)
    w = lambda n, cap=False: WORDS[n].capitalize() if cap else WORDS[n]
    v = [
        # v30 paragraph 73 (Results)
        (c["patients"], "stat_counts"),
        (u.rho[0], "stat_a_unadjusted"), (u.rho[1], "stat_a_unadjusted"), (u.rho[2], "stat_a_unadjusted"),
        (u.n[1], "stat_a_unadjusted"),
        (dep.rho.iloc[0], "stat_a_adjustment_grid"), (dep.rho.iloc[1], "stat_a_adjustment_grid"),
        (dep.rho.iloc[2], "stat_a_adjustment_grid"),
        # v30 paragraph 145 (Methods)
        (c["overview_rows"], "stat_counts"), (c["original_features"], "stat_counts"),
        (c["patients"], "stat_counts"), (w(n_pairs), "PAIRS"),
        (int(min(u.n)), "stat_a_unadjusted"), (int(max(u.n)), "stat_a_unadjusted"),
        (c["min_n_tested"], "stat_counts"), (fam["pairs_KO_native_unadjusted"], "stat_family_sizes"),
        (clip, "stat_clip"),
        # v30 paragraph 239 (legend)
        (w(n_pairs, True), "PAIRS"), (c["patients"], "stat_counts"), (clip, "stat_clip"),
        (c["original_features"], "stat_counts"), (c["original_observations"], "stat_counts"),
        # v30 paragraph 409 (Supplementary Methods)
        (c["original_features"], "stat_counts"), (c["patients"], "stat_counts"),
        (c["min_n_tested"], "stat_counts"),
        (len(cand), "stat_candidates"), (max_q, "stat_candidates"),
        (n_models, "stat_a_adjustment_grid"),
        (u.n[0], "stat_a_unadjusted"), (u.n[1], "stat_a_unadjusted"), (u.n[2], "stat_a_unadjusted"),
        (common, "stat_a_adjustment_grid"), (len(g), "stat_a_adjustment_grid"),
        # v30 paragraph 410 (Supplementary Methods)
        (w(c["focused_rows"]), "stat_counts"), (w(n_pairs), "PAIRS"),
        (c["original_features"], "stat_counts"), (c["original_observations"], "stat_counts"),
    ]
    return [(i + 1, val.item() if hasattr(val, "item") else val, f) for i, (val, f) in enumerate(v)]
