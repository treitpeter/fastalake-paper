"""Supplementary Figure 11 (v30): pathway-level cross-omics profiles and sensitivity analyses.

Panels
------
a, b  KEGG pathway distance-correlation matrices; a PREDICT (77 patients), b CapScan (13
      participants x 5 sampling groups). The native, eggNOG OG and KEGG KO matrices of the
      v28f/v29 panels a, b are removed in v30 (they are the values of Fig. 4i and 5g); their
      tests are still recomputed here and kept in STATS.tsv and Source Data.
c, d  Native-profile correlations, all observed features versus >=10% prevalence filter.
e, f  Controlled native-distance grid M1-M8; f CapScan region adjusted.

``render`` ports the v30 build script (build_sfig11.py): one 180 x 190 mm canvas in the
repository fl_style at 6 pt, the drawing code and row selection unchanged, paths replaced.

The statistics functions recompute every plotted coefficient, permutation P value and
Benjamini-Hochberg q value from the saved distance matrices, using the upstream code
(``crossomics.py`` -> ``check_prevalence.evaluate`` for functional profiles,
``native_pyramids.correlations`` for native profiles, ``experiment.correlate`` for the
controlled grid) with the arithmetic unchanged.

Source data
-----------
source_data/SFig11/functional_correlations_PREDICT.tsv     PREDICT rows of the producer input
source_data/Fig5/functional_correlations.tsv               CapScan rows of the producer input
source_data/SFig11/native_correlations_PREDICT.tsv         PREDICT rows of the producer input
source_data/Fig5/native_correlations.tsv                   CapScan rows of the producer input
source_data/SFig11/controlled_176_comparisons.tsv          producer input, panels e, f
source_data/SFig11/distances_PREDICT_{functional,native,controlled}.tsv
source_data/SFig11/distances_CAPSCAN_native_prevalence10.tsv
source_data/SFig11/distances_CAPSCAN_controlled.tsv
source_data/Fig5/functional_distances_long.tsv             CapScan functional distances
source_data/Fig5/native_distances_long.tsv                 CapScan native distances, all observed
source_data/SFig11/predict_native_retained_signal.tsv      Supplementary Methods (v30 paragraph 416)
"""

import hashlib
import io
import os
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

# ---- choices (each a value the original code used) ---------------------------------------
SEED = 20260907                 # permutation seed; Supplementary Methods 8, v28d line 815
N_PERMUTATIONS = 1999           # v28d lines 453, 815, 821, 827, 835
P_TOLERANCE_NATIVE = 1e-12      # |null| >= |rho| - 1e-12 in native and controlled tests (code only)
FUNCTIONAL_CHUNK = 200          # permutation chunk in check_prevalence.evaluate (no effect on P)
Q_CUTOFF = 0.05                 # asterisks; v28d line 453
REGIONS = ["duodenum", "jejunum", "ileum", "colon", "stool"]   # v28d line 161 (five sampling groups)
CAPSCAN_DONORS = [1, 2, 3, 4, 6, 8, 9, 10, 11, 12, 13, 14, 15]  # 13 DNA-complete participants; v28d lines 453, 835
FUNCTIONAL_PREVALENCE = 0.0     # all observed features in a, b; v28d lines 279, 813
PREDICT_SCOPE = "displayed77"   # 77 PREDICT patients; v28d line 453
CAPSCAN_SCOPE = "matched12"     # producer label for the DNA-matched family (13 participants at 347)
FIG_MM = (180.0, 190.0)         # v30 build: authored at the placed size
FONT_PT = 6.0                   # v30 build: 6 pt text
MATRIX_CELL_MM = 7.6            # v30 build: matrix cell size in a, b
PNG_DPI = 300                   # v30 build: PNG resolution

PREDICT_FUNCTIONAL_LAYERS = ["metaG", "metaT", "microbial_metaP", "host_metaP", "metabolites"]
CAPSCAN_FUNCTIONAL_LAYERS = ["metaG", "microbial_metaP", "host_metaP", "metabolites"]
PREDICT_NATIVE_LAYERS = ["DNA", "RNA", "Non-host protein", "Host protein", "Metabolites"]
CAPSCAN_NATIVE_LAYERS = ["DNA", "Non-host protein", "Host protein", "Metabolites"]
NUMBER_WORDS = {1: "one", 2: "two", 3: "three", 4: "four", 5: "five", 6: "six", 7: "seven",
                8: "eight", 9: "nine", 10: "ten", 11: "eleven", 12: "twelve", 13: "thirteen",
                14: "fourteen", 15: "fifteen", 16: "sixteen"}

#: sha256 of the shared tables written by the Fig5 unit that this figure reads.
FIG5_SHA256 = {
    "functional_correlations.tsv": "455bfa9d3185aadaba071d640f77241fe9e11c898768ab8c6b86513b0cb1cbf4",
    "native_correlations.tsv": "0004f35d522422580e1d7cce3bbcbc2d2967a8fbc8af28668b73408d52305375",
    "functional_distances_long.tsv": "b0bd1978b64a802fb1874ec0a5e1c684fc286550b0ce172ef83443707d2b4fa4",
    "native_distances_long.tsv": "07bda50dc99acf5fc5e6a77ec91d9fd71bc90451a924c2413cac193de6de4530",
}


def sha256(path):
    """Return the hex sha256 of a file."""
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def rd(p, **kw):
    """Read a tab-separated table exactly as the producer did (``common.rd``)."""
    return pd.read_csv(p, sep="\t", **kw)


# ---- loading ------------------------------------------------------------------------------

def load_inputs(sd, fig5):
    """Load the three producer input tables.

    Parameters
    ----------
    sd, fig5 : pathlib.Path
        ``source_data/SFig11`` and ``source_data/Fig5``.

    Returns
    -------
    dict
        ``functional``, ``native`` (PREDICT rows then CapScan rows, as in the producer
        inputs) and ``grid`` (176 controlled estimates), as pandas DataFrames.
    """
    for name, h in FIG5_SHA256.items():
        assert sha256(fig5 / name) == h, f"source_data/Fig5/{name} changed"

    def joined(name):
        # Header + PREDICT lines + CapScan lines is byte-identical to the producer input.
        a = (sd / name.replace(".tsv", "_PREDICT.tsv")).read_text(encoding="utf-8").splitlines(keepends=True)
        b = (fig5 / name).read_text(encoding="utf-8").splitlines(keepends=True)
        assert a[0] == b[0]
        return rd(io.StringIO("".join(a + b[1:])))

    return dict(functional=joined("functional_correlations.tsv"), native=joined("native_correlations.tsv"),
                grid=rd(sd / "controlled_176_comparisons.tsv"))


def load_distances(path, keycols, capscan, filters=None):
    """Rebuild square distance matrices from a long upper-triangle table.

    Parameters
    ----------
    path : pathlib.Path
        Long table with key columns, profile columns and ``distance``.
    keycols : list of str
        Columns that identify one matrix.
    capscan : bool
        True if profiles are (participant, region) pairs in ``CAPSCAN_DONORS`` x ``REGIONS``
        order; False for PREDICT positions ``profile_i``/``profile_j``.
    filters : dict or None
        Column -> value restrictions applied before grouping.

    Returns
    -------
    dict
        {key tuple: (n, n) symmetric matrix with zero diagonal}. The upper triangle in
        row-major order is the stored vector, exactly.
    """
    from scipy.spatial.distance import squareform

    t = pd.read_csv(path, sep="\t", float_precision="round_trip")
    for c, v in (filters or {}).items():
        t = t[t[c] == v]
    out = {}
    for key, g in t.groupby(keycols, sort=False):
        key = key if isinstance(key, tuple) else (key,)
        if capscan:
            prof = [(d, r) for d in CAPSCAN_DONORS for r in REGIONS]
            n = len(prof)
            iu = np.triu_indices(n, 1)
            assert list(zip(g.participant_i, g.region_i)) == [prof[i] for i in iu[0]]
            assert list(zip(g.participant_j, g.region_j)) == [prof[j] for j in iu[1]]
        else:
            n = int(g.profile_j.max())
            iu = np.triu_indices(n, 1)
            assert np.array_equal(g.profile_i.to_numpy() - 1, iu[0]) and np.array_equal(g.profile_j.to_numpy() - 1, iu[1])
        out[key] = squareform(g.distance.to_numpy(), checks=False)
    return out


# ---- v30 figure ---------------------------------------------------------------------------

SHORT = {'DNA': 'DNA', 'RNA': 'RNA', 'Non-host protein': 'Protein', 'Host protein': 'Human', 'Metabolites': 'Metabolites'}


def render(inputs, out_stem):
    """Draw Supplementary Fig. 11 (v30) on one 180 x 190 mm canvas, as build_sfig11.py did.

    Ported from the v30 build script with paths replaced; drawing and value selection are
    unchanged. a, b: KEGG pathway distance-correlation matrices (PREDICT, CapScan; all
    observed features). c, d: native-profile correlations, all observed versus >=10%
    prevalence. e, f: controlled native-distance grid M1-M8 (CapScan region adjusted). The
    v30 build's duplicate test against the v29 build tables (native, OG and KO matrices
    equal Fig. 4i and 5g) needs files outside this repository and is not repeated here.

    Parameters
    ----------
    inputs : dict
        Output of ``load_inputs`` (``functional``, ``native``, ``grid``).
    out_stem : pathlib.Path
        Output path without extension; ``.pdf`` and ``.png`` (``PNG_DPI``) are written.

    Returns
    -------
    pandas.DataFrame
        Plotted values, one row per drawn number, as the v30 build wrote plotted_values.tsv.
    """
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap
    from matplotlib.patches import Rectangle
    from style import fl_style as fl

    functional, native, grid = inputs["functional"], inputs["native"], inputs["grid"]
    FS = FONT_PT
    fl.init()
    cm = LinearSegmentedColormap.from_list('fl', [fl.BLUE, 'white', fl.ORANGE])
    W, H = FIG_MM
    fig = plt.figure(figsize=fl.mm(W, H))

    def ax_mm(x, y, w, h):  # x,y = left, top in mm
        return fig.add_axes([x / W, 1 - (y + h) / H, w / W, h / H])

    def letter(x, y, s):
        fig.text(x / W, 1 - y / H, s, fontsize=7, fontweight='bold', color=fl.INK, va='top', ha='left')

    # ---- a,b: KEGG pathway matrices ----
    records = []

    def pyramid(ax, co):
        layers = ['DNA', 'RNA', 'Protein', 'Human protein', 'Metabolites'] if co == 'PREDICT' else ['DNA', 'Protein', 'Human protein', 'Metabolites']
        km = {'DNA': 'metaG', 'RNA': 'metaT', 'Protein': 'microbial_metaP', 'Human protein': 'host_metaP', 'Metabolites': 'metabolites'}
        t = functional[(functional.cohort == co) & (functional.level == 'Pathway') & (functional.prevalence == 0) & (functional.scope == ('displayed77' if co == 'PREDICT' else 'matched12'))]
        lookup = {frozenset((r.layer_a, r.layer_b)): r for r in t.itertuples()}
        n = len(layers)
        for i, la in enumerate(layers):
            for j, lb in enumerate(layers[:i + 1]):
                r = lookup.get(frozenset((km[la], km[lb])))
                if i == j: rho, lab = 1., '1.00'
                elif r is None: rho, lab = 0., '—'; assert {la, lb} == {'Human protein', 'Metabolites'}
                else:
                    rho, lab = r.rho, f'{r.rho:.2f}' + ('*' if r.q < .05 else '')
                    records.append(dict(panel='a' if co == 'PREDICT' else 'b', cohort=co, level='Pathway', layer_a=la, layer_b=lb, rho=rho, q=r.q, label=lab))
                ax.add_patch(Rectangle((j, i), 1, 1, facecolor=cm((rho + 1) / 2) if lab != '—' else '#f5f5f5', edgecolor='white', lw=1))
                ax.text(j + .5, i + .5, lab, ha='center', va='center', fontsize=FS)
        ax.set(xlim=(0, n), ylim=(n, 0), xticks=np.arange(n) + .5, yticks=np.arange(n) + .5)
        ax.set_xticklabels([str(k + 1) for k in range(n)], fontsize=FS); ax.set_yticklabels([f'{k + 1} {l}' for k, l in enumerate(layers)], fontsize=FS)
        ax.tick_params(length=0, pad=2)
        for s in ax.spines.values(): s.set_visible(False)
    C = MATRIX_CELL_MM
    pyramid(ax_mm(26, 6, 5 * C, 5 * C), 'PREDICT'); letter(2, 4, 'a')
    pyramid(ax_mm(118, 6, 4 * C, 4 * C), 'CAPSCAN'); letter(94, 4, 'b')

    # ---- c,d: native prevalence ----
    short = SHORT
    for (x, w, l, lx), co in zip([(26, 58, 'c', 2), (118, 56, 'd', 94)], ['PREDICT', 'CAPSCAN']):
        ax = ax_mm(x, 56, w, 46); letter(lx, 54, l)
        t = native[native.cohort.eq(co)]; keys = list(dict.fromkeys(zip(t.layer_a, t.layer_b))); y = np.arange(len(keys))
        for analysis, color, off, lab in [('sensitivity_all_features', fl.BLUE, -.13, 'All observed'), ('primary_prevalence10', fl.ORANGE, .13, '10% prevalence')]:
            s = t[t.analysis.eq(analysis)].set_index(['layer_a', 'layer_b']).loc[keys]
            ax.scatter(s.rho, y + off, s=9, c=color, lw=0, label=lab, zorder=3)
            records.extend(dict(panel=l, cohort=co, analysis=analysis, layer_a=a, layer_b=b, rho=v) for (a, b), v in s.rho.items())
        ax.set(yticks=y, xlim=(-.1, .85), xlabel='Spearman distance correlation')
        ax.set_yticklabels([short[a] + '–' + short[b] for a, b in keys]); ax.invert_yaxis(); ax.axvline(0, c='#bbbbbb', lw=1, zorder=1)
        fl.style(ax); ax.tick_params(labelsize=FS); ax.xaxis.label.set_size(FS)
        if l == 'c': handles = ax.get_legend_handles_labels()
    fig.legend(*handles, loc='upper center', bbox_to_anchor=(.5, 1 - 112 / H), ncol=2, frameon=False, fontsize=FS, handletextpad=.3, columnspacing=1.5)

    # ---- e,f: controlled native-distance methods ----
    methods = [(v, d, o) for v in ['prevalence10', 'all_observed'] for d in ['hellinger', 'rank'] for o in ['pearson', 'spearman']]
    for (x, w, l, lx), co in zip([(26, 60, 'e', 2), (112, 56, 'f', 90)], ['PREDICT', 'CAPSCAN']):
        t = grid[(grid.cohort == co) & (grid.adjustment == ('unadjusted' if co == 'PREDICT' else 'region_adjusted'))]
        pairs = list(dict.fromkeys(zip(t.layer_a, t.layer_b))); v = np.zeros((len(pairs), 8)); q = np.zeros_like(v)
        for i, (a, b) in enumerate(pairs):
            for j, (f, d, o) in enumerate(methods):
                row = t[(t.layer_a == a) & (t.layer_b == b) & (t.features == f) & (t.distance == d) & (t.outer_correlation == o)].iloc[0]
                v[i, j], q[i, j] = row.coefficient, row.q
                records.append(dict(panel=l, cohort=co, layer_a=a, layer_b=b, features=f, distance=d, outer_correlation=o, coefficient=row.coefficient, q=row.q))
        hh = 5.4 * len(pairs)
        ax = ax_mm(x, 126, w, hh); letter(lx, 124, l)
        im = ax.imshow(v, vmin=-1, vmax=1, cmap=cm, aspect='auto', interpolation='nearest')
        ax.set(yticks=range(len(pairs)), xticks=range(8)); ax.set_yticklabels([short[a] + '–' + short[b] for a, b in pairs], fontsize=FS); ax.set_xticklabels(['M' + str(k) for k in range(1, 9)], fontsize=FS)
        ax.tick_params(length=0, pad=2)
        for s in ax.spines.values(): s.set_visible(False)
        for i in range(len(pairs)):
            for j in range(8): ax.text(j, i, f'{v[i, j]:.2f}' + ('*' if q[i, j] < .05 else ''), ha='center', va='center', fontsize=FS)
    cax = ax_mm(171, 126, 2.2, 32); cb = fig.colorbar(im, cax=cax, ticks=[-1, 0, 1]); cb.ax.tick_params(labelsize=FS, length=2, pad=1)
    cb.set_label('Distance correlation', fontsize=FS, labelpad=2)

    fl.save(fig, out_stem, dpi=PNG_DPI, formats=('pdf', 'png'))
    plt.close(fig)
    return pd.DataFrame(records)


# ---- statistics ---------------------------------------------------------------------------

def _permutations(n, block):
    """Permutations exactly as the upstream code draws them (one fresh generator per call)."""
    rng = np.random.default_rng(SEED)
    if block:
        blocks = np.arange(n).reshape(-1, 5)
        return np.array([blocks[rng.permutation(len(blocks))].ravel() for _ in range(N_PERMUTATIONS)])
    return np.array([rng.permutation(n) for _ in range(N_PERMUTATIONS)])


def bh(values):
    """Benjamini-Hochberg q values (moved from ``native_pyramids.bh``)."""
    v=np.asarray(values);o=np.argsort(v);q=np.empty(len(v));q[o]=np.minimum(1,np.minimum.accumulate((v[o]*len(v)/np.arange(1,len(v)+1))[::-1])[::-1]);return q


def stat_functional_family(mats, cohort, layers):
    """Recompute one functional-profile family (panels a, b; v28d lines 139, 835).

    Port of ``check_prevalence.evaluate`` as run by ``crossomics.py`` (NPERM 1,999, seed
    20260907): Spearman correlation between upper triangles of two Euclidean distance
    matrices of Hellinger-transformed (or standardized log10 metabolite) profiles,
    two-sided permutation P = (1 + #{|null| >= |rho|}) / (1 + 1,999), rows and columns of
    one matrix permuted together (PREDICT: patient labels; CapScan: whole five-group
    participant blocks). Host protein-metabolite is tested only at OG. Benjamini-Hochberg
    (statsmodels ``fdr_bh``) over the whole family (all levels).

    Statistical unit: patient (PREDICT, n = 77) or participant (CapScan, n = 13, 65
    participant-by-group profiles).

    Parameters
    ----------
    mats : dict
        {(level, layer): square distance matrix}.
    cohort : str
        ``"PREDICT"`` or ``"CAPSCAN"``.
    layers : list of str
        Layer order of the upstream frames dict (defines pair order).

    Returns
    -------
    pandas.DataFrame
        cohort, level, layer_a, layer_b, n_profiles, rho, p, q.
    """
    from scipy.stats import rankdata
    from scipy.spatial.distance import squareform
    from statsmodels.stats.multitest import multipletests

    rows = []
    for lev in dict.fromkeys(k[0] for k in mats):
        n = len(mats[(lev, layers[0])]); ix = np.triu_indices(n, 1)
        perms = _permutations(n, cohort == "CAPSCAN")
        ranks = {}
        for k in layers:
            r = rankdata(mats[(lev, k)][ix]); r -= r.mean(); r /= np.linalg.norm(r); ranks[k] = r
        for a, b in combinations(layers, 2):
            if lev != 'OG' and {a, b} == {'host_metaP', 'metabolites'}:
                continue
            rho = float(ranks[a] @ ranks[b])
            bm = squareform(ranks[b]); ext = 0
            for start in range(0, N_PERMUTATIONS, FUNCTIONAL_CHUNK):
                pp = perms[start:start + FUNCTIONAL_CHUNK]; null = bm[pp[:, ix[0]], pp[:, ix[1]]] @ ranks[a]
                ext += np.count_nonzero(np.abs(null) >= abs(rho))
            rows.append(dict(cohort=cohort, level=lev, layer_a=a, layer_b=b, n_profiles=n, rho=rho, p=(1 + ext) / (1 + N_PERMUTATIONS)))
    t = pd.DataFrame(rows)
    t["q"] = multipletests(t.p, method="fdr_bh")[1]
    return t


def stat_native(mats, cohort):
    """Recompute one native-profile family (panels a-d; legend v28d line 453).

    Port of ``native_pyramids.correlations``: Spearman correlation of distance-matrix upper
    triangles; two-sided permutation P = (1 + #{|null| >= |rho| - 1e-12}) / 2,000 with
    1,999 permutations (seed 20260907; PREDICT patient labels, CapScan whole participant
    blocks); Benjamini-Hochberg within the family (10 PREDICT or 6 CapScan comparisons).

    Statistical unit: patient (PREDICT, n = 77) or participant (CapScan, n = 13).

    Parameters
    ----------
    mats : dict
        {layer: square distance matrix}, in upstream layer order.
    cohort : str
        ``"PREDICT"`` or ``"CAPSCAN"``.

    Returns
    -------
    pandas.DataFrame
        cohort, layer_a, layer_b, rho, p, n_profiles, n_participants, q.
    """
    from scipy.stats import rankdata, spearmanr
    from scipy.spatial.distance import squareform

    distances = mats
    n=len(next(iter(distances.values())));ix=np.triu_indices(n,1)
    perms = _permutations(n, cohort == "CAPSCAN")
    ranks={k:rankdata(v[ix]) for k,v in distances.items()};norm={}
    for k,r in ranks.items():
        centered=r-r.mean();scale=np.linalg.norm(centered);assert scale>1e-10
        norm[k]=centered/scale
    rows=[]
    for a,b in combinations(distances,2):
        rho=float(norm[a]@norm[b]);checkrho=float(spearmanr(distances[a][ix],distances[b][ix]).statistic)
        assert abs(rho-checkrho)<1e-12
        bmat=squareform(norm[b]);extreme=0
        for chunk in np.array_split(perms,20):
            null=bmat[chunk[:,ix[0]],chunk[:,ix[1]]]@norm[a]
            extreme+=int((abs(null)>=abs(rho)-P_TOLERANCE_NATIVE).sum())
        rows.append(dict(cohort=cohort,layer_a=a,layer_b=b,rho=rho,p=(1+extreme)/(1+N_PERMUTATIONS),n_profiles=n,n_participants=n//5 if cohort=='CAPSCAN' else n))
    table=pd.DataFrame(rows);table['q']=bh(table.p);return table


def stat_controlled(mats, cohort, outer):
    """Recompute one controlled-grid family (panels e, f; v28d lines 161, 453, 823).

    Port of ``experiment.correlate``: outer Pearson or Spearman correlation between the
    upper triangles of two distance matrices; two-sided P = (1 + #{|null| >= |r| - 1e-12})
    / 2,000 and positive-tail P = (1 + #{null >= r - 1e-12}) / 2,000 with 1,999
    permutations (seed 20260907, same permutations as the native test); Benjamini-Hochberg
    on the two-sided P within the family (10 PREDICT or 6 CapScan comparisons).

    Statistical unit: patient (PREDICT, n = 77) or participant (CapScan, n = 13).

    Parameters
    ----------
    mats : dict
        {layer: square distance matrix}, in upstream layer order.
    cohort : str
        ``"PREDICT"`` or ``"CAPSCAN"``.
    outer : str
        ``"pearson"`` or ``"spearman"``.

    Returns
    -------
    pandas.DataFrame
        cohort, layer_a, layer_b, coefficient, p, p_positive_tail, n_profiles, n_participants, q.
    """
    from scipy.stats import pearsonr, rankdata, spearmanr
    from scipy.spatial.distance import squareform

    ds = mats; co = cohort
    n=len(next(iter(ds.values())));ix=np.triu_indices(n,1)
    perms = _permutations(n, co == "CAPSCAN")
    norm={}
    for key,d in ds.items():
        v=rankdata(d[ix]) if outer=='spearman' else d[ix].copy();v-=v.mean();norm[key]=v/np.linalg.norm(v)
    rows=[]
    for a,b in combinations(ds,2):
        r=float(norm[a]@norm[b]);ind=(spearmanr if outer=='spearman' else pearsonr)(ds[a][ix],ds[b][ix]).statistic;assert abs(ind-r)<1e-12
        mat=squareform(norm[b]);null=mat[perms[:,ix[0]],perms[:,ix[1]]]@norm[a]
        rows.append(dict(cohort=co,layer_a=a,layer_b=b,coefficient=r,p=(1+int((abs(null)>=abs(r)-P_TOLERANCE_NATIVE).sum()))/(1+N_PERMUTATIONS),
            p_positive_tail=(1+int((null>=r-P_TOLERANCE_NATIVE).sum()))/(1+N_PERMUTATIONS),n_profiles=n,n_participants=n//5 if co=='CAPSCAN' else n))
    t=pd.DataFrame(rows);t['q']=bh(t.p);return t


def recompute_all(sd, fig5):
    """Recompute every plotted test from the distance matrices.

    Returns
    -------
    dict
        ``functional`` (PREDICT 28 + CapScan 16 rows), ``native`` (32 rows) and
        ``controlled`` (176 rows) DataFrames with the key columns of the producer tables.
    """
    pf = load_distances(sd / "distances_PREDICT_functional.tsv", ["level", "layer"], False)
    cf = load_distances(fig5 / "functional_distances_long.tsv", ["level", "layer"], True,
                        {"scope": CAPSCAN_SCOPE, "prevalence": FUNCTIONAL_PREVALENCE})
    f = pd.concat([stat_functional_family(pf, "PREDICT", PREDICT_FUNCTIONAL_LAYERS).assign(scope=PREDICT_SCOPE),
                   stat_functional_family(cf, "CAPSCAN", CAPSCAN_FUNCTIONAL_LAYERS).assign(scope=CAPSCAN_SCOPE)],
                  ignore_index=True).assign(prevalence=FUNCTIONAL_PREVALENCE)

    pn = load_distances(sd / "distances_PREDICT_native.tsv", ["analysis", "layer"], False)
    cn = load_distances(sd / "distances_CAPSCAN_native_prevalence10.tsv", ["analysis", "layer"], True)
    cn.update(load_distances(fig5 / "native_distances_long.tsv", ["analysis", "layer"], True))
    nat = []
    for co, mats, layers in [("PREDICT", pn, PREDICT_NATIVE_LAYERS), ("CAPSCAN", cn, CAPSCAN_NATIVE_LAYERS)]:
        for an in ["primary_prevalence10", "sensitivity_all_features"]:
            t = stat_native({k: mats[(an, k)] for k in layers}, co); t.insert(1, "analysis", an); nat.append(t)
    nat = pd.concat(nat, ignore_index=True)

    pc = load_distances(sd / "distances_PREDICT_controlled.tsv", ["features", "distance_kind", "adjustment", "layer"], False)
    cc = load_distances(sd / "distances_CAPSCAN_controlled.tsv", ["features", "distance_kind", "adjustment", "layer"], True)
    rows = []
    for co, mats, layers, adjs in [("PREDICT", pc, PREDICT_NATIVE_LAYERS, ["unadjusted"]),
                                   ("CAPSCAN", cc, CAPSCAN_NATIVE_LAYERS, ["unadjusted", "region_adjusted"])]:
        for feat in ["all_observed", "prevalence10"]:
            for metric in ["hellinger", "rank"]:
                for adj in adjs:
                    for outer in ["pearson", "spearman"]:
                        t = stat_controlled({k: mats[(feat, metric, adj, k)] for k in layers}, co, outer)
                        t.insert(1, "features", feat); t.insert(2, "distance", metric); t.insert(3, "adjustment", adj); t.insert(4, "outer_correlation", outer)
                        rows.append(t)
    return dict(functional=f, native=nat, controlled=pd.concat(rows, ignore_index=True))


def compare_to_inputs(rec, inputs):
    """Maximum absolute difference between recomputed and producer-table values.

    Returns
    -------
    pandas.DataFrame
        One row per (table, column): rows compared, max |difference|, and the number of
        q < 0.05 classifications that differ.
    """
    out = []
    specs = [("functional", ["cohort", "scope", "level", "layer_a", "layer_b", "prevalence"], ["rho", "p", "q"]),
             ("native", ["cohort", "analysis", "layer_a", "layer_b"], ["rho", "p", "q"]),
             ("controlled", ["cohort", "features", "distance", "adjustment", "outer_correlation", "layer_a", "layer_b"],
              ["coefficient", "p", "p_positive_tail", "q"])]
    for name, keys, cols in specs:
        src = inputs["grid"] if name == "controlled" else inputs[name]
        m = rec[name].merge(src, on=keys, suffixes=("_recomputed", "_table"), validate="one_to_one")
        assert len(m) == len(rec[name])
        for c in cols:
            out.append(dict(table=name, column=c, rows=len(m),
                            max_abs_difference=float(np.max(np.abs(m[c + "_recomputed"] - m[c + "_table"]))),
                            q_class_changes=int(((m["q_recomputed"] < Q_CUTOFF) != (m["q_table"] < Q_CUTOFF)).sum())))
    return pd.DataFrame(out)


def stat_line139(functional):
    """PREDICT pathway RNA-protein agreement (v30 paragraph 74): rho and BH q of the 28-test family.

    Statistical unit: patient, n = 77.

    Returns
    -------
    dict
        rho, p, q.
    """
    r = functional[(functional.cohort == "PREDICT") & (functional.scope == PREDICT_SCOPE) & (functional.prevalence == 0)
                   & (functional.level == "Pathway") & (functional.layer_a == "metaT") & (functional.layer_b == "microbial_metaP")]
    assert len(r) == 1
    r = r.iloc[0]
    return dict(rho=float(r.rho), p=float(r.p), q=float(r.q))


def stat_panel_f_counts(controlled):
    """Counts behind v30 paragraphs 83 and 421: region-adjusted CapScan controlled comparisons.

    Statistical unit: participant, n = 13.

    Returns
    -------
    dict
        Number of region-adjusted rank and Hellinger comparisons, rank comparisons with
        q >= 0.05 (overall, by layer pair, and involving metabolites) and Hellinger
        comparisons with q < 0.05.
    """
    t = controlled[(controlled.cohort == "CAPSCAN") & (controlled.adjustment == "region_adjusted")]
    rk, hl = t[t.distance == "rank"], t[t.distance == "hellinger"]
    ns = rk[rk.q >= Q_CUTOFF]
    by = {f"{a}–{b}": int(((ns.layer_a == a) & (ns.layer_b == b)).sum()) for a, b in dict.fromkeys(zip(rk.layer_a, rk.layer_b))}
    return dict(rank_n=len(rk), rank_not_significant=len(ns), rank_not_significant_by_pair=by,
                rank_not_significant_with_metabolites=int(((ns.layer_a == "Metabolites") | (ns.layer_b == "Metabolites")).sum()),
                hellinger_n=len(hl), hellinger_significant=int((hl.q < Q_CUTOFF).sum()))


def stat_retention(t):
    """Median retained sequence-resolved abundance under the >=10% prevalence filter (v30 paragraph 416).

    Statistical unit: patient, n = 77 per layer.

    Returns
    -------
    dict
        {layer: median retained fraction} for DNA and RNA.
    """
    t = t[(t.cohort == "PREDICT") & (t.analysis == "primary_prevalence10")]
    return {k: float(np.median(g.retained_signal_fraction.to_numpy())) for k, g in t.groupby("layer")}


def stat_design(inputs, rec):
    """Population sizes and family sizes printed in the legend and Methods (v30 paragraphs 241, 427).

    Returns
    -------
    dict
        PREDICT patients, CapScan participants and profiles per participant, permutations
        implied by the smallest P, controlled estimates, CapScan functional family size.
    """
    nat, grid, fun = inputs["native"], inputs["grid"], inputs["functional"]
    pred = set(nat.loc[nat.cohort == "PREDICT", "n_participants"]) | set(grid.loc[grid.cohort == "PREDICT", "n_participants"]) \
        | set(fun.loc[(fun.cohort == "PREDICT") & (fun.scope == PREDICT_SCOPE), "n_donors"])
    cap = set(nat.loc[nat.cohort == "CAPSCAN", "n_participants"]) | set(grid.loc[grid.cohort == "CAPSCAN", "n_participants"]) \
        | set(fun.loc[(fun.cohort == "CAPSCAN") & (fun.scope == CAPSCAN_SCOPE), "n_donors"])
    prof = set(nat.loc[nat.cohort == "CAPSCAN", "n_profiles"]) | set(grid.loc[grid.cohort == "CAPSCAN", "n_profiles"])
    assert len(pred) == len(cap) == len(prof) == 1
    minp = min(rec["functional"].p.min(), rec["native"].p.min(), rec["controlled"].p.min())
    fam = rec["functional"]
    return dict(predict_patients=pred.pop(), capscan_participants=cap.pop(),
                capscan_groups_per_participant=prof.pop() // len(CAPSCAN_DONORS),
                permutations_from_min_p=round(1 / minp) - 1, controlled_estimates=len(grid),
                capscan_functional_family=int((fam.cohort == "CAPSCAN").sum()),
                predict_functional_family=int((fam.cohort == "PREDICT").sum()))


def plotted_values(records_pyramid, records_supp):
    """Plotted-value tables in the producer's format (``all_plotted_crossomics_values.tsv`` rows of
    this figure and ``supplementary_plotted_values.tsv``)."""
    return pd.DataFrame(records_pyramid), pd.DataFrame(records_supp)


def stats_table(rec, versions):
    """Build STATS.tsv rows (one per test) from the recomputed results.

    Parameters
    ----------
    rec : dict
        Output of ``recompute_all``.
    versions : str
        Software versions, e.g. ``"numpy 2.2.6; scipy 1.16.3; statsmodels 0.14.5"``.

    Returns
    -------
    pandas.DataFrame
        Columns as CONVENTIONS.md section 4.5.
    """
    rows = []
    unit = {"PREDICT": "patient", "CAPSCAN": "participant"}
    n = {"PREDICT": "77", "CAPSCAN": "13 (65 participant-by-group profiles)"}
    rep = {"PREDICT": "distinct", "CAPSCAN": "repeated: five sampling groups within participant"}
    perm = {"PREDICT": "patient labels", "CAPSCAN": "whole five-group participant blocks"}

    def add(claim, line, test, co, cov, family, fsize, est_name, est, p, q, code):
        rows.append(dict(claim=claim, v28d_line=line, test=test, sided="two-sided", unit=unit[co], n=n[co],
                         distinct_or_repeated=rep[co], covariates=cov, correction="Benjamini-Hochberg",
                         family=family, family_size=fsize, estimate=f"{est_name} {est!r}; BH q {q!r}",
                         uncertainty="not computed", uncertainty_construction="NA",
                         statistic=f"{est_name} = {est!r}", df="NA", exact_P=repr(p), effect_size=repr(est),
                         effect_size_how=f"{est_name} between the upper triangles of the two distance matrices",
                         draws=str(N_PERMUTATIONS), seed=str(SEED), software_version=versions, code=code))

    f = rec["functional"]
    for r in f.itertuples():
        co = r.cohort
        line = "74" if (co == "PREDICT" and r.level == "Pathway" and r.layer_a == "metaT" and r.layer_b == "microbial_metaP") else ("241" if co == "PREDICT" else "241; 427")
        cov = "none" if co == "PREDICT" else "sampling group (intercept and four indicators), residualized before distances"
        fam = (f"{co} functional profiles, all observed features, OG/KO/pathway, "
               + ("77 patients" if co == "PREDICT" else "13 DNA-complete participants"))
        where = (f"SF11{'a' if co == 'PREDICT' else 'b'}" if r.level == "Pathway"
                 else f"Source Data (shown in Fig. {'4i' if co == 'PREDICT' else '5g'})")
        add(f"{where} {r.level} {r.layer_a} vs {r.layer_b}", line,
            f"Spearman correlation of Euclidean distance matrices (Mantel-type), permutation of {perm[co]}",
            co, cov, fam, int((f.cohort == co).sum()), "Spearman rho", r.rho, r.p, r.q,
            "figures/SFig11/panels.py:stat_functional_family")
    nat = rec["native"]
    for r in nat.itertuples():
        co = r.cohort
        cov = "none" if co == "PREDICT" else "sampling group (residualized before distances)"
        lab = "c" if co == "PREDICT" else "d"
        add(f"SF11{lab} native {r.analysis} {r.layer_a} vs {r.layer_b}", "241",
            f"Spearman correlation of Euclidean distance matrices (Mantel-type), permutation of {perm[co]}",
            co, cov, f"{co} native profiles, {r.analysis}", int(((nat.cohort == co) & (nat.analysis == r.analysis)).sum()),
            "Spearman rho", r.rho, r.p, r.q, "figures/SFig11/panels.py:stat_native")
    c = rec["controlled"]
    for r in c.itertuples():
        co = r.cohort
        cov = "sampling group (region means removed)" if r.adjustment == "region_adjusted" else "none"
        name = "Spearman rho" if r.outer_correlation == "spearman" else "Pearson r"
        line = "241; 83; 421" if (co == "CAPSCAN" and r.adjustment == "region_adjusted") else "241"
        fam = f"{co} controlled grid, {r.features}/{r.distance}/{r.adjustment}/{r.outer_correlation}"
        size = int(((c.cohort == co) & (c.features == r.features) & (c.distance == r.distance)
                    & (c.adjustment == r.adjustment) & (c.outer_correlation == r.outer_correlation)).sum())
        add(f"SF11{'e' if co == 'PREDICT' else 'f'} {r.features} {r.distance} {r.adjustment} {r.outer_correlation} {r.layer_a} vs {r.layer_b}",
            line, f"{name.split()[0]} correlation of {r.distance} distance matrices, permutation of {perm[co]}",
            co, cov, fam, size, name, r.coefficient, r.p, r.q, "figures/SFig11/panels.py:stat_controlled")
    return pd.DataFrame(rows)
