"""Figure 5: personal functional profiles across intestinal sampling regions (CapScan, n=15).

Panels
------
a  study design (artwork, ``artwork/panel_a.pdf``; typed numbers checked against
   ``design_347.tsv`` by :func:`stat_a_design`).
b  human-assigned, microbial and unresolved LFQ shares by region (``origin_fraction_cells.tsv``,
   ``four_region_tests.tsv``; recomputed from ``taxonomy_signal.tsv`` by :func:`stat_b_region_tests`).
c  own versus other stool distance (``own_stool_participants.tsv``, ``similarity_summary.tsv``;
   recomputed from the Zenodo OG profile matrix by :func:`stat_cd_profiles`).
d  variance partition (``orthogonal_partition.tsv``; recomputed by :func:`stat_cd_profiles`).
e  leading-species mosaic of 24 functions (``mosaic_donor_regions.tsv``, ``mosaic_abundance.tsv``,
   ``heatmap_function_key.tsv``; agreement counts by :func:`stat_e_mosaic`).
f  leading-species agreement under direct peptide assignment, all eligible and strict subset
   (``participant_leaders_*.tsv``; v30, the original protein grouping is no longer shown).
g  region-adjusted cross-omics Spearman correlations (``native_correlations.tsv``,
   ``functional_correlations.tsv``; recomputed from ``*_distances_long.tsv`` by
   :func:`stat_g_native` and :func:`stat_g_functional`).

The drawing code is moved from the original producer (``build_figure5_n15.py``, see
PROVENANCE.md) as it ran: one shared page canvas of 183 x 242.8 mm, panels b to g drawn on
it in this order, then composed with the panel-a artwork. Only file locations changed (they
are now repository-relative) and the producer's single script was split into one function
per panel. The page size, formerly read from an older Figure 5 PDF, is a constant here.

Text statistics in the Results that belong to Figure 5 (leader-share median, evenness) are
also computed here.
"""
import hashlib
import json
import math
import os
import re
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- choices (every value used below)
REG = ['duodenum', 'jejunum', 'ileum', 'colon', 'stool']   # internal ids of regions 1-4 and stool (Methods, v28d line 843)
REGL = ['Region 1', 'Region 2', 'Region 3', 'Region 4', 'Stool']
NP = 15                       # participants with all five sampling groups (Results, v28d line 149)
N_DNA_PARTICIPANTS = 13       # participants with complete DNA, protein and metabolite profiles (v28d line 161)
BOOT_SEED = 20261003          # participant bootstrap seed (not stated in the text; capscan_fig5bcd.py boot())
BOOT_DRAWS = 2000             # participant bootstrap resamples (Methods, v28d line 843)
OWN_STOOL_SEED = 20260910     # stool-label permutation seed (Supplementary Methods 8, v28d line 849)
OWN_STOOL_PERMUTATIONS = 99_999  # (v28d line 849)
PREVALENCES = [0.0, 0.1, 0.2]   # OG prevalence filters for c and d (v28d line 847)
CROSSOMICS_SEED = 20260907    # block-permutation seed, native and functional (not stated in the text; crossomics.py SEED)
CROSSOMICS_PERMUTATIONS = 1999   # whole-participant block permutations (Methods, v28d line 835)
FUNCTION_SCOPE = 'matched12'  # producer label of the 13-participant matched functional family (label predates participant 15)
NATIVE_ANALYSIS = 'sensitivity_all_features'   # native profiles, all features (Figure 5g)
EVEN_BAD_LABELS = ['Host', 'Mixed host/microbial', 'Conflicting species', 'No species annotation']  # excluded from named species (v28d line 849)
PANEL_F_MEAN_LABELS = ['46.4', '57.6']  # Fig. 5f means, all eligible and strict subset (v30 Results paragraph 81)
JITTER_SEED = 5               # horizontal jitter of panel b points (display only)

# page geometry of the producer (points); the producer read it from the v26b Figure 5 PDF
PAGE_W_PT = 518.7401733398438
PAGE_H_PT = 688.2833251953125
Y0 = 140.0                    # panel-a band ends here; panels b-g start below
# v31 layout (v33 figure): panel e region labels 1.0 pt under the heatmap, legend raised, panel f two-line x labels
E_LABEL_GAP_PT = 1.0
E_LEGEND_LY = 483.0
F_XLABELS = ['All\nfunctions', 'Well-resolved\nfunctions']

INK = '#4e4f51'
HOST, BLUE, ORANGE, SKY, GREY, LG = '#D55E00', '#0072B2', '#E69F00', '#56B4E9', '#7f7f7f', '#cfcfcf'
GREEN, VERM, PURPLE, YELLOW = '#009E73', '#D55E00', '#CC79A7', '#F0E442'
UNRES, RESID = '#d0d0d0', '#cccccc'

#: crop boxes for per-panel pixel checks, page points (x0, y0, x1, y1) on the full page
PANEL_BOXES_PT = {
    'b': (0.0, 140.0, 292.0, 292.0),
    'c': (292.0, 140.0, 405.0, 292.0),
    'd': (405.0, 140.0, PAGE_W_PT, 292.0),
    'e': (0.0, 292.0, PAGE_W_PT, 522.0),
    'f': (0.0, 522.0, 166.0, PAGE_H_PT),
    'g': (166.0, 522.0, PAGE_W_PT, PAGE_H_PT),
}

WORDS = {4: 'four', 5: 'five', 13: 'thirteen', 14: 'fourteen', 15: 'fifteen'}

#: Files above 20 MB read from build/zenodo/Fig5/ (see source_data/Fig5/README.md), with their sha256.
ZENODO_FILES = {
    'OG_other_unresolved.npz': 'e23d55a44eeb588caa9d4964ca6c707d934a567748b551eb3c3a326a353392ed',
    'OG_other_unresolved_axes.json': 'bf88b28c6e715b3ba2a2334f6257726e80287b73448dd9eff9ec9941116954ad',
    'function_species.npz': '4182b86aeb96de12e95c0572c039445a00e3fd70abe91b3addcf7c585e76817d',
}


def zenodo_ready(zdir):
    """True if every Zenodo file is present in ``zdir`` with the recorded sha256; raises on a wrong hash."""
    zdir = Path(zdir)
    if not all((zdir / f).is_file() for f in ZENODO_FILES):
        return False
    for f, h in ZENODO_FILES.items():
        assert sha256(zdir / f) == h, f'sha256 mismatch for {f}'
    return True


def sha256(path):
    """Return the hex sha256 of a file."""
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def as_word(n):
    """Return the English word the manuscript uses for a small count (e.g. 15 -> 'fifteen')."""
    return WORDS[int(n)]


# ---------------------------------------------------------------- fonts and style
def font_dir():
    """Return the Söhne font directory from FASTALAKE_FONT_DIR (or FL_STYLE_DIR), or None."""
    d = os.environ.get('FASTALAKE_FONT_DIR') or os.environ.get('FL_STYLE_DIR')
    return Path(d) if d else None


def fonts_available():
    """True if the two TrueType files the producer used are present (soehne-buch.ttf, soehne-halbfett.ttf)."""
    d = font_dir()
    return bool(d) and (d / 'soehne-buch.ttf').is_file() and (d / 'soehne-halbfett.ttf').is_file()


def setup_style():
    """Reset matplotlib to its defaults and apply the producer's rcParams and fonts.

    The producer did not use ``style/fl_style``; it started from matplotlib defaults and set
    the parameters below. Resetting first keeps the panels independent of any earlier style.

    Returns
    -------
    dict
        Font handles: ``HALB`` (letters), ``OBL`` (italic species names), ``have_fonts``.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager as fm
    matplotlib.rcdefaults()
    have = fonts_available()
    if have:
        d = font_dir()
        for f in ['soehne-buch.ttf', 'soehne-halbfett.ttf']:
            fm.fontManager.addfont(str(d / f))
        buch = fm.FontProperties(fname=str(d / 'soehne-buch.ttf')).get_name()
        halb = fm.FontProperties(fname=str(d / 'soehne-halbfett.ttf'))
    else:
        buch = 'DejaVu Sans'
        halb = fm.FontProperties(family='DejaVu Sans')
    obl = fm.FontProperties(fname=str(Path(matplotlib.get_data_path()) / 'fonts/ttf/DejaVuSans-Oblique.ttf'))
    plt.rcParams.update({'font.family': buch, 'pdf.fonttype': 42, 'font.size': 7.12, 'axes.unicode_minus': True,
                         'axes.linewidth': .8, 'xtick.major.width': .8, 'ytick.major.width': .8,
                         'xtick.major.size': 2.6, 'ytick.major.size': 2.6, 'text.color': INK,
                         'axes.labelcolor': INK, 'xtick.color': INK, 'ytick.color': INK, 'axes.edgecolor': INK})
    return dict(HALB=halb, OBL=obl, have_fonts=have)


class Canvas:
    """The shared page on which panels b-g are drawn (page points, origin top left).

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/Fig5``.
    """

    def __init__(self, sd):
        import matplotlib.pyplot as plt
        self.sd = Path(sd)
        self.fonts = setup_style()
        self.W, self.H = PAGE_W_PT, PAGE_H_PT
        self.h = self.H - Y0
        self.fig = plt.figure(figsize=(self.W / 72, self.h / 72))
        self.NUM, self.SRCS, self.data = {}, {}, {}

    def read(self, name):
        """Read a source_data table (as the producer did, pandas defaults) and record its sha256."""
        p = self.sd / name
        self.SRCS[name] = sha256(p)
        return pd.read_csv(p, sep='\t')

    def axpt(self, x0, y0, x1, y1):
        """Add axes at page coordinates (pt)."""
        W, h = self.W, self.h
        return self.fig.add_axes([x0 / W, 1 - (y1 - Y0) / h, (x1 - x0) / W, (y1 - y0) / h])

    def ftext(self, x, y, s, **kw):
        """Figure text at page coordinates."""
        return self.fig.text(x / self.W, 1 - (y - Y0) / self.h, s, **kw)

    def letter(self, x, y, s):
        """Panel letter as the producer drew it."""
        self.ftext(x, y, s, fontproperties=self.fonts['HALB'], fontsize=8.13, color='#262626', va='top', ha='left')


# ---------------------------------------------------------------- panels
def panel_b(cv):
    """Panel b: stacked mean LFQ shares per sampling group, participant points, region-vs-stool q.

    Parameters
    ----------
    cv : Canvas

    Returns
    -------
    dict
        q labels and plotted mean percentages.
    """
    o = cv.read('origin_fraction_cells.tsv')
    o['Microbial'] = o['Microbial above species'] + o['Microbial species']
    cells = o.rename(columns={'Human-assigned': 'Human', 'Unresolved origin': 'Unresolved'})
    assert sorted(cells.participant.unique()) == list(range(1, NP + 1)) and len(cells) == 5 * NP
    classes = [('Human', 'Human-assigned', HOST), ('Microbial', 'Microbial', BLUE), ('Unresolved', 'Unresolved', UNRES)]
    tests = cv.read('four_region_tests.tsv').set_index('region')
    assert len(tests) == 4 and (tests.n == NP).all()
    ax = cv.axpt(30.5, 158.4, 175, 253.45)
    mean = cells.groupby('region')[[c for c, _, _ in classes]].mean().loc[REG] * 100
    bottom = np.zeros(5)
    for c, lab, col in classes:
        ax.bar(range(5), mean[c], bottom=bottom, color=col, width=.68, label=lab); bottom += mean[c].values
    assert np.allclose(bottom, 100, atol=1e-6)
    rng = np.random.default_rng(JITTER_SEED)
    for i, r in enumerate(REG):
        v = cells[cells.region.eq(r)].sort_values('participant').Human * 100
        ax.scatter(i + rng.uniform(-.2, .2, len(v)), v, s=5, facecolor='white', edgecolor=INK, lw=.4, zorder=3)
    qlab = []
    for i, r in enumerate(REG[:4]):
        q = float(tests.loc[r, 'q'])
        s = 'q<0.001' if q < 0.001 else (f'q={q:.4f}' if q < 0.01 else f'q={q:.3f}')
        ax.text(i, 103, s, ha='center', va='bottom', fontsize=5.59); qlab.append(s)
    ax.text(4, 103, 'ref.', ha='center', va='bottom', fontsize=5.59)
    ax.set_xlim(-.55, 4.55); ax.set_ylim(0, 112); ax.set_yticks([0, 25, 50, 75, 100])
    ax.set_xticks(range(5), REGL, rotation=30, ha='right', rotation_mode='anchor'); ax.tick_params(axis='x', length=0)
    ax.spines[['top', 'right']].set_visible(False)
    ax.set_ylabel('Share of LFQ intensity (%)', labelpad=2)
    ax.legend(frameon=False, loc='upper left', bbox_to_anchor=(1.06, .9), handlelength=.9, labelspacing=.3, fontsize=6.1)
    cv.letter(34.6, 149.0, 'b')
    assert qlab == ['q=0.0028', 'q=0.0028', 'q=0.0028', 'q=0.038'], qlab
    cv.NUM['b'] = dict(q=qlab, mean_pct={c: [round(float(x), 2) for x in mean[c]] for c, _, _ in classes})
    cv.data.update(mean=mean, tests=tests, classes=classes)
    return cv.NUM['b']


def panel_c(cv):
    """Panel c: per-participant mean distance of the four capsule profiles to own versus other stools.

    Returns
    -------
    dict
        n_closer, means, P and feature count as plotted.
    """
    own = cv.read('own_stool_participants.tsv').query('prevalence==0').sort_values('participant')
    sim = cv.read('similarity_summary.tsv').query('prevalence==0').iloc[0]
    assert len(own) == NP
    ax = cv.axpt(296, 158.4, 400, 253.4)
    a_, b_ = own.own_stool_distance.to_numpy(), own.other_stool_mean_distance.to_numpy()
    for u, v in zip(a_, b_):
        ax.plot([0, 1], [u, v], color=LG, lw=.5, zorder=1)
    ax.scatter([0] * NP, a_, s=6, color=BLUE, lw=0, zorder=2); ax.scatter([1] * NP, b_, s=6, color=GREY, lw=0, zorder=2)
    for x, v in ((0, a_.mean()), (1, b_.mean())):
        ax.plot([x - .28, x + .28], [v, v], color='k', lw=1.1, zorder=3)
    ax.set_xlim(-.7, 1.7); ax.set_xticks([0, 1], ['Own\nstool', 'Other\nstools']); ax.tick_params(axis='x', length=0)
    assert a_.min() > .25 and b_.max() < .85
    ax.set_ylim(.25, .85); ax.set_yticks([.3, .4, .5, .6, .7, .8]); ax.spines[['top', 'right']].set_visible(False)
    ax.set_ylabel('Capsule–stool functional\nprofile distance', labelpad=2)
    cv.letter(299.7, 149.0, 'c')
    cv.NUM['c'] = dict(n_closer=int((a_ < b_).sum()), own_mean=float(a_.mean()), other_mean=float(b_.mean()), p=float(sim.p),
                       features=int(sim.features))
    assert cv.NUM['c']['n_closer'] == NP and abs(cv.NUM['c']['p'] - 1e-5) < 1e-12 and cv.NUM['c']['features'] == 73228
    cv.data.update(own=own)
    return cv.NUM['c']


def panel_d(cv):
    """Panel d: variance partition with and without stool (prevalence 0).

    Returns
    -------
    list of str
        The percentage labels drawn in the bars.
    """
    import matplotlib.pyplot as plt
    part = cv.read('orthogonal_partition.tsv').query('prevalence==0').set_index(['population', 'component']).fraction * 100
    ax = cv.axpt(433, 154.6, 515, 250.1)
    comps = [('Participant', BLUE, 'white'), ('Capsules versus stool', ORANGE, '#222222'),
             ('Between capsule regions', SKY, '#222222'), ('Residual', RESID, None)]
    bottom = np.zeros(2); dl = []
    for comp, col, tc in comps:
        vals = np.array([part[(pop, comp)] for pop in ['All five regions', 'Capsules only']])
        ax.bar([0, 1], vals, bottom=bottom, width=.62, color=col)
        for i, v in enumerate(vals):
            if v > 0 and tc:
                ax.text(i, bottom[i] + v / 2, f'{v:.1f}%', ha='center', va='center', fontsize=5.4, color=tc); dl.append(f'{v:.1f}%')
        bottom += vals
    assert np.allclose(bottom, 100, atol=1e-6)
    ax.set_ylim(0, 100); ax.set_yticks([0, 20, 40, 60, 80, 100]); ax.set_xlim(-.5, 1.5)
    ax.set_xticks([0, 1], ['Capsules\n+ stool', 'Capsules\nonly']); ax.tick_params(axis='x', length=0, pad=2)
    ax.tick_params(labelsize=5.7); ax.spines[['top', 'right']].set_visible(False)
    ax.set_ylabel('Share of profile variation (%)', labelpad=2, fontsize=6.0)
    W, h = cv.W, cv.h
    for x, y, lab, col in [(412, 274.0, 'Participant', BLUE), (468, 274.0, 'Residual', RESID),
                           (412, 284.0, 'Capsules vs stool', ORANGE), (468, 284.0, 'Capsule regions', SKY)]:
        cv.fig.add_artist(plt.Rectangle((x / W, 1 - (y + 2.6 - Y0) / h), 5 / W, 5 / h, transform=cv.fig.transFigure,
                                        facecolor=col, edgecolor='none'))
        cv.ftext(x + 7.5, y, lab, fontsize=5.4, va='center', color=INK)
    cv.letter(438.0, 149.0, 'd')
    assert dl == ['44.9%', '56.9%', '10.8%', '4.2%', '5.3%'], dl
    cv.NUM['d'] = dl
    cv.data.update(part=part)
    return dl


def spname(x):
    """Abbreviate a GTDB species label as the producer did: (italic part, roman suffix)."""
    g, _, sp = x.partition(' '); g = re.sub(r'_[A-Z]+$', '', g); sp = re.sub(r'_[A-Z]+$', '', sp)
    return ((g[0] + '. ' + sp), '') if sp and not sp.startswith('sp') else (g, ' ' + sp)


def panel_e(cv):
    """Panel e: leading annotated species per function, participant and region (24 functions).

    Returns
    -------
    dict
        The twelve coloured species and the agreement counts (as in :func:`stat_e_mosaic`).
    """
    import matplotlib.pyplot as plt
    fig, W, h = cv.fig, cv.W, cv.h
    mr = cv.read('mosaic_donor_regions.tsv'); key = cv.read('heatmap_function_key.tsv')
    ab = cv.read('mosaic_abundance.tsv').sort_values('row')
    assert set(key.KO) == set(ab.KO) and len(ab) == 24, 'the fixed 24 functions changed identity'
    key = key.set_index('KO').loc[ab.KO].reset_index()
    res = mr[mr.status.eq('Resolved')]; top = list(res.species.value_counts().index[:12])
    PAL = [BLUE, ORANGE, GREEN, VERM, SKY, PURPLE, YELLOW, '#1d2433', '#8c6d46', '#6fa8a1', '#b56576', '#9b9b3a']
    col = {sp: c for sp, c in zip(top, PAL)}; OTHER = '#c9c9c9'; NA = '#ffffff'
    MX0, MX1, MY0, MY1 = 151.0, 464.0, 311.0, 467.0
    axF = cv.axpt(MX0, MY0, MX1, MY1); axF.axis('off')
    cw, gap = 1.0, 0.9; x = 0; centers = {}; look = mr.set_index(['KO', 'participant', 'region'])
    for pnum in range(1, NP + 1):
        xs = []
        for rg in REG:
            for i, ko in enumerate(key.KO):
                row = look.loc[(ko, pnum, rg)]; c = NA if row.status != 'Resolved' else col.get(row.species, OTHER)
                axF.add_patch(plt.Rectangle((x, i - .5), cw, 1, facecolor=c, edgecolor='white', lw=.15))
            xs.append(x + cw / 2); x += cw
        centers[pnum] = float(np.mean(xs)); x += gap
    x -= gap
    axF.set_xlim(0, x); axF.set_ylim(23.5, -.5)
    for k, c in centers.items():
        axF.text(c, -1.0, str(k), ha='center', va='bottom', fontsize=5.5)
    for i, name in enumerate(key.display_name):
        axF.text(-.8, i, name, ha='right', va='center', fontsize=5.5)
    for j, ch in enumerate('1234S'):
        axF.text(j * cw + cw / 2, 23.5 + E_LABEL_GAP_PT / ((MY1 - MY0) / 24), ch, ha='center', va='top', fontsize=5.0, color='#8a8a8a')
    axFb = cv.axpt(477, MY0, 508, MY1)
    axFb.barh(range(24), ab.mean_share_percent, color=GREY, height=.6); axFb.set_ylim(23.5, -.5); axFb.set_yticks([])
    axFb.spines[['top', 'right', 'left']].set_visible(False)
    axFb.set_xlabel('Mean share of\ntotal LFQ (%)', fontsize=5.5, labelpad=2); axFb.set_xticks([0, 1, 2]); axFb.tick_params(labelsize=5.5)
    assert ab.mean_share_percent.max() < 2.2
    cv.letter(36.8, 296.0, 'e')
    # legend (3 rows, column-major); italic species names
    items = [(col[sp], spname(sp)) for sp in top] + [(OTHER, ('', 'Other species')), (NA, ('', 'Not resolved or not sampled'))]
    fig.canvas.draw(); rend = fig.canvas.get_renderer()

    def twidth(t):
        return t.get_window_extent(rend).width * 72 / fig.dpi
    LX, LY, LS = 82.0, E_LEGEND_LY, 8.2
    for cidx in range(5):
        widths = []
        for ridx in range(3):
            n = cidx * 3 + ridx
            if n >= len(items): break
            c, (it, ro) = items[n]; y = LY + ridx * LS
            fig.add_artist(plt.Rectangle((LX / W, 1 - (y + 2.6 - Y0) / h), 5 / W, 5 / h, transform=fig.transFigure,
                                         facecolor=c, edgecolor='#bdbdbd' if c == NA else 'none', lw=.4))
            xx = LX + 7.0; wsum = 0
            if it:
                t1 = cv.ftext(xx, y, it, fontproperties=cv.fonts['OBL'], fontsize=5.3, va='center', color=INK); w1 = twidth(t1); wsum += w1
            if ro:
                t2 = cv.ftext(xx + wsum, y, ro, fontsize=5.5, va='center', color=INK); wsum += twidth(t2)
            widths.append(wsum)
        LX += 7.0 + max(widths) + 9.0
    assert LX < W, LX
    st = stat_e_mosaic(mr, list(key.KO))
    cv.NUM['e'] = dict(top12=[''.join(spname(s)) for s in top], persistence=[st['persistence'], st['persistence_profiles']],
                       same_leader=[st['same_leader'], st['comparisons']])
    assert cv.NUM['e']['persistence'] == [111, 360] and cv.NUM['e']['same_leader'] == [726, 1440], cv.NUM['e']
    return cv.NUM['e']


def panel_f(cv):
    """Panel f (v30, production only): per-participant leading-species agreement under direct peptide assignment.

    One dot per participant (n = 15) for all eligible functions and for the strict subset, with
    the mean as a bar. Reads ``participant_leaders_resolved.tsv`` and
    ``participant_leaders_strict_both_methods.tsv`` (column ``agreement_simplified``).

    Parameters
    ----------
    cv : Canvas
        The shared page canvas.

    Returns
    -------
    list of str
        Mean labels drawn (all eligible, strict subset).
    """
    lr = cv.read('participant_leaders_resolved.tsv'); ls = cv.read('participant_leaders_strict_both_methods.tsv')
    assert len(lr) == NP and len(ls) == NP
    axE = cv.axpt(50, 542, 152, 605); Xe = [0, 1]; fmeans = []
    for (tab, x0) in ((lr, 0), (ls, 1)):
        axE.scatter([x0] * NP, tab.agreement_simplified * 100, s=6, color=BLUE, lw=0, zorder=2)
        v = tab.agreement_simplified.mean() * 100
        axE.plot([x0 - .28, x0 + .28], [v, v], color='k', lw=1.1, zorder=3)
        axE.text(x0 + .32, v, f'{v:.1f}', ha='left', va='center', fontsize=5.5, bbox=dict(facecolor='white', edgecolor='none', pad=.3))
        fmeans.append(f'{v:.1f}')
    axE.set_xlim(-.6, 1.75); axE.set_xticks(Xe, F_XLABELS, fontsize=6.5); axE.tick_params(axis='x', length=0)
    axE.set_ylabel('Same leading species,\ncapsule and own stool (%)', fontsize=6.5); axE.set_ylim(0, 100)
    axE.set_yticks([0, 25, 50, 75, 100]); axE.spines[['top', 'right']].set_visible(False)
    cv.letter(36.8, 527.0, 'f')
    assert fmeans == PANEL_F_MEAN_LABELS, fmeans
    cv.NUM['f'] = fmeans
    cv.data.update(lr=lr, ls=ls)
    return fmeans


LAYERS = ['DNA', 'Protein', 'Host protein', 'Metabolites']
KN = {'DNA': 'DNA', 'Protein': 'Non-host protein', 'Host protein': 'Host protein', 'Metabolites': 'Metabolites'}
KM = {'DNA': 'metaG', 'Protein': 'microbial_metaP', 'Host protein': 'host_metaP', 'Metabolites': 'metabolites'}


def panel_g(cv):
    """Panel g: lower-triangle Spearman correlations of region-adjusted distances (native, OG, KO).

    Returns
    -------
    list of str
        Cell labels in drawing order ('0.51*', ..., '—' for the untested KO host-metabolite pair).
    """
    from matplotlib.patches import Rectangle
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.colorbar import ColorbarBase
    W = cv.W
    nat = cv.read('native_correlations.tsv'); nat = nat[nat.cohort.eq('CAPSCAN') & nat.analysis.eq(NATIVE_ANALYSIS)]
    fun = cv.read('functional_correlations.tsv'); fun = fun[fun.cohort.eq('CAPSCAN') & fun.scope.eq(FUNCTION_SCOPE) & (fun.prevalence == 0)]
    assert len(nat) == 6 and len(fun) == 16 and (nat.n_participants == N_DNA_PARTICIPANTS).all() and (fun.n_donors == N_DNA_PARTICIPANTS).all()
    grows = []
    for lev in ['native', 'OG', 'KO']:
        src = nat if lev == 'native' else fun[fun.level.eq(lev)]
        maps = KN if lev == 'native' else KM
        tab = {frozenset((r.layer_a, r.layer_b)): r for r in src.itertuples()}
        for i, la in enumerate(LAYERS):
            for j, lb in enumerate(LAYERS[:i]):
                rec = tab.get(frozenset((maps[la], maps[lb])))
                if rec is None:
                    assert lev == 'KO' and {la, lb} == {'Host protein', 'Metabolites'}
                    grows.append(dict(representation=lev, row_layer=la, column_layer=lb, rho='', q='', label='—')); continue
                grows.append(dict(representation=lev, row_layer=la, column_layer=lb, rho=float(rec.rho), q=float(rec.q),
                                  label=f'{rec.rho:.2f}' + ('*' if rec.q < .05 else '')))
    gt = {(r['representation'], r['row_layer'], r['column_layer']): r for r in grows}
    s_ = W / 510.2362060546875
    GX0, GY0 = 169 * s_, 513 * s_; gw, gh = (510.2362060546875 - 169), 164
    axG = cv.axpt(GX0, GY0, GX0 + gw * s_, GY0 + gh * s_); axG.set_xlim(0, gw); axG.set_ylim(gh, 0); axG.axis('off')
    cmg = LinearSegmentedColormap.from_list('fl', ['#56B4E9', 'white', '#E69F00']); normg = Normalize(-1, 1)
    cell, n, gp, gx, gy = 18.3, 4, 14, 64, 35
    for c_, (lev, title) in enumerate([('native', 'Native'), ('OG', 'eggNOG OG'), ('KO', 'KEGG KO')]):
        x = gx + c_ * (n * cell + gp)
        axG.text(x + n * cell / 2, 30, title, fontsize=7 * s_, ha='center', va='center')
        for i, la in enumerate(LAYERS):
            if c_ == 0:
                axG.text(gx - 6, gy + (i + .5) * cell, f'{i+1}  ' + la, ha='right', va='center', fontsize=6.1 * s_)
            axG.text(x + (i + .5) * cell, gy + n * cell + 6, str(i + 1), ha='center', va='top', fontsize=6 * s_)
            for j, lb in enumerate(LAYERS[:i + 1]):
                if i == j:
                    rho, label = 1., '1.00'
                else:
                    r = gt[(lev, la, lb)]; label = r['label']; rho = 0. if label == '—' else r['rho']
                axG.add_patch(Rectangle((x + j * cell, gy + i * cell), cell, cell, facecolor=cmg(normg(rho)) if label != '—' else '#f5f5f5',
                                        edgecolor='white', lw=.5))
                axG.text(x + (j + .5) * cell, gy + (i + .5) * cell, label, ha='center', va='center', fontsize=6.0 * s_, color='#202020')
    cv.letter(GX0 + 3 * s_ - 1, GY0 + 5.5 * s_, 'g')
    bx, by = gw - 99, 8
    cax = cv.axpt(GX0 + bx * s_, GY0 + by * s_, GX0 + (bx + 85) * s_, GY0 + (by + 4) * s_)
    cb = ColorbarBase(cax, cmap=cmg, norm=normg, orientation='horizontal', ticks=[-1, 0, 1])
    cb.outline.set_visible(False); cb.ax.tick_params(labelsize=5.5, length=1.5, pad=1)
    glabels = [r['label'] for r in grows]
    assert [gt[('native', 'Protein', 'DNA')]['label'], gt[('native', 'Host protein', 'Protein')]['label'],
            gt[('native', 'Metabolites', 'Protein')]['label']] == ['0.51*', '0.55*', '0.38*'], glabels
    cv.NUM['g'] = glabels
    cv.data.update(grows=grows)
    return glabels


def save_body(cv, path):
    """Save panels b-g (the producer's ``panels_bg_n15.pdf``) and close the canvas figure."""
    import matplotlib.pyplot as plt
    cv.fig.savefig(path); plt.close(cv.fig)


def compose(panel_a_pdf, body_pdf, out_pdf, out_png):
    """Compose the full page: panel-a artwork in the top band, panels b-g below (producer step).

    Parameters
    ----------
    panel_a_pdf, body_pdf : pathlib.Path
        ``artwork/panel_a.pdf`` and the output of :func:`save_body`.
    out_pdf, out_png : pathlib.Path
        Composite PDF and its 300-dpi rasterisation.
    """
    import fitz
    from PIL import Image
    out = fitz.open(); op = out.new_page(width=PAGE_W_PT, height=PAGE_H_PT)
    with fitz.open(body_pdf) as b: op.show_pdf_page(fitz.Rect(0, Y0, PAGE_W_PT, PAGE_H_PT), b, 0)
    with fitz.open(panel_a_pdf) as a: op.show_pdf_page(fitz.Rect(0, 0, PAGE_W_PT, Y0), a, 0)
    out.save(out_pdf, garbage=4, deflate=True)
    pm = op.get_pixmap(dpi=300, alpha=False); pm.save(out_png)
    im = Image.open(out_png); im.save(out_png, dpi=(300, 300))
    return dict(png_px=[pm.width, pm.height])


def check_composite(pdf):
    """Producer's text checks on the composite: minimum font size, no stale and all required strings, overlaps.

    Returns
    -------
    dict
        Font sizes, font names and overlapping text spans.
    """
    import fitz
    d = fitz.open(pdf)[0]
    spans = [sp for b in d.get_text('dict')['blocks'] for l in b.get('lines', []) for sp in l['spans'] if sp['text'].strip()]
    sizes = sorted({round(sp['size'], 2) for sp in spans}); fonts = sorted({sp['font'] for sp in spans})
    assert min(sizes) >= 5.0, sizes
    alltext = ' '.join(sp['text'] for sp in spans)
    for stale in ['331 metaproteomes', '274 capsule', '14/14', 'q=0.069', '12 participants', 'Rows and columns', 'CapScan · region',
                  'Own stool closer', '70 acquisitions', '61 acquisitions']:
        assert stale not in alltext, stale
    for need in ['347 metaproteomes', '290 capsule · 57 stool', '74 acquisitions', '77 acquisitions', '65 acquisitions', 'q=0.038']:
        assert need in alltext, need
    boxes = [fitz.Rect(sp['bbox']) for sp in spans]
    over = [(spans[i]['text'], spans[j]['text']) for i in range(len(boxes)) for j in range(i + 1, len(boxes))
            if (boxes[i] & boxes[j]).get_area() > 0.25 * min(boxes[i].get_area(), boxes[j].get_area()) > 0]
    return dict(font_sizes_pt=sizes, fonts=fonts, text_overlaps=over)


def plotted_values(cv):
    """Every plotted value of panels b, c, d, f and g as one long table (the producer's plotted-values file)."""
    mean, tests, classes = cv.data['mean'], cv.data['tests'], cv.data['classes']
    own, part, lr, ls, grows = cv.data['own'], cv.data['part'], cv.data['lr'], cv.data['ls'], cv.data['grows']
    rows = [dict(panel='b', item='mean_share_pct', group=g_, component=c, value=mean.loc[g_, c]) for g_ in REG for c, _, _ in classes]
    rows += [dict(panel='b', item='q_vs_stool', group=g_, component='Human', value=tests.loc[g_, 'q']) for g_ in REG[:4]]
    rows += [dict(panel='c', item=k, group=int(r.participant), component='', value=getattr(r, k)) for r in own.itertuples()
             for k in ('own_stool_distance', 'other_stool_mean_distance')]
    rows += [dict(panel='d', item='variation_share_pct', group=p, component=c, value=v) for (p, c), v in part.items()]
    rows += [dict(panel='f', item=t, group=int(r.participant), component=m, value=getattr(r, m)) for t, T in (('resolved', lr), ('strict', ls))
             for r in T.itertuples() for m in ('agreement_simplified',)]
    rows += [dict(panel='g', item=r['representation'], group=r['row_layer'], component=r['column_layer'], value=r['rho'], q=r['q']) for r in grows]
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- pixel-check helpers
def rasterise(pdf, dpi=300):
    """Rasterise page 1 of a PDF to an RGB array (PyMuPDF, as tools/check_pixels.py does)."""
    import fitz
    with fitz.open(pdf) as doc:
        pix = doc[0].get_pixmap(dpi=dpi, alpha=False, colorspace=fitz.csRGB)
        return np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, 3).copy()


def crop_panels(body_pdf, out_dir, dpi=300):
    """Crop each panel b-g from a rasterised panels-b-g PDF and save ``panel_<x>.png``.

    The PDF page starts at page y = Y0, so box y coordinates are shifted by Y0.

    Returns
    -------
    dict
        panel -> written path.
    """
    from PIL import Image
    img = rasterise(body_pdf, dpi); s = dpi / 72; out = {}
    for k, (x0, y0, x1, y1) in PANEL_BOXES_PT.items():
        c = img[int(round((y0 - Y0) * s)):int(round((y1 - Y0) * s)), int(round(x0 * s)):int(round(x1 * s))]
        p = Path(out_dir) / f'panel_{k}.png'; Image.fromarray(c).save(p); out[k] = p
    return out


# ---------------------------------------------------------------- statistics
def load_design(sd):
    """Public CapScan design: one row per acquisition (347), indexed by sample_id."""
    return pd.read_csv(Path(sd) / 'design_347.tsv', sep='\t').set_index('sample_id')


def stat_a_design(sd):
    """Counts typed in panel a and printed in the Results (v28d line 147).

    Statistical unit: acquisition (counts), specimen (308), participant (15).

    Returns
    -------
    dict
        acquisitions, per-region counts, capsule and stool totals, specimens, participants.
    """
    d = load_design(sd)
    vc = d.region.value_counts()
    return dict(acquisitions=int(len(d)), counts={r: int(vc[r]) for r in REG}, capsule=int(vc[REG[:4]].sum()), stool=int(vc['stool']),
                specimens=int(d.specimen_id.nunique()), participants=int(d.participant.nunique()), capsule_regions=len(REG) - 1,
                og_positive_acquisitions=int(d.og_profile_positive.sum()))


def boot(x, seed=BOOT_SEED, B=BOOT_DRAWS):
    """Percentile bootstrap 95% interval of the mean, resampling participants (capscan_fig5bcd.py boot()).

    Parameters
    ----------
    x : array-like
        One value per participant (n = 15).

    Returns
    -------
    tuple of float
        2.5th and 97.5th percentiles of 2,000 resampled means. The same seed is used for
        every call, so all intervals share one set of resampled participant indices.
    """
    x = np.asarray(x); rng = np.random.default_rng(seed); idx = rng.integers(0, len(x), (B, len(x)))
    m = x[idx].mean(1); return float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def stat_b_region_tests(sd):
    """Human-assigned LFQ share, each capsule region versus own stool (Figure 5b, v28d line 149).

    Statistical unit: participant, n = 15 (paired). Acquisition shares (Human-assigned /
    total LFQ) are averaged within specimen, then within participant and region. Test: exact
    two-sided paired sign-flip over all 2^15 = 32,768 sign patterns of the per-participant
    differences (statistic: sum of differences); 95% CI: percentile bootstrap of the mean
    difference (2,000 participant resamples, seed BOOT_SEED); Benjamini-Hochberg over the four
    regions (scipy.stats.false_discovery_control).

    Returns
    -------
    tuple of pandas.DataFrame
        (cells, tests): 75 participant-region origin fractions, and the four region tests.
    """
    from scipy.stats import false_discovery_control
    design = load_design(sd)
    roster = design.index[design.og_profile_positive.astype(bool)]
    obs = design.loc[roster].groupby('participant').region.nunique()
    parts = sorted(obs.index[obs.eq(5)]); n = len(parts)
    ix = pd.MultiIndex.from_product([parts, REG], names=['participant', 'region'])
    d = design.loc[roster]
    tx = pd.read_csv(Path(sd) / 'taxonomy_signal.tsv', sep='\t').set_index('sample')
    tf = tx.loc[roster].drop(columns='total_LFQ').div(tx.loc[roster, 'total_LFQ'], axis=0)
    tc = tf.groupby([d.participant, d.region, d.specimen_id]).mean().groupby(level=[0, 1]).mean().loc[ix]
    w = tc['Human-assigned'].unstack().loc[parts, REG] * 100
    signs = 2 * ((np.arange(2**n, dtype=np.uint64)[:, None] >> np.arange(n, dtype=np.uint64)) & 1).astype(np.int8) - 1
    tests = []
    for reg in REG[:4]:
        delta = (w[reg] - w.stool).to_numpy(); p = float(np.mean(abs(signs @ delta) >= abs(delta.sum()) - 1e-12)); lo, hi = boot(delta)
        tests.append(dict(region=reg, n=n, mean_region=float(w[reg].mean()), mean_stool=float(w.stool.mean()), mean_difference=float(delta.mean()),
                          ci_low=lo, ci_high=hi, p=p, sum_of_differences=float(delta.sum()), patterns=int(2**n)))
    for r, q in zip(tests, false_discovery_control([t['p'] for t in tests])): r['q'] = float(q)
    return tc, pd.DataFrame(tests)


def frame(folder, key):
    """Load a sparse profile matrix (``<key>.npz`` + ``<key>_axes.json``) as a dense DataFrame."""
    from scipy import sparse
    a = json.loads((Path(folder) / (key + '_axes.json')).read_text()); v = sparse.load_npz(Path(folder) / (key + '.npz'))
    return pd.DataFrame(v.toarray(), index=a['samples'], columns=a['features'])


def partition(h, n):
    """Orthogonal sums-of-squares partition (participant, capsule-vs-stool, between capsule regions, residual).

    Statistical unit: participant-region profile (75 with stool, 60 capsule only). Descriptive.
    """
    cube = np.asarray(h).reshape(n, 5, -1); rows = []
    for label, xx in [('All five regions', cube), ('Capsules only', cube[:, :4])]:
        nn, r, k = xx.shape; grand = xx.mean((0, 1)); pm = xx.mean(1) - grand; rm = xx.mean(0) - grand
        res = xx - grand - pm[:, None] - rm[None]; ss = ((xx - grand)**2).sum()
        d = {'Participant': r*(pm**2).sum()/ss, 'Residual': (res**2).sum()/ss}
        if r == 5:
            means = xx.mean(0); contrast = np.array([1, 1, 1, 1, -4])/np.sqrt(20)
            d['Capsules versus stool'] = nn*((contrast[:, None]*means).sum(0)**2).sum()/ss
            d['Between capsule regions'] = nn*((means[:4] - means[:4].mean(0))**2).sum()/ss
        else: d.update({'Capsules versus stool': 0., 'Between capsule regions': nn*(rm**2).sum()/ss})
        assert np.isclose(sum(d.values()), 1, atol=1e-10)
        rows.extend(dict(population=label, component=k, fraction=float(v)) for k, v in d.items())
    return pd.DataFrame(rows)


def stat_cd_profiles(sd, zenodo_dir):
    """Own-stool test (Figure 5c) and variance partition (Figure 5d) from the acquisition OG profiles.

    Statistical unit: participant, n = 15; 75 participant-region profiles. Acquisition OG
    profiles (non-host OG allocation, 347 acquisitions; 346 positive) are normalised to sum
    1, averaged within specimen then participant-region, filtered by prevalence, square-root
    transformed (Hellinger). Own distance = mean Euclidean distance from the four capsule
    profiles to own stool; other = mean over the 14 other stools. Test: two-sided
    permutation of stool-to-participant assignments (99,999 permutations, seed
    OWN_STOOL_SEED, P = (1 + #|null| >= |obs|) / 100,000). 95% CI: percentile bootstrap of
    the mean advantage (2,000 participant resamples, seed BOOT_SEED).

    Parameters
    ----------
    sd : pathlib.Path
    zenodo_dir : pathlib.Path
        Folder holding ``OG_other_unresolved.npz`` and ``OG_other_unresolved_axes.json``.

    Returns
    -------
    tuple of pandas.DataFrame
        (similarity_summary, own_stool_participants, orthogonal_partition, roster_check)
    """
    from scipy.spatial.distance import cdist
    design = load_design(sd)
    og = frame(zenodo_dir, 'OG_other_unresolved'); assert set(design.index) == set(og.index), 'profile/design mismatch'
    pos = og.loc[design.index].sum(1) > 0
    roster_check = bool((pos.to_numpy() == design.og_profile_positive.astype(bool).to_numpy()).all())
    roster = design.index[pos]
    obs = design.loc[roster].groupby('participant').region.nunique()
    parts = sorted(obs.index[obs.eq(5)]); n = len(parts)
    ix = pd.MultiIndex.from_product([parts, REG], names=['participant', 'region'])
    d = design.loc[roster]
    v_ = og.loc[roster]; rel = v_.div(v_.sum(1), axis=0)
    cc = rel.groupby([d.participant, d.region, d.specimen_id]).mean().groupby(level=[0, 1]).mean().loc[ix]
    del og, v_, rel
    x = cc.to_numpy()
    perms = np.argsort(np.random.default_rng(OWN_STOOL_SEED).random((OWN_STOOL_PERMUTATIONS, n)), axis=1)
    rows, don, parts_out = [], [], []
    for prev in PREVALENCES:
        keep = (x > 0).mean(0) >= prev; v = x[:, keep]; h = np.sqrt(v/v.sum(1, keepdims=True)); cube = h.reshape(n, 5, -1)
        dd = np.stack([cdist(cube[i, :4], cube[:, 4]).mean(0) for i in range(n)]); own = np.diag(dd); other = (dd.sum(1) - own)/(n - 1)
        obs_ = float((other - own).mean()); pdg = dd[np.arange(n)[None, :], perms].mean(1); null = (dd.sum(1).mean() - n*pdg)/(n - 1)
        p = (1 + int((abs(null) >= abs(obs_) - 1e-14).sum()))/(OWN_STOOL_PERMUTATIONS + 1)
        lo, hi = boot(other - own)
        rows.append(dict(prevalence=prev, participants=n, features=int(keep.sum()), positive_features=int((v > 0).any(0).sum()), own_closer=int((own < other).sum()),
                         mean_own=float(own.mean()), mean_other=float(other.mean()), mean_advantage=obs_, ci_low=lo, ci_high=hi, p=p,
                         extreme_permutations=int((abs(null) >= abs(obs_) - 1e-14).sum())))
        don.extend(dict(prevalence=prev, participant=parts[i], own_stool_distance=own[i], other_stool_mean_distance=other[i]) for i in range(n))
        pt = partition(h, n); pt['prevalence'] = prev; parts_out.append(pt)
    return pd.DataFrame(rows), pd.DataFrame(don), pd.concat(parts_out, ignore_index=True), roster_check


def stat_e_mosaic(mr, kos):
    """Leading-species agreement in the 24-function mosaic (Figure 5e, v28d lines 159 and 787).

    Statistical unit: function-participant-region comparison (same_leader, n = 1,440) and
    participant-function profile (persistence, n = 360). Descriptive counts.

    Parameters
    ----------
    mr : pandas.DataFrame
        ``mosaic_donor_regions.tsv``.
    kos : list of str
        The 24 displayed KOs.
    """
    look = mr.set_index(['KO', 'participant', 'region'])
    same = []
    for ko in kos:
        for pn in range(1, NP + 1):
            st = look.loc[(ko, pn, 'stool')]
            if st.status != 'Resolved': continue
            caps = [rg for rg in REG[:4] if look.loc[(ko, pn, rg)].status == 'Resolved']
            if caps: same.append(all(look.loc[(ko, pn, rg)].species == st.species for rg in caps))
    pairs = [(look.loc[(ko, pn, rg)].species == look.loc[(ko, pn, 'stool')].species)
             for ko in kos for pn in range(1, NP + 1) for rg in REG[:4]
             if look.loc[(ko, pn, rg)].status == 'Resolved' and look.loc[(ko, pn, 'stool')].status == 'Resolved']
    return dict(same_leader=int(sum(pairs)), comparisons=len(pairs), same_leader_fraction=sum(pairs) / len(pairs),
                persistence=int(sum(same)), persistence_profiles=len(same), persistence_fraction=sum(same) / len(same),
                functions=len(kos))


def stat_f_leaders(sd):
    """Mean participant leading-species agreement under direct peptide assignment (Figure 5f, v30).

    Unit: participant, n = 15. Descriptive means, no test.

    Parameters
    ----------
    sd : Path
        ``source_data/Fig5``.

    Returns
    -------
    dict
        ``participants``, ``resolved_direct`` (all eligible) and ``strict_direct`` (strict subset), as fractions.
    """
    lr = pd.read_csv(Path(sd) / 'participant_leaders_resolved.tsv', sep='\t')
    ls = pd.read_csv(Path(sd) / 'participant_leaders_strict_both_methods.tsv', sep='\t')
    return dict(participants=int(len(lr)), resolved_direct=float(lr.agreement_simplified.mean()),
                strict_direct=float(ls.agreement_simplified.mean()))


def stat_species_label_coverage(sd):
    """Median species-resolved fraction of non-host functional intensity over eligible function-participant cells (v28d line 159).

    Unit: function-participant cell (21,168 cells, 1,467 functions, 15 participants). Descriptive.
    """
    c = pd.read_csv(Path(sd) / 'contributor_cells.tsv', sep='\t')
    return dict(cells=int(len(c)), functions=int(c.KO.nunique()), participants=int(c.participant.nunique()),
                median_resolved_fraction=float(c.species_resolved_fraction.median()))


def evenness_metrics(v):
    """Leading-species share and Pielou evenness for each column of a species x sample matrix (C4 metrics())."""
    v = np.asarray(v, float); tot = v.sum(0)
    p = np.divide(v, tot[None, :], out=np.zeros_like(v), where=tot[None, :] > 0)
    S = (p > 0).sum(0); logs = np.zeros_like(p); np.log(p, out=logs, where=p > 0)
    entropy = -(p * logs).sum(0)
    even = np.divide(entropy, np.log(np.maximum(S, 1)), out=np.full(len(tot), np.nan), where=S > 1)
    top1 = np.where(tot > 0, np.sort(p, axis=0)[-1], np.nan)
    return dict(top1=top1, pielou=even, richness=S.astype(float), resolved_mass=tot)


def stat_evenness_per_acquisition(sd, zenodo_dir):
    """Recompute the per-acquisition evenness table from ``function_species.npz`` (C4 load_a33 + per_acquisition).

    Unit: acquisition (346 with species signal; one roster acquisition without signal kept as undefined).
    """
    design = load_design(sd).reset_index()
    z = np.load(Path(zenodo_dir) / 'function_species.npz', allow_pickle=True)
    vals, sp, samples = z['values'], z['species'], z['samples']
    named = ~np.isin(sp, EVEN_BAD_LABELS); nonhost = ~np.isin(sp, ['Host', 'Mixed host/microbial'])
    names, inv = np.unique(sp[named], return_inverse=True)
    taxa = np.zeros((len(names), vals.shape[1])); np.add.at(taxa, inv, vals[named])
    denom = vals[nonhost].sum(0)
    m = evenness_metrics(taxa)
    d = design.set_index('sample_id').loc[list(samples)].reset_index()
    out = d[['sample_id', 'participant', 'region', 'specimen_id']].copy()
    out['source'] = 'a33'
    out['compartment'] = np.where(out.region == 'stool', 'stool', 'fluid')
    for k in ['top1', 'pielou', 'richness', 'resolved_mass']:
        out[k] = m[k]
    out['species_resolved_fraction'] = np.divide(m['resolved_mass'], denom, out=np.full(len(denom), np.nan), where=denom > 0)
    missing = sorted(set(design.sample_id) - set(samples))
    if missing:
        add = design[design.sample_id.isin(missing)][['sample_id', 'participant', 'region', 'specimen_id']].copy()
        add['source'] = 'a33'; add['compartment'] = np.where(add.region == 'stool', 'stool', 'fluid')
        for k in ['top1', 'pielou', 'species_resolved_fraction']:
            add[k] = np.nan
        add['richness'] = 0.0; add['resolved_mass'] = 0.0
        out = pd.concat([out, add], ignore_index=True)
    return out


def signflip(v):
    """Exact two-sided paired sign-flip on the sum (C4 signflip())."""
    v = np.asarray(v, float); v = v[np.isfinite(v)]; n = len(v)
    signs = (((np.arange(2 ** n, dtype=np.uint32)[:, None] >> np.arange(n)) & 1).astype(np.int8) * 2 - 1)
    return float(np.mean(np.abs(signs @ v) >= abs(v.sum()) - 1e-12))


def stat_evenness(acq):
    """Intestinal fluid versus stool: leading-species share and Pielou evenness (v28d line 159).

    Statistical unit: participant, n = 15 paired. Per acquisition, metrics over named-species
    signal; per participant, the median over fluid acquisitions and over stool acquisitions.
    Test: Wilcoxon signed-rank, two-sided, exact (scipy.stats.wilcoxon, method='exact',
    zero_method='wilcox'); exact sign-flip reported alongside. No multiplicity correction.

    Returns
    -------
    tuple of pandas.DataFrame
        (summary, participant_medians)
    """
    from scipy.stats import wilcoxon
    rows = []; part_rows = []
    for metric in ['top1', 'pielou']:
        good = acq[np.isfinite(acq[metric])]
        pm = good.groupby(['participant', 'compartment'])[metric].median().unstack()
        both = pm.dropna(subset=['fluid', 'stool'])
        diff = (both.fluid - both.stool).to_numpy()
        w = wilcoxon(both.fluid, both.stool, alternative='two-sided', method='exact')
        for p_, r in pm.iterrows():
            part_rows.append(dict(metric=metric, participant=p_, fluid_median=r.get('fluid', np.nan), stool_median=r.get('stool', np.nan)))
        rows.append(dict(metric=metric, acq_n_fluid=int((good.compartment == 'fluid').sum()), acq_n_stool=int((good.compartment == 'stool').sum()),
                         acq_n_undefined=int(len(acq) - len(good)), paired_participants=int(len(both)),
                         participant_median_fluid=float(both.fluid.median()), participant_median_stool=float(both.stool.median()),
                         median_paired_difference=float(np.median(diff)), mean_paired_difference=float(diff.mean()),
                         participants_fluid_gt_stool=int((diff > 0).sum()), participants_fluid_lt_stool=int((diff < 0).sum()),
                         wilcoxon_W=float(w.statistic), wilcoxon_p_exact=float(w.pvalue),
                         signflip_p_exact=signflip(diff), signflip_patterns=2 ** len(diff)))
    return pd.DataFrame(rows), pd.DataFrame(part_rows)


def distances_from_long(df, n_participants=N_DNA_PARTICIPANTS):
    """Rebuild square distance matrices from a long upper-triangle table.

    Returns
    -------
    numpy.ndarray
        Symmetric (5n x 5n) matrix with zero diagonal, profile order = participants x REG.
    """
    n = 5 * n_participants; iu = np.triu_indices(n, 1)
    assert len(df) == len(iu[0])
    m = np.zeros((n, n)); m[iu] = df.distance.to_numpy(); m = m + m.T
    return m


def bh(values):
    """Benjamini-Hochberg q values (native_pyramids.py bh())."""
    v = np.asarray(values); o = np.argsort(v); q = np.empty(len(v)); q[o] = np.minimum(1, np.minimum.accumulate((v[o]*len(v)/np.arange(1, len(v)+1))[::-1])[::-1]); return q


def stat_g_native(sd):
    """Native profile concordance (Figure 5g left, v28d line 161).

    Statistical unit: participant, n = 13 (65 participant-region profiles). Spearman
    correlation of the upper triangles of two region-adjusted Euclidean distance matrices;
    two-sided P from 1,999 whole-participant block permutations (seed CROSSOMICS_SEED,
    P = (1 + #|null| >= |rho| - 1e-12) / 2,000); Benjamini-Hochberg over the six layer pairs.
    Ported from native_pyramids.py correlations().
    """
    from scipy.spatial.distance import squareform
    from scipy.stats import rankdata, spearmanr
    long = pd.read_csv(Path(sd) / 'native_distances_long.tsv', sep='\t', float_precision='round_trip')
    long = long[long.analysis.eq(NATIVE_ANALYSIS)]
    distances = {k: distances_from_long(long[long.layer.eq(k)]) for k in ['DNA', 'Non-host protein', 'Host protein', 'Metabolites']}
    n = len(next(iter(distances.values()))); ix = np.triu_indices(n, 1)
    rng = np.random.default_rng(CROSSOMICS_SEED)
    blocks = np.arange(n).reshape(-1, 5); perms = np.array([blocks[rng.permutation(len(blocks))].ravel() for _ in range(CROSSOMICS_PERMUTATIONS)])
    ranks = {k: rankdata(v[ix]) for k, v in distances.items()}; norm = {}
    for k, r in ranks.items():
        centered = r - r.mean(); scale = np.linalg.norm(centered); assert scale > 1e-10
        norm[k] = centered/scale
    rows = []
    for a, b in combinations(distances, 2):
        rho = float(norm[a] @ norm[b]); checkrho = float(spearmanr(distances[a][ix], distances[b][ix]).statistic)
        assert abs(rho - checkrho) < 1e-12
        bmat = squareform(norm[b]); extreme = 0
        for chunk in np.array_split(perms, 20):
            null = bmat[chunk[:, ix[0]], chunk[:, ix[1]]] @ norm[a]
            extreme += int((abs(null) >= abs(rho) - 1e-12).sum())
        rows.append(dict(cohort='CAPSCAN', analysis=NATIVE_ANALYSIS, layer_a=a, layer_b=b, rho=rho, p=(1 + extreme)/(CROSSOMICS_PERMUTATIONS + 1),
                         extreme=extreme, permutations=CROSSOMICS_PERMUTATIONS, n_profiles=n, n_participants=n//5))
    table = pd.DataFrame(rows); table['q'] = bh(table.p)
    return table


def stat_g_functional(sd):
    """Functional profile concordance, matched 13-participant family (Figure 5g OG and KO; v28d lines 161, 835).

    Statistical unit: participant, n = 13 (65 profiles). Spearman correlation of distance
    upper triangles (ranks of pdist of region-adjusted Hellinger profiles); two-sided P
    from 1,999 whole-participant block permutations (seed CROSSOMICS_SEED, fresh generator
    per representation, P = (1 + #|null| >= |rho|) / 2,000); Benjamini-Hochberg
    (statsmodels multipletests 'fdr_bh') over the 16 OG, KO and pathway comparisons. Ported
    from check_prevalence.py evaluate(); the KO and pathway host-metabolite pairs are not tested.
    """
    from scipy.spatial.distance import squareform
    from scipy.stats import rankdata
    from statsmodels.stats.multitest import multipletests
    long = pd.read_csv(Path(sd) / 'functional_distances_long.tsv', sep='\t', float_precision='round_trip')
    rows = []
    for lev in ['OG', 'KO', 'Pathway']:
        sub = long[long.level.eq(lev)]
        mats = {k: distances_from_long(sub[sub.layer.eq(k)]) for k in ['metaG', 'microbial_metaP', 'host_metaP', 'metabolites']}
        n = len(next(iter(mats.values()))); rng = np.random.default_rng(CROSSOMICS_SEED)
        blocks = np.arange(n).reshape(n//5, 5); perms = np.array([blocks[rng.permutation(n//5)].ravel() for _ in range(CROSSOMICS_PERMUTATIONS)])
        ix = np.triu_indices(n, 1); ranks = {}
        for k, m in mats.items():
            dist = m[ix]; r = rankdata(dist); r -= r.mean(); r /= np.linalg.norm(r); ranks[k] = r
        for a, b in combinations(mats, 2):
            if lev != 'OG' and {a, b} == {'host_metaP', 'metabolites'}: continue
            rho = float(ranks[a] @ ranks[b])
            bm = squareform(ranks[b]); ext = 0
            for start in range(0, CROSSOMICS_PERMUTATIONS, 200):
                pp = perms[start:start+200]; null = bm[pp[:, ix[0]], pp[:, ix[1]]] @ ranks[a]; ext += np.count_nonzero(np.abs(null) >= abs(rho))
            rows.append(dict(cohort='CAPSCAN', scope=FUNCTION_SCOPE, level=lev, layer_a=a, layer_b=b, prevalence=0.0, n_units=n,
                             n_donors=n//5, rho=rho, p=(1+ext)/(1+CROSSOMICS_PERMUTATIONS), extreme=int(ext)))
    res = pd.DataFrame(rows); assert len(res) == 16
    res['q'] = multipletests(res.p, method='fdr_bh')[1]
    return res


def compare_tables(new, ref, keys, cols):
    """Return the largest absolute difference of numeric ``cols`` between two tables joined on ``keys``."""
    m = new.merge(ref, on=keys, suffixes=('_new', '_ref'), how='outer', indicator=True)
    assert (m['_merge'] == 'both').all(), 'row sets differ'
    return {c: float(np.max(np.abs(m[c + '_new'].astype(float) - m[c + '_ref'].astype(float)))) for c in cols}


#: NUMBERS.tsv rows in order: key into the notebook's results dictionary ``R``.
NUMBER_KEYS = [
    'L147_acquisitions', 'L147_specimens', 'L147_participants', 'L147_capsule_regions_word',
    'L149_region1_share', 'L149_region4_share', 'L149_stool_share', 'L149_participants', 'L149_groups_word',
    'L149_diff_r1', 'L149_diff_r2', 'L149_diff_r3', 'L149_ci_lo_r1', 'L149_ci_hi_r1', 'L149_ci_lo_r2', 'L149_ci_hi_r2',
    'L149_ci_lo_r3', 'L149_ci_hi_r3', 'L149_q_r1', 'L149_q_r2', 'L149_q_r3', 'L149_diff_r4', 'L149_ci_lo_r4', 'L149_ci_hi_r4', 'L149_q_r4',
    'L153_participants_b', 'L153_features', 'L153_functions', 'L153_participants_g',
    'L155_closer', 'L155_advantage', 'L155_ci_lo', 'L155_ci_hi', 'L155_p',
    'L155_participant', 'L155_sampling_group', 'L155_capsule_stool', 'L155_among_region', 'L155_caps_participant', 'L155_caps_region',
    'L157_functions', 'L157_resolved_direct', 'L157_participants', 'L157_strict_direct',
    'L159_functions', 'L159_same_leader_pct', 'L159_same_leader', 'L159_comparisons', 'L159_persistence', 'L159_persistence_profiles',
    'L159_persistence_pct', 'L159_label_coverage', 'L159_top1_fluid', 'L159_top1_stool', 'L159_pielou_fluid', 'L159_pielou_stool',
    'L159_expected_direction', 'L159_paired_n', 'L159_wilcoxon_top1', 'L159_wilcoxon_pielou',
    'L161_participants', 'L161_rho_dna', 'L161_rho_host', 'L161_rho_met', 'L161_q_dna', 'L161_q_host', 'L161_q_met',
    'L281_features', 'L281_profiles',
    'L787_persistence_pct',
    'L847_profiles', 'L847_features', 'L847_features_10', 'L847_features_20', 'L847_p_sensitivity',
    'L849_capsule_stool_exact', 'L849_among_region_exact', 'L849_caps_profiles', 'L849_caps_region_exact', 'L849_sampling_group_exact',
]
