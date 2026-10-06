"""Supplementary Figure 6: plotting and statistics functions.

Panels
------
a, b, c
    Spectrum-split experiment. Protein databases selected from all predictions (FULL),
    odd-scan predictions (JOINT, called HALF_SELF or "self" in the tables) or even-scan
    predictions (DISJOINT), each searched against the same odd-scan MS2 spectra.
    a is a schematic drawn in code on the same canvas as b and c (one producer figure).
    b: disjoint/joint peptide-count ratio; c: disjoint/full ratio; 710 acquisitions.
d, e
    Same-spectrum benchmark against the matched metagenome search (former Fig. 4b,c of
    v28d): peptides of 5-50 residues per acquisition (d) and the fraction of
    metagenome-search peptide identities also found by FastaLake (e); 296 PREDICT and
    282 CapScan acquisitions. One producer figure.

Source data (``source_data/SFig06/``)
-------------------------------------
split_half_acquisitions_710.tsv
    panels b, c and the spectrum-split numbers of v28d lines 101, 123 and 433.
split_half_acquisitions_726.tsv
    the same experiment after the 347-acquisition CapScan re-run; sensitivity only.
same_spectrum_benchmark_0_predict.tsv, same_spectrum_benchmark_1_capscan.tsv
    panels d, e and the same-spectrum numbers of v28d line 123.

How the panels are drawn
------------------------
The drawing code of ``_render_abc`` and ``_render_de`` is the original producer code,
moved unchanged except for input and output paths (see PROVENANCE.md). The originals
were standalone scripts that started from matplotlib defaults. ``panel_abc`` and
``panel_de`` therefore run them in a fresh Python interpreter, so that the notebook's
own style settings cannot change a pixel. ``compose`` repeats the two PyMuPDF steps that
built the v28f composite page.
"""

import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

UNIT = Path(((os.environ.get('FASTALAKE_PAPER_REPO') or str(Path(__file__).resolve().parents[3])) + '/figures/SFig06'))  # len950 copy: repo unit, read-only
REPO = UNIT.parents[1]
SD = REPO / "source_data" / UNIT.name

# ---------------------------------------------------------------- choices
MIN_LENGTH = int(os.environ.get('LEN950_MIN', '5'))  # len950 copy: 5 (control) or 9; peptides of 5-50 residues; v28d line 267 and legend line 433
MAX_LENGTH = 50                # v28d line 267; applied upstream (record only, not re-applied here)
Q_CUTOFF = 0.01                # Sage peptide q <= 0.01; applied upstream (record only); v28d line 267
COHORT_ORDER = ["CAMPI_MOCK", "HSTOOL", "PREDICT", "CAPSCAN"]   # panel b,c x order (producer)
SPLIT_ROSTER_N = 710           # acquisitions in the v28f panels b,c (v28d lines 101, 433)
JITTER_SEED_BC = 20260927      # producer: default_rng(20260927 + panel index), b = +0, c = +1
JITTER_SEED_DE = 4             # producer: default_rng(4) for d and e
ABC_SAVE_DPI = 220             # producer savefig dpi for the a-c canvas
ABC_MATPLOTLIB = "3.7.0"       # matplotlib version that wrote the original a-c PDF (its /Creator)
BENCH_TEST_ARM = "FL_public"   # public-sequence FastaLake; v28d line 123
BENCH_REF_ARM = "MG"           # matched metagenome search; v28d line 123
BENCH_N = {"PREDICT": 296, "CAPSCAN": 282}   # matched acquisitions; v28d line 123
COMPOSE_Y_CUT_PT = 440.0       # composite keeps a-c above this y (pt) and puts d,e below
COMPOSITE_H_PT = 690.0         # composite page height (pt), 440 + 250
# The producer drawing code below is kept as it ran and uses literals. These asserts tie
# the constants above to those literals, so a changed constant stops the build instead of
# having no effect. To change a choice, change the literal in the producer code as well.
ROSTER = int(os.environ.get('LEN950_ROSTER', '710'))  # len950 copy: 710 (control) or 726 (v28f roster)
N_CAPSCAN = {710: 331, 726: 347}[ROSTER]
assert MIN_LENGTH in (5, 9) and (SPLIT_ROSTER_N, JITTER_SEED_BC, ABC_SAVE_DPI) == (710, 20260927, 220)
assert (JITTER_SEED_DE, COMPOSE_Y_CUT_PT, COMPOSITE_H_PT) == (4, 440.0, 690.0)
PLACED_W_MM = 180.34           # placed width of the figure in the manuscript (mm)
PNG_DPI_AT_PLACED = 300        # composite PNG resolution at the placed width


# ---------------------------------------------------------------- data
def load_split(roster=710):
    """Read the spectrum-split acquisition table.

    Parameters
    ----------
    roster : int
        710 (the v28f panels) or 726 (347-acquisition CapScan re-run, sensitivity).

    Returns
    -------
    pandas.DataFrame
        One row per acquisition and length view (``min_length`` 5 or 9).
    """
    return pd.read_csv(SD / f"split_half_acquisitions_{roster}.tsv", sep="\t")


def load_structure():
    """Read acquisitions, participants and specimens per cohort of the spectrum-split run.

    Returns
    -------
    pandas.DataFrame
        One row per cohort (710 roster, 5-50 residues), indexed by cohort; counted by
        ``prepare_source_data.py`` before the identifier columns were removed.
    """
    return pd.read_csv(SD / "split_half_structure.tsv", sep="\t").set_index("cohort")


def load_benchmark():
    """Read the same-spectrum benchmark rows used by panels d and e.

    Returns
    -------
    pandas.DataFrame
        Rows with ``min_length`` 5, ``test_arm`` FL_public and ``reference_arm`` MG, one
        per acquisition (296 PREDICT, 282 CapScan), exactly the producer's filter.
    """
    pa = pd.concat([pd.read_csv(SD / f, sep="\t") for f in
                    ("same_spectrum_benchmark_0_predict.tsv", "same_spectrum_benchmark_1_capscan.tsv")])
    pa = pa[pa.min_length.eq(MIN_LENGTH) & pa.test_arm.eq(BENCH_TEST_ARM) & pa.reference_arm.eq(BENCH_REF_ARM)]
    assert pa.groupby("cohort").size().to_dict() == {"CAPSCAN": BENCH_N["CAPSCAN"], "PREDICT": BENCH_N["PREDICT"]}
    return pa


def font_dir():
    """Return the Söhne font directory from ``FASTALAKE_FONT_DIR`` (or ``FL_STYLE_DIR``), or None.

    Returns
    -------
    pathlib.Path or None
        The directory if it holds ``soehne-buch.ttf``, else None.
    """
    d = os.environ.get("FASTALAKE_FONT_DIR") or os.environ.get("FL_STYLE_DIR")
    return Path(d) if d and (Path(d) / "soehne-buch.ttf").is_file() else None


# ---------------------------------------------------------------- panels a-c
def _render_abc(out_dir):
    """Draw panels a-c exactly as the original producer did (fresh interpreter only).

    Original: ``build_completed_figure.py`` (PROVENANCE.md). Changes: input path,
    output path and file stem, font path taken from ``FASTALAKE_FONT_DIR``, and the
    producer's FIGURE_CHECK.json (which recorded lab paths) is not written.

    Parameters
    ----------
    out_dir : pathlib.Path
        Output directory; receives ``panel_abc.png``, ``.pdf`` and ``.svg``.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
    sys.path.insert(0, str(REPO))
    from style._house import init_plotting
    if os.environ.get("FASTALAKE_PY_MPL370"):
        # the exact a-c check is only meaningful with the version that wrote the original
        assert matplotlib.__version__ == ABC_MATPLOTLIB, matplotlib.__version__
    fd = font_dir()
    init_plotting(font_path=str(fd / 'soehne-buch.ttf') if fd else None)
    plt.rcParams.update({'font.size':9,'axes.labelsize':10,'xtick.labelsize':9,'ytick.labelsize':9,'svg.fonttype':'none','pdf.fonttype':42})
    D = Path(out_dir)
    df=load_split(ROSTER);df=df[df.min_length==MIN_LENGTH];assert len(df)==ROSTER;NZERO=int((df[df.cohort=='CAPSCAN'].disjoint==0).sum())
    order=['CAMPI_MOCK','HSTOOL','PREDICT','CAPSCAN'];labels=['CAMPI\n15 runs','Healthy stool\n58 runs','PREDICT\n306 runs',f'CapScan\n{N_CAPSCAN} runs']
    fig=plt.figure(figsize=(9.3,6.4));gs=fig.add_gridspec(2,2,height_ratios=[1.5,3.25],hspace=.48,wspace=.31,left=.09,right=.98,bottom=.16,top=.94)
    a=fig.add_subplot(gs[0,:]);a.set(xlim=(0,10),ylim=(0,3.2));a.axis('off');a.text(-.58,3.2,'a',weight='bold',fontsize=13,va='top');a.text(0,3.15,'Completed experiment: select with different inputs, score the same spectra',fontsize=11,weight='bold',va='top')
    colors=['#6D7680','#4774A3','#B77A45']
    for y,title,db,col in zip([2.1,1.2,.3],['All predictions','Odd-scan predictions','Even-scan predictions'],['Full database','Joint database','Disjoint database'],colors):
     for x,width,text in [(0,2.45,title),(3.15,2.35,db)]:
      a.add_patch(FancyBboxPatch((x,y-.27),width,.59,boxstyle='round,pad=0.025,rounding_size=0.06',facecolor=col,edgecolor='none',alpha=.13));a.text(x+width/2,y+.025,text,ha='center',va='center',fontsize=9.5,color=col)
     a.add_patch(FancyArrowPatch((2.5,y+.025),(3.1,y+.025),arrowstyle='-|>',mutation_scale=10,color=col,linewidth=.9))
     a.add_patch(FancyArrowPatch((5.56,y+.025),(7.04,1.225),arrowstyle='-|>',mutation_scale=10,color=col,linewidth=.9))
    a.add_patch(FancyBboxPatch((7.1,.76),2.55,.94,boxstyle='round,pad=0.035,rounding_size=.06',facecolor='#EEF1F3',edgecolor='#B7BDC3',linewidth=.6));a.text(8.375,1.23,'Search identical\nodd-scan spectra',ha='center',va='center',fontsize=10)
    for panel,col,title,letter,color in [(0,'D_over_S','Disjoint / joint yield','b','#4774A3'),(1,'D_over_F','Disjoint / full yield','c','#B77A45')]:
     ax=fig.add_subplot(gs[1,panel]);ax.text(-.18,1.11,letter,transform=ax.transAxes,fontsize=13,weight='bold');ax.set_title(title,loc='left',fontsize=11,weight='bold',pad=14)
     rng=np.random.default_rng(20260927+panel)
     for i,co in enumerate(order):
      v=df[df.cohort==co][col].to_numpy()*100
      ax.scatter(i+rng.uniform(-.16,.16,len(v)),v,s=9,c=color,alpha=.27,edgecolors='none',rasterized=False,zorder=2)
      q1,med,q3=np.quantile(v,[.25,.5,.75]);ax.plot([i,i],[q1,q3],color='#30383E',lw=3.5,zorder=3);ax.plot([i-.12,i+.12],[med,med],color='white',lw=1.6,zorder=4)
      ax.text(i,105,f'{med:.1f}%',ha='center',fontsize=9,weight='bold',color=color)
     ax.axhline(100,color='#8B9399',lw=.8,ls=(0,(3,3)),zorder=1);ax.set(xlim=(-.5,3.5),ylim=(-3,112),xticks=range(4),xticklabels=labels,yticks=[0,25,50,75,100],ylabel='Peptide-count ratio (%)');ax.spines[['top','right']].set_visible(False);ax.tick_params(length=3)
     if panel==0:ax.annotate({1:'One zero-yield run\nretained',2:'Two zero-yield runs\nretained'}[NZERO],xy=(3,0),xytext=(2,19),fontsize=8,ha='center',arrowprops=dict(arrowstyle='-',lw=.7,color='#606970'))
    fig.text(.09,.055,f'Each point is an acquisition; bars show median and interquartile range. Accepted peptides: {MIN_LENGTH}–50 residues, native peptide q at most 0.01.\nCount ratios do not measure identical-peptide retention or whole-workflow FDR. Balanced reciprocal follow-up is separate.',fontsize=8,color='#555B60',va='bottom')
    for suffix in ['png','pdf','svg']:fig.savefig(D/f'panel_abc.{suffix}',dpi=220)
    plt.close(fig)


def abc_python():
    """Interpreter for panels a-c: ``FASTALAKE_PY_MPL370`` if set, else the running one.

    Returns
    -------
    str
        Path of the Python executable.

    The original a-c PDF was written by matplotlib 3.7.0 (Python 3.10.9, NumPy 1.23.5,
    pandas 1.5.3). Matplotlib 3.10 places some text differently (README, Known
    limitations), so the pixel check of a-c is exact only with that version.
    """
    return os.environ.get("FASTALAKE_PY_MPL370") or sys.executable


def _fresh(kind, out_dir):
    """Run ``_render_<kind>`` in a new interpreter (clean matplotlib state, as the producer had).

    Parameters
    ----------
    kind : str
        ``"abc"`` or ``"de"``.
    out_dir : pathlib.Path
        Output directory passed to the renderer.
    """
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONNOUSERSITE="1",
               PYTHONPATH=os.pathsep.join([str(REPO), str(UNIT), os.environ.get("PYTHONPATH", "")]))
    py = abc_python() if kind == "abc" else sys.executable
    r = subprocess.run([py, str(Path(__file__).resolve()), kind, str(out_dir)],
                       env=env, cwd=str(UNIT), capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(f"render {kind} failed:\n{r.stderr[-3000:]}")


def panel_abc(out_dir):
    """Render panels a-c (one canvas) into ``out_dir``.

    Parameters
    ----------
    out_dir : pathlib.Path
        Output directory; receives ``panel_abc.png`` (220 dpi), ``.pdf`` and ``.svg``.

    Returns
    -------
    pathlib.Path
        Path of ``panel_abc.pdf``.
    """
    _fresh("abc", out_dir)
    return Path(out_dir) / "panel_abc.pdf"


# ---------------------------------------------------------------- panels d,e
def _render_de(out_dir):
    """Draw panels d,e exactly as the original producer did (fresh interpreter only).

    Original: the "Supplementary Fig. 6: old Fig. 4b,c as panels d,e" block of
    ``fig4bc_crossfit/build/build.py`` (PROVENANCE.md). The rcParams that block inherited
    from the Fig. 4 block run before it in the same script are set first, with the same
    values. Changes: input and output paths, font path from ``FASTALAKE_FONT_DIR``, and
    the page width read from the regenerated a-c PDF (the producer read it from the
    composite page that carried the same a-c PDF).

    Parameters
    ----------
    out_dir : pathlib.Path
        Output directory holding ``panel_abc.pdf``; receives ``panel_de.pdf`` and
        ``panel_de_plotted_values.tsv``.
    """
    import matplotlib; matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    import fitz
    sys.path.insert(0, str(REPO))
    from style import fl_style as fl
    D = Path(out_dir)
    fd = font_dir()
    if fd:
        FONT = fd / 'soehne-buch.ttf'
        font_manager.fontManager.addfont(str(FONT)); BUCH = font_manager.FontProperties(fname=str(FONT)).get_name()
    else:
        BUCH = 'sans-serif'
    LG = '#cfcfcf'
    # rcParams set by the Fig. 4 block of build.py, which ran first in the producer process
    plt.rcParams.update({'font.family': BUCH, 'font.size': 6.0, 'axes.labelsize': 6.0, 'xtick.labelsize': 6.0, 'ytick.labelsize': 6.0,
                         'axes.linewidth': .46, 'xtick.major.width': .46, 'ytick.major.width': .46, 'xtick.major.size': 3.2, 'ytick.major.size': 3.2,
                         'axes.spines.top': False, 'axes.spines.right': False, 'pdf.fonttype': 42, 'savefig.bbox': None,
                         'text.color': fl.INK, 'axes.labelcolor': fl.INK, 'axes.edgecolor': fl.INK, 'xtick.color': fl.INK, 'ytick.color': fl.INK,
                         'axes.unicode_minus': True})
    pa = load_benchmark()
    # =============================== Supplementary Fig. 6: old Fig. 4b,c as panels d,e ===============================
    W6 = fitz.open(D / 'panel_abc.pdf')[0].rect.width; H6 = COMPOSITE_H_PT
    assert round(W6, 1) == 669.6
    Y6 = 440.0; HD = H6 - Y6
    INK6 = '#4e4f51'
    plt.rcParams.update({'font.size': 9.0, 'axes.labelsize': 9.0, 'xtick.labelsize': 9.0, 'ytick.labelsize': 9.0,
                         'axes.linewidth': .6, 'xtick.major.width': .6, 'ytick.major.width': .6, 'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
                         'text.color': INK6, 'axes.labelcolor': INK6, 'axes.edgecolor': INK6, 'xtick.color': INK6, 'ytick.color': INK6})
    fig = plt.figure(figsize=(W6 / 72, HD / 72))
    def ax6(x0, y0, x1, y1):
        """Add axes at a box given in pt from the top-left of the d,e canvas (producer code, unchanged)."""
        return fig.add_axes([x0 / W6, 1 - y1 / HD, (x1 - x0) / W6, (y1 - y0) / HD])
    dY0, dY1 = 40.0, 190.0
    axD = ax6(80.0, dY0, 400.0, dY1); axE = ax6(500.0, dY0, 650.0, dY1)
    rng = np.random.default_rng(4); xs_ = {'PREDICT': (0, 1), 'CAPSCAN': (2.6, 3.6)}; pv6 = []
    for coh, (x0, x1) in xs_.items():
        d = pa[pa.cohort.eq(coh)]; j0 = rng.uniform(-.14, .14, len(d)); j1 = rng.uniform(-.14, .14, len(d))
        for a_, b_, ja, jb in zip(d.n_reference / 1e3, d.n_test / 1e3, j0, j1):
            axD.plot([x0 + ja, x1 + jb], [a_, b_], color=LG, lw=.35, alpha=.6, zorder=1)
        axD.scatter(x0 + j0, d.n_reference / 1e3, s=5, color=fl.METAG, alpha=.7, lw=0, zorder=2)
        axD.scatter(x1 + j1, d.n_test / 1e3, s=5, color=fl.FL, alpha=.7, lw=0, zorder=2)
        for xx, v in ((x0, d.n_reference.median()), (x1, d.n_test.median())):
            axD.plot([xx - .3, xx + .3], [v / 1e3] * 2, color='k', lw=1.4, zorder=3)
            axD.text(xx + .34, v / 1e3, f'{v / 1e3:.1f}k', ha='left', va='center', fontsize=8.0, zorder=4)
        assert max(d.n_reference.max(), d.n_test.max()) / 1e3 < 55
        for _, r in d.iterrows():
            pv6.append(dict(panel='SuppFig6d,e', cohort=coh, sample=r['sample'], participant=r.subject_id, **{f'mg_peptides_{MIN_LENGTH}_50':int(r.n_reference)},
                            **{f'fl_peptides_{MIN_LENGTH}_50':int(r.n_test), f'recovery_pct_{MIN_LENGTH}_50':100 * r.reference_recovery},
                            source='Figure4/benchmark_0_per_acquisition.tsv + benchmark_1_per_acquisition.tsv (min_length '+str(MIN_LENGTH)+', FL_public vs MG)'))
    axD.set_xticks([0, 1, 2.6, 3.6], ['MG', 'FL', 'MG', 'FL']); axD.set_ylabel(f'Peptides per acquisition\n(thousands; same spectra,\n{MIN_LENGTH}–50 residues)', labelpad=3)
    axD.set_ylim(0, 55); axD.set_xlim(-.6, 4.2)
    axD.text(.5, -.13, 'PREDICT (n = 296)', ha='center', va='top', transform=axD.get_xaxis_transform())
    axD.text(3.1, -.13, 'CapScan (n = 282)', ha='center', va='top', transform=axD.get_xaxis_transform())
    for i, coh in enumerate(['PREDICT', 'CAPSCAN']):
        d = pa[pa.cohort.eq(coh)]; j = rng.uniform(-.18, .18, len(d))
        axE.scatter(i + j, 100 * d.reference_recovery, s=5, color=fl.FL, alpha=.7, lw=0)
        m = 100 * d.reference_recovery.median(); axE.plot([i - .3, i + .3], [m, m], color='k', lw=1.4, zorder=3)
        axE.text(i + .34, m, f'{m:.1f}%', ha='left', va='center', fontsize=8.0, zorder=4)
    axE.set_xticks([0, 1], ['PREDICT', 'CapScan']); axE.set_ylabel(f'MG peptides found by FL (%)\nsame spectra, {MIN_LENGTH}–50 residues', labelpad=3)
    axE.set_ylim(0, 105); axE.set_xlim(-.6, 1.9)
    fig.text(13.8 / W6, 1 - 26 / HD, 'd', fontsize=13, color=INK6, va='baseline')
    fig.text(440.0 / W6, 1 - 26 / HD, 'e', fontsize=13, color=INK6, va='baseline')
    SF6P = D / 'panel_de.pdf'; fig.savefig(SF6P); plt.close(fig)
    pd.DataFrame(pv6).to_csv(D / 'panel_de_plotted_values.tsv', sep='\t', index=False)


def panel_de(out_dir):
    """Render panels d,e (one canvas) into ``out_dir``; needs ``panel_abc.pdf`` there.

    Parameters
    ----------
    out_dir : pathlib.Path
        Output directory; receives ``panel_de.pdf`` and ``panel_de_plotted_values.tsv``.

    Returns
    -------
    pathlib.Path
        Path of ``panel_de.pdf``.
    """
    _fresh("de", out_dir)
    return Path(out_dir) / "panel_de.pdf"


def rasterise(pdf, png, dpi=300):
    """Save page 0 of ``pdf`` as a PNG at ``dpi`` (for display in the notebook).

    Parameters
    ----------
    pdf : pathlib.Path
        PDF to rasterise (page 0 only).
    png : pathlib.Path
        Output PNG path.
    dpi : int
        Resolution of the PNG.

    Returns
    -------
    pathlib.Path
        Path of the PNG.
    """
    import fitz
    pix = fitz.open(pdf)[0].get_pixmap(dpi=dpi, alpha=False)
    pix.save(png)
    return Path(png)


# ---------------------------------------------------------------- composite
def compose(out_dir):
    """Assemble the v28f composite from ``panel_abc.pdf`` and ``panel_de.pdf``.

    Repeats the producer steps (ASSEMBLY.md): (1) ``sf6d.py`` placed the a-c PDF at the
    top of a 690-pt page and checked that it has no ink below y = 440 pt; (2)
    ``build.py`` kept that page above y = 440 pt, placed the d,e PDF below it, and
    rasterised the page so that the placed width of 180.34 mm is 300 dpi. In step 1 the
    producer also drew the earlier panel d below 440 pt; it is cut away in step 2, so it
    is not drawn here.

    Parameters
    ----------
    out_dir : pathlib.Path
        Directory holding ``panel_abc.pdf`` and ``panel_de.pdf``; receives
        ``step1_abc_page.pdf``, ``SFig06_v28f.pdf`` and ``SFig06_v28f.png``.

    Returns
    -------
    dict
        Paths of ``SFig06_v28f.pdf`` and ``SFig06_v28f.png`` and the PNG size in pixels.
    """
    import fitz
    from PIL import Image  # noqa: F401  (pixmap save uses PyMuPDF)
    D = Path(out_dir)
    # step 1 (sf6d.py, compose)
    src = fitz.open(D / 'panel_abc.pdf'); W = src[0].rect.width; H0 = src[0].rect.height
    Y0 = COMPOSE_Y_CUT_PT; HD = COMPOSITE_H_PT - Y0
    blank = src[0].get_pixmap(dpi=150, alpha=False, clip=fitz.Rect(0, Y0, W, H0))
    assert np.frombuffer(blank.samples, np.uint8).min() == 255, 'ink in the overlap strip'
    doc = fitz.open(); page = doc.new_page(width=W, height=Y0 + HD)
    page.show_pdf_page(fitz.Rect(0, 0, W, H0), src, 0)
    step1 = D / 'step1_abc_page.pdf'; doc.save(step1, garbage=4, deflate=True); doc.close()
    # step 2 (build.py, Supplementary Fig. 6 block)
    s6 = fitz.open(step1); W6, H6 = s6[0].rect.width, s6[0].rect.height
    assert (round(W6, 1), round(H6, 1)) == (669.6, 690.0)
    Y6 = 440.0
    doc = fitz.open(); page = doc.new_page(width=W6, height=H6)
    page.show_pdf_page(fitz.Rect(0, 0, W6, Y6), s6, 0, clip=fitz.Rect(0, 0, W6, Y6))
    page.show_pdf_page(fitz.Rect(0, Y6, W6, H6), fitz.open(D / 'panel_de.pdf'), 0)
    S6 = D / 'SFig06_v28f.pdf'; doc.save(S6, garbage=4, deflate=True); doc.close()
    pg6 = fitz.open(S6)[0]
    k6 = PLACED_W_MM * 72 / 25.4 / W6
    dpi6 = PNG_DPI_AT_PLACED * k6
    pix6 = pg6.get_pixmap(matrix=fitz.Matrix(dpi6 / 72, dpi6 / 72), alpha=False); pix6.set_dpi(300, 300)
    pix6.save(D / 'SFig06_v28f.png')
    return {"pdf": S6, "png": D / 'SFig06_v28f.png', "png_px": [pix6.width, pix6.height]}


# ---------------------------------------------------------------- statistics
def stat_split_summary(df):
    """Cohort summaries of the spectrum-split experiment at 5-50 residues.

    Statistical unit: acquisition (n = 15 CAMPI, 58 healthy stool, 306 PREDICT and 331
    CapScan for the 710 roster; 347 CapScan for the 726 roster). Descriptive only: no
    test, no interval beyond the interquartile range. Medians and quartiles are
    ``pandas.Series.median`` and ``quantile`` (linear interpolation), as in the
    producer's ``review.py``; the plotted bars use ``numpy.quantile`` (same definition).

    Parameters
    ----------
    df : pandas.DataFrame
        Output of ``load_split``.

    Returns
    -------
    pandas.DataFrame
        One row per cohort with n, zero-yield counts, number of acquisitions with fewer
        disjoint than joint peptides, and median/q25/q75 of D/S, D/F and of the fraction
        of joint-search peptide identities also accepted with disjoint selection, all as
        fractions.
    """
    g5 = df[df.min_length == MIN_LENGTH]
    rows = []
    for co in COHORT_ORDER:
        g = g5[g5.cohort == co]
        r = dict(cohort=co, n=len(g), zero_disjoint=int((g.disjoint == 0).sum()),
                 D_lower_S=int((g.disjoint < g.self).sum()))
        for col in ("D_over_S", "D_over_F", "fraction_self_identities_retained"):
            v = g[col]
            r[f"median_{col}"] = float(v.median()); r[f"q25_{col}"] = float(v.quantile(.25)); r[f"q75_{col}"] = float(v.quantile(.75))
        rows.append(r)
    return pd.DataFrame(rows).set_index("cohort")


def stat_split_text(df):
    """Numbers printed for the spectrum-split experiment (v28d lines 101, 123, 433).

    Statistical unit: acquisition; n = 710 (15 + 58 + 306 + 331). Descriptive medians.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of ``load_split`` (710 or 726 roster).

    Returns
    -------
    dict
        ``n_total``; ``median_DS_pct`` per cohort (100 x median D/S); ``loss_pct`` per
        cohort (100 x (1 - median D/S)); ``retention_pct`` per cohort (100 x median
        fraction of joint identities retained) and its minimum and maximum over cohorts;
        ``n_capscan_zero_disjoint``; ``all_fewer`` (True if every acquisition has fewer
        disjoint than joint peptides).
    """
    s = stat_split_summary(df)
    ret = 100 * s["median_fraction_self_identities_retained"]
    return dict(n_total=int(s.n.sum()),
                n_by_cohort=s.n.to_dict(),
                median_DS_pct=(100 * s.median_D_over_S).to_dict(),
                median_DF_pct=(100 * s.median_D_over_F).to_dict(),
                loss_pct=(100 * (1 - s.median_D_over_S)).to_dict(),
                retention_pct=ret.to_dict(), retention_min_pct=float(ret.min()), retention_max_pct=float(ret.max()),
                n_capscan_zero_disjoint=int(s.loc["CAPSCAN", "zero_disjoint"]),
                all_fewer=bool((s.D_lower_S == s.n).all()))


def stat_bench(pa):
    """Same-spectrum benchmark summaries for panels d,e and v28d line 123.

    Statistical unit: acquisition (PREDICT n = 296, CapScan n = 282). Descriptive
    medians; no test and no interval. Gain is ``100 * (median(count_ratio) - 1)`` with
    ``count_ratio`` = FastaLake / metagenome peptides in the same acquisition, as in the
    producer. Recovery is the median over acquisitions of |FL and MG| / |MG|.

    Parameters
    ----------
    pa : pandas.DataFrame
        Output of ``load_benchmark`` (one row per acquisition, 5-50 residues,
        FL_public vs MG).

    Returns
    -------
    dict
        Per cohort: n, median MG and FL peptides, median recovery (%), number of
        acquisitions with more FL than MG peptides, median gain (%).
    """
    out = {}
    for coh in ("PREDICT", "CAPSCAN"):
        d = pa[pa.cohort.eq(coh)]
        out[coh] = dict(n=len(d), median_MG=float(d.n_reference.median()), median_FL=float(d.n_test.median()),
                        median_recovery_pct=float(100 * d.reference_recovery.median()),
                        n_FL_gt_MG=int((d.n_test > d.n_reference).sum()),
                        median_gain_pct=float(100 * (d.count_ratio.median() - 1)))
    return out


def stat_726_vs_710(df710, df726):
    """Compare the 710 and 726 spectrum-split tables (sensitivity, not plotted).

    Statistical unit: acquisition. For the 331 CapScan acquisitions present in both,
    counts how many differ in the joint (``self``), disjoint and full peptide counts at
    5-50 residues, and the largest absolute difference.

    Parameters
    ----------
    df710, df726 : pandas.DataFrame
        Outputs of ``load_split(710)`` and ``load_split(726)``.

    Returns
    -------
    dict
        ``text_726`` (output of ``stat_split_text`` on the 726 table) and per-column
        difference counts for the shared CapScan acquisitions.
    """
    a = df710[(df710.cohort == "CAPSCAN") & (df710.min_length == MIN_LENGTH)].set_index("sample")
    b = df726[(df726.cohort == "CAPSCAN") & (df726.min_length == MIN_LENGTH)].set_index("sample")
    common = a.index.intersection(b.index)
    diff = {}
    for col in ("self", "disjoint", "full"):
        dd = (a.loc[common, col] - b.loc[common, col]).abs()
        diff[col] = dict(n_differ=int((dd > 0).sum()), max_abs_diff=int(dd.max()))
    return dict(text_726=stat_split_text(df726), n_capscan_shared=len(common),
                n_capscan_only_726=int(len(b.index.difference(a.index))), capscan_count_differences=diff)


NUMBER_WORDS = {0: "zero", 1: "one", 2: "two", 3: "three"}


def as_word(n):
    """Return the English word the manuscript uses for a small count (``1`` -> ``'one'``).

    Parameters
    ----------
    n : int
        Count from 0 to 3.

    Returns
    -------
    str
        The word.
    """
    return NUMBER_WORDS[int(n)]


if __name__ == "__main__":
    kind, out = sys.argv[1], Path(sys.argv[2])
    out.mkdir(parents=True, exist_ok=True)
    {"abc": _render_abc, "de": _render_de}[kind](out)
