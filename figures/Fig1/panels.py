"""Figure 1: plotting and statistics functions.

Panels
------
a  schematic; shipped as ``artwork/panel_a.pdf`` (panel-a crop of the v30 composite
   ``artwork/Figure1_v30.pdf``). Its FastaLake column is rebuilt by
   ``artwork/build_fig1a_funnel.py`` (:func:`build_panel_a_funnel`); the rest is
   artwork (see ASSEMBLY.md and ARTWORK_NUMBERS.md)
b  selected FASTA sizes, 726 acquisitions            ``acquisition_sizes_726.tsv``, ``human_reference.json``
c  exact reservoir match by peptide length            ``length_by_arm.tsv``, ``five_seed_summary.tsv``
d  proteins supported by exactly one production query ``acquisition_support_summary_726.tsv``
e  donor-selected entrapment, 58 acquisitions         ``donor_seed_means.tsv`` (and ``donor_recount.tsv``)

Text statistics of Results lines 51-65 also use ``scramble_end_to_end.tsv``,
``selection_grid_acquisitions.tsv``, ``stage1_stats.json`` and
``worked_acquisition_prediction_support.json``. All files are in ``source_data/Fig1/``.

Origin of the plot code
-----------------------
Panels b, d and e are drawn by one 180 x 52 mm four-axis figure, moved unchanged from the
original producer ``build_panels_bd_n15.py`` (mode ``n15``; it is itself a copy of
``build_panels_bcd_v10.py`` with a data switch). Panels b and d are saved as crops of that
figure exactly as the producer did. Panel e was placed in the published figure from the
v10 run of the same plot code and is checked against that render. Panel c is moved
unchanged from ``build_FAR12.py``, which replaced the v10 panel c. Only file paths were
changed; see PROVENANCE.md.
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd

# ---- choices (each value is used by the code below) --------------------------------
Q_CUTOFF = 0.01                    # peptide q, Results v28d line 53 and Methods line 727
ENTRAPMENT_TIER = "T1_within_group"  # donor-selected transfer within design, Methods line 727
ENTRAPMENT_LENGTHS = (5, 9)        # 5-50 and 9-50 residues, legend line 71, Methods line 729
ENTRAPMENT_SEEDS = (20260910, 20260911, 20260912)  # three seeds, Methods line 729 (values not printed)
SCRAMBLE_SEEDS = (20260901, 20260913, 20260914, 20260915, 20260916)  # Methods line 665
WORKED_SEED = 20260901             # seed propagated through the full workflow, Methods line 665
WORKED_ACQUISITION = "20181231_QX7_PhGe_LC95_SA_WorkflowReplicates_2"  # representative acquisition
TOP_FRACTION_ARM = "tp0.30"        # top-30% lookup fraction, Methods line 557
SUPPORT_POPULATION = "all_selected"  # all selected proteins, Methods line 677
COHORTS = ("CAMPI", "HSTOOL", "PREDICT", "CapScan")  # panel order of b and d
ROSTER = {"CAMPI": 15, "HSTOOL": 58, "PREDICT": 306, "CapScan": 347}  # 726 acquisitions, legend line 71
JITTER_SEED_BASE = 610             # point jitter: default_rng(610 + column index), from the producer
SCRAMBLED_MIN_LENGTH_NO_MATCH = 13  # ">= 13 residues", Results line 59
DPI = 300                          # producer PNG resolution and pixel-check resolution

_ROW_CACHE = {}
_BASE_RC = {}

# ---- the producers' style module -------------------------------------------------------
# The original producers imported the lab copy of fl_style.py (sha256
# 58b64c343f095cbace09c3e9766e3e78060c1c8cf522b075c4047300322164c3; PROVENANCE.md).
# The repository's style/fl_style.py is a later revision of that module that sets 1 pt
# lines, bold panel letters and savefig.bbox None. To regenerate the panels as they were
# produced, these functions apply the producer-era values; colours, fonts and every other
# setting come from style/fl_style.py, where they are identical. Values copied verbatim
# from the producer module's init(), style() and panel().
PRODUCER_RC = {
    "lines.linewidth": 0.6, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.minor.width": 0.6, "ytick.minor.width": 0.6,
    "grid.linewidth": 0.8, "patch.linewidth": 1.0,   # matplotlib defaults; the producer module did not set them
    "savefig.bbox": "tight", "savefig.pad_inches": 0.02,
}


def _style(ax):
    """``style()`` of the producers' fl_style module (0.6 pt ticks)."""
    for sp in ('top', 'right'):
        ax.spines[sp].set_visible(False)
    ax.tick_params(length=2.6, width=0.6, labelsize=7)
    ax.xaxis.label.set_size(7)
    ax.yaxis.label.set_size(7)
    return ax


def _panel(ax, letter, fontsize=8):
    """``panel()`` of the producers' fl_style module (regular-weight panel letter)."""
    from style import fl_style as fl
    ax.set_title(str(letter).strip(), loc='left', fontsize=fontsize, color=fl.INK,
                 fontweight='normal', pad=4)


# ---- loading ------------------------------------------------------------------------
def load(sd):
    """Read every Figure 1 source table.

    Parameters
    ----------
    sd : pathlib.Path
        ``source_data/Fig1``.

    Returns
    -------
    dict
        Keys ``sizes``, ``human``, ``lba``, ``fss``, ``sup``, ``ent``, ``recount``,
        ``e2e``, ``grid``, ``stage1``, ``worked``.
    """
    import matplotlib as mpl
    if not _BASE_RC:  # rcParams as the producers' fl.init() left them; each panel starts from these
        _BASE_RC.update(mpl.rcParams.copy())
        _BASE_RC.update(PRODUCER_RC)
    sd = Path(sd)
    rd = lambda n: pd.read_csv(sd / n, sep="\t")  # noqa: E731
    return dict(
        sizes=rd("acquisition_sizes_726.tsv"),
        human=json.loads((sd / "human_reference.json").read_text()),
        lba=rd("length_by_arm.tsv"),
        fss=rd("five_seed_summary.tsv"),
        sup=rd("acquisition_support_summary_726.tsv"),
        ent=rd("donor_seed_means.tsv"),
        recount=rd("donor_recount.tsv"),
        e2e=rd("scramble_end_to_end.tsv"),
        grid=rd("selection_grid_acquisitions.tsv"),
        stage1=json.loads((sd / "stage1_stats.json").read_text()),
        worked=json.loads((sd / "worked_acquisition_prediction_support.json").read_text()),
    )


# ---- panel a: artwork ------------------------------------------------------------------
def panel_a(unit, out, dpi=150):
    """Panel a: rasterise the schematic ``artwork/panel_a.pdf`` for display.

    The schematic is not regenerated. Its typed numbers are recomputed by
    :func:`artwork_numbers` and listed in ARTWORK_NUMBERS.md.

    Parameters
    ----------
    unit : pathlib.Path
        ``figures/Fig1``.
    out : pathlib.Path
        Directory for ``panel_a.png``.
    dpi : int
        Display resolution.

    Returns
    -------
    pathlib.Path
        The PNG written.
    """
    import fitz
    dest = Path(out) / "panel_a.png"
    with fitz.open(Path(unit) / "artwork" / "panel_a.pdf") as doc:
        doc[0].get_pixmap(dpi=dpi, alpha=False).save(dest)
    return dest


def build_panel_a_funnel(unit, out, python=None):
    """Panel a: rebuild the v33 FastaLake column with ``artwork/build_fig1a_funnel.py``.

    Runs the ported producer in a subprocess (it is a script, not a module) with
    ``FASTALAKE_OUT_DIR = out / "funnel"``. Needs ``FASTALAKE_FONT_DIR`` (Soehne). The
    script asserts its own layout checks and 0 changed pixels outside the column.

    Parameters
    ----------
    unit : pathlib.Path
        ``figures/Fig1``.
    out : pathlib.Path
        Build directory; the result is ``out/funnel/Figure1.pdf`` and ``Figure1.png``.
    python : str, optional
        Interpreter; defaults to the running one.

    Returns
    -------
    pathlib.Path
        The rebuilt ``Figure1.pdf``.
    """
    import os
    import subprocess
    import sys
    dest = Path(out) / "funnel"
    env = dict(os.environ, FASTALAKE_OUT_DIR=str(dest), MPLBACKEND="Agg")
    subprocess.run([python or sys.executable, str(Path(unit) / "artwork" / "build_fig1a_funnel.py")],
                   env=env, check=True, stdout=subprocess.DEVNULL)
    return dest / "Figure1.pdf"


def artwork_numbers(d):
    """Recompute every number typed in panel a, as printed there.

    Statistical unit: one acquisition (``WORKED_ACQUISITION``) except the reservoir and
    the study-level lookup (58 healthy-stool acquisitions pooled).

    Returns
    -------
    dict
        Typed text -> value formatted as in the artwork.
    """
    a, w, h = stat_a_reservoir(d), stat_a_worked(d), d["human"]["entries"]
    return {
        "20,431 reference proteins": f"{h:,}",
        "233 M proteins (reservoir, and table Real and Scrambled)": f"{a['reservoir'] / 1e6:.0f} M",
        "18 M study-supported proteins (table, Real)": f"{a['study_supported'] / 1e6:.0f} M",
        "5,841 study-supported proteins (table, Scrambled; v30)": f"{w['scr_candidates']:,}",
        "24,560 unique peptides": f"{w['queries']:,}",
        "36,588 unique peptides (scrambled)": f"{w['scr_queries']:,}",
        "6.67 M candidates, this acquisition": f"{w['candidates'] / 1e6:.2f} M",
        "5,841 candidates (scrambled)": f"{w['scr_candidates']:,}",
        "12,606 sample FASTA after parsimony": f"{w['fasta']:,}",
        "525 sample FASTA (scrambled)": f"{w['scr_fasta']:,}",
        "29,663 peptides": f"{w['peptides']:,}",
        "0 peptides (scrambled)": f"{w['scr_peptides']:,}",
    }


def show(path):
    """Return an IPython image of a rendered panel for display in the notebook."""
    from IPython.display import Image
    return Image(filename=str(path))


# ---- panels b, d, e: the producer's four-axis figure ---------------------------------
def _row(d):
    """Draw the producer's 180 x 52 mm four-axis figure (b, v10 c, d, e) once and cache it.

    The body is ``build_panels_bd_n15.py`` (mode n15) moved unchanged; only the data
    arguments replace its file reads. The v10 axis c is drawn because the layout of b, d
    and e depends on it; it is not the published panel c.

    Returns
    -------
    tuple
        ``(fig, axs, MED)`` where MED holds the panel-e medians and quartiles (percent).
    """
    if "row" in _ROW_CACHE:
        return _ROW_CACHE["row"]
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    from style import fl_style as fl

    sizes, human, lba, fss = d["sizes"], d["human"], d["lba"], d["fss"]
    sup = d["sup"][d["sup"].population.eq(SUPPORT_POPULATION)]
    ent = d["ent"][d["ent"].tier.eq(ENTRAPMENT_TIER) & d["ent"].q_threshold.eq(Q_CUTOFF)]
    plt.rcParams.update(_BASE_RC)
    plt.rcParams.update({'axes.labelcolor':fl.INK,'text.color':fl.INK,'xtick.color':fl.INK,'ytick.color':fl.INK,'font.size':6.5,'axes.labelsize':6.5,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'legend.fontsize':6})
    assert sizes.groupby('cohort').size().to_dict()=={'CAMPI':15,'HSTOOL':58,'PREDICT':306,'CapScan':347} and sizes.entries.lt(human['entries']).all()
    W,H=180,52
    fig,axs=plt.subplots(1,4,figsize=fl.mm(W,H),layout='constrained');fig.set_constrained_layout_pads(w_pad=.04,h_pad=.03,wspace=.09,hspace=.0)
    def points(ax,values,labels,color=fl.FL):
        for i,v in enumerate(values):
            v=np.asarray(v,float);v=v[np.isfinite(v)];j=np.random.default_rng(JITTER_SEED_BASE+i).uniform(-.14,.14,len(v))
            ax.scatter(i+j,v,s=7 if len(v)<80 else 3,color=color,alpha=.6,lw=0);ax.plot([i-.22,i+.22],[np.median(v)]*2,color=fl.INK,lw=1.1,zorder=5)
        ax.set_xticks(range(len(labels)),labels);ax.set_xlim(-.55,len(labels)-.45);ax.tick_params(axis='x',length=0,labelsize=6.0)
        if len(labels)>2:
            for lab in ax.get_xticklabels():lab.set_rotation(35);lab.set_ha('right');lab.set_rotation_mode('anchor')
    for ax,l in zip(axs,'bcde'):_style(ax);ax.tick_params(labelsize=6.5);ax.xaxis.label.set_size(6.5);ax.yaxis.label.set_size(6.5);_panel(ax,l)
    co=list(COHORTS);NAME={'HSTOOL':'Healthy stool','CAPSCAN':'CapScan'};disp=lambda c:NAME.get(c,c)
    ax=axs[0];points(ax,[sizes.loc[sizes.cohort.eq(c),'entries'] for c in co],[f'{disp(c)} (n = {int(sizes.cohort.eq(c).sum())})' for c in co])
    ax.axhline(human['entries'],color=fl.GREY,ls='--',lw=.8);ax.text(3.45,human['entries']+450,f"Human proteome {human['entries']:,}",ha='right',va='bottom',fontsize=5.5,color=fl.INK)
    ax.set_ylim(0,24500);ax.set_ylabel('Selected FASTA entries');ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x,p:f'{x:,.0f}'))
    ax=axs[1];real=lba[lba.arm.eq('REAL')].sort_values('length');scr=lba[~lba.arm.eq('REAL')]
    for seed,g in scr.groupby('seed'):ax.plot(g.length,100*g.fraction_matched,color=fl.ORANGE,lw=.6,alpha=.55)
    ax.plot(fss.length,100*fss.scrambled_median,'o-',color=fl.ORANGE,ms=2.6,lw=.9,label='Scrambled, five seeds')
    ax.plot(real.length,100*real.fraction_matched,'o-',color=fl.FL,ms=2.6,lw=.9,label='Authentic de novo')
    ax.set_yscale('symlog',linthresh=.01,linscale=.6);ax.set_yticks([0,.01,.1,1,10,100]);ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x,p:f'{x:g}'))
    ax.set_ylabel('Distinct peptides with an\nexact reservoir match (%)');ax.set_xlabel('Peptide length (residues)');ax.legend(frameon=False,loc='center right',fontsize=6,handlelength=1.6)
    ax=axs[2];key={c.upper():c for c in sup.cohort.unique()};cc=[key[c.upper()] for c in co]
    points(ax,[100*sup.loc[sup.cohort.eq(c),'one_top_prediction_fraction'] for c in cc],[f'{disp(c)} (n = {int(sup.cohort.eq(c).sum())})' for c in cc])
    ax.set_ylim(0,100);ax.set_ylabel('Selected entries supported by\nexactly one prediction (%)')
    ax=axs[3];e5=100*ent.loc[ent.min_length.eq(5),'paired_estimate'];e9=100*ent.loc[ent.min_length.eq(9),'paired_estimate'];assert len(e5)==len(e9)==58
    points(ax,[e5,e9],[f'5–50\nresidues',f'9–50\nresidues'])
    ax.axhline(1.0,color=fl.GREY,ls='--',lw=.8);ax.text(1.45,1.04,'nominal 1%',ha='right',va='bottom',fontsize=6,color=fl.INK)
    ax.set_ylim(0,1.35);ax.set_ylabel('Donor-selected entrapment\nestimate (%), n = 58')
    MED=dict(e5_median=float(np.median(e5)),e9_median=float(np.median(e9)),e5_q25=float(np.percentile(e5,25)),e5_q75=float(np.percentile(e5,75)),e9_q25=float(np.percentile(e9,25)),e9_q75=float(np.percentile(e9,75)))
    labs_b=[t.get_text() for t in axs[0].get_xticklabels()];labs_d=[t.get_text() for t in axs[2].get_xticklabels()];assert 'Healthy stool (n = 58)' in labs_b and 'Healthy stool (n = 58)' in labs_d and not any('HSTOOL' in l for l in labs_b+labs_d),(labs_b,labs_d)
    _ROW_CACHE["row"] = (fig, axs, MED)
    return _ROW_CACHE["row"]


def _crop(fig, axs, i):
    """Return the producer's crop box (inches) around axis ``i`` (b = 0, d = 2, e = 3)."""
    from matplotlib.transforms import Bbox
    fig.canvas.draw();r=fig.canvas.get_renderer();tb=[a.get_tightbbox(r).transformed(fig.dpi_scale_trans.inverted()) for a in axs];FW,FH=fig.get_size_inches()
    cuts=[0]+[(tb[k].x1+tb[k+1].x0)/2 for k in range(3)]+[FW]
    PAD=.04
    return Bbox([[max(cuts[i],tb[i].x0-PAD),0],[min(cuts[i+1],tb[i].x1+PAD),FH]])


def _save_crop(d, out, letter, i):
    fig, axs, _ = _row(d)
    bb = _crop(fig, axs, i)
    for ext in ('pdf', 'png'):
        fig.savefig(Path(out) / f'panel_{letter}.{ext}', dpi=DPI, bbox_inches=bb)
    return fig


def panel_b(d, out):
    """Panel b: entries of each acquisition-specific FASTA, 726 acquisitions, by cohort.

    Parameters
    ----------
    d : dict
        Output of :func:`load`.
    out : pathlib.Path
        Directory for ``panel_b.pdf`` and ``panel_b.png`` (producer crop of the row).

    Returns
    -------
    matplotlib.figure.Figure
        The four-axis row the crop was taken from.
    """
    return _save_crop(d, out, 'b', 0)


def panel_d(d, out):
    """Panel d: fraction of selected proteins with exactly one production query, per acquisition.

    Parameters
    ----------
    d : dict
        Output of :func:`load`.
    out : pathlib.Path
        Directory for ``panel_d.pdf`` and ``panel_d.png``.

    Returns
    -------
    matplotlib.figure.Figure
    """
    return _save_crop(d, out, 'd', 2)


def panel_e(d, out):
    """Panel e: donor-selected entrapment estimates, 58 acquisitions, 5-50 and 9-50 residues.

    The published panel e is the v10 render of the same plot code, which saved only the
    whole row. This function saves the regenerated row (``panels_bcde_row.pdf/.png``) and the
    panel-e crop (``panel_e.png``) taken from the row raster at 300 dpi with the producer's
    crop rule (see :func:`crop_panel_e`).

    Returns
    -------
    matplotlib.figure.Figure
    """
    fig, axs, _ = _row(d)
    out = Path(out)
    for ext in ('pdf', 'png'):  # as the producer's fl.save(): tight bbox, white face
        fig.savefig(out / f'panels_bcde_row.{ext}', dpi=DPI, bbox_inches='tight', facecolor='white', edgecolor='none')
    crop_panel_e(out / 'panels_bcde_row.pdf', out / 'panel_e.png', fig, axs)
    return fig


def crop_panel_e(row_pdf, dest_png, fig=None, axs=None):
    """Rasterise a four-axis row PDF at 300 dpi and cut out panel e by the producer's crop rule.

    The crop box is the one :func:`_crop` gives for axis e, converted to whole pixels
    (floor of the left edge, full height). No resampling is done. The same function makes
    ``expected/panel_e.png`` from ``expected/panels_bcde_row_v10.pdf``.

    Parameters
    ----------
    row_pdf, dest_png : pathlib.Path
    fig, axs : matplotlib objects
        The regenerated row (from :func:`_row`); it defines the crop box.

    Returns
    -------
    list of int
        Pixel box ``[x0, y0, x1, y1]`` that was cut out.
    """
    import fitz
    from PIL import Image
    if fig is None:
        raise ValueError("pass the regenerated row (fig, axs) to define the crop box")
    bb = _crop(fig, axs, 3)
    with fitz.open(row_pdf) as doc:
        pix = doc[0].get_pixmap(dpi=DPI, alpha=False, colorspace=fitz.csRGB)
        arr = np.frombuffer(pix.samples, dtype=np.uint8).reshape(pix.height, pix.width, 3)
    x0 = int(np.floor(bb.x0 * DPI))
    x1 = min(int(np.ceil(bb.x1 * DPI)), arr.shape[1])
    Image.fromarray(arr[:, x0:x1].copy()).save(dest_png)
    return [x0, 0, x1, arr.shape[0]]


def make_expected_panel_e(d, unit):
    """Write ``expected/panel_e.png`` from the copied v10 row ``expected/panels_bcde_row_v10.pdf``.

    Uses the crop box of the regenerated row; the notebook re-derives it and checks that
    the committed file is identical.

    Returns
    -------
    list of int
        Pixel box used.
    """
    fig, axs, _ = _row(d)
    unit = Path(unit)
    return crop_panel_e(unit / "expected" / "panels_bcde_row_v10.pdf", unit / "expected" / "panel_e.png", fig, axs)


# ---- panel c: build_FAR12.py ---------------------------------------------------------
def panel_c(d, out):
    """Panel c: fraction of distinct authentic or scrambled queries with an exact reservoir match.

    Moved unchanged from ``build_FAR12.py`` (132 x 143 pt axes block). Only the input reads
    and the output directory changed; the composite-overlay part of that script is
    described in ASSEMBLY.md and not repeated here.

    Parameters
    ----------
    d : dict
        Output of :func:`load`.
    out : pathlib.Path
        Directory for ``panel_c.pdf`` and ``panel_c.png``.

    Returns
    -------
    matplotlib.figure.Figure
    """
    import matplotlib.pyplot as plt
    from matplotlib import ticker
    from style import fl_style as fl
    out = Path(out)
    with plt.rc_context(_BASE_RC or None):
        plt.rcParams.update({'font.size':6.5,'axes.labelsize':6.5,'xtick.labelsize':6,
                             'ytick.labelsize':6,'legend.fontsize':5.7,'savefig.bbox':None})
        lba, fss = d["lba"], d["fss"]
        real=lba[lba.arm.eq('REAL')].sort_values('length');scr=lba[~lba.arm.eq('REAL')]
        fig=plt.figure(figsize=(132/72,143/72));ax=fig.add_axes([40/132,35/143,84/132,78/143]);_style(ax)
        series=[]
        def plot_checked(x,y,**kw):
            line,=ax.plot(x,y,**kw)
            assert np.array_equal(line.get_xdata(),x,equal_nan=True)
            assert np.array_equal(line.get_ydata(),y,equal_nan=True)
            series.append(dict(label=kw.get('label'),points=len(x)))
        for seed,g in scr.groupby('seed'):
            plot_checked(g.length,100*g.fraction_matched,color=fl.ORANGE,lw=.6,alpha=.55)
        plot_checked(fss.length,100*fss.scrambled_median,marker='o',color=fl.ORANGE,ms=2.6,lw=.9,label='Scrambled, five seeds')
        plot_checked(real.length,100*real.fraction_matched,marker='o',color=fl.FL,ms=2.6,lw=.9,label='Authentic de novo')
        ax.set_yscale('symlog',linthresh=.01,linscale=.6)
        ax.set_yticks([0,.01,.1,1,10,100]);ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x,p:f'{x:g}'))
        ax.set_ylabel('Distinct peptides with an\nexact reservoir match (%)')
        ax.set_xlabel('Peptide length (residues)')
        ax.set_xticks([10,20,30,40])
        legend=ax.legend(frameon=False,loc='lower left',bbox_to_anchor=(-.02,1.025),
                         fontsize=5.7,handlelength=1.6,borderaxespad=0)
        fig.text(40/132,1-2/143,'c',fontsize=8,color=fl.INK,va='top')
        fig.canvas.draw();renderer=fig.canvas.get_renderer()
        assert legend.get_window_extent(renderer).y0>ax.get_window_extent(renderer).y1
        assert legend.get_window_extent(renderer).y1<fig.bbox.y1
        fig.savefig(out/'panel_c.pdf');fig.savefig(out/'panel_c.png',dpi=DPI)
    return fig


# ---- statistics ----------------------------------------------------------------------
def stat_b_sizes(d):
    """Selected FASTA sizes per cohort (Fig. 1b; Results line 55).

    Statistical unit: acquisition; n = 15 CAMPI, 58 healthy stool, 306 PREDICT,
    347 CapScan (726). Descriptive only, no test.

    Returns
    -------
    dict
        ``n`` per cohort, ``total``, ``median`` per cohort, ``max``, ``n_below_human``,
        ``human_entries``.
    """
    s, h = d["sizes"], d["human"]["entries"]
    n = s.groupby("cohort").size().to_dict()
    assert n == ROSTER, n
    return dict(n=n, total=int(len(s)), median={c: float(s.loc[s.cohort.eq(c), "entries"].median()) for c in COHORTS},
                max=int(s.entries.max()), min=int(s.entries.min()),
                n_below_human=int(s.entries.lt(h).sum()), human_entries=int(h))


def stat_a_reservoir(d):
    """Reservoir size, study-level lookup and the reservoir/human ratio (Results lines 51, 53).

    Statistical unit: sequence (protein entries); pooled over 58 healthy-stool acquisitions.

    Returns
    -------
    dict
        ``reservoir``, ``reservoir_million``, ``pooled_queries``, ``study_supported``,
        ``study_supported_million``, ``reservoir_over_human``, ``hstool_acquisitions``.
    """
    st, h = d["stage1"], d["human"]["entries"]
    res = int(st["proteins_searched"])
    return dict(reservoir=res, reservoir_million=res / 1e6, pooled_queries=int(st["total_peptides"]),
                study_supported=int(st["proteins_matched"]), study_supported_million=st["proteins_matched"] / 1e6,
                reservoir_over_human=res / h,
                hstool_acquisitions=int(d["sizes"].cohort.eq("HSTOOL").sum()))


def stat_a_worked(d):
    """Worked acquisition, authentic and scrambled arms (Fig. 1a; Results lines 53, 59, 61).

    Statistical unit: one acquisition (``WORKED_ACQUISITION``); counts of distinct queries,
    proteins and peptides. Peptides accepted at peptide q <= 0.01, 5-50 residues.

    Returns
    -------
    dict
        Authentic: ``queries`` (selection grid s3_seeds and summed length_by_arm REAL),
        ``candidates``, ``candidates_million``, ``fasta``, ``peptides``, ``reservoir_over_fasta``.
        Scrambled seed 20260901: ``scr_queries``, ``scr_candidates``, ``scr_fasta``,
        ``scr_peptides``. Controls: ``n_control_workflows``, ``n_controls_zero``,
        ``n_controls_zero_all_gates``, ``controls_per_design``.
    """
    e, g, lba = d["e2e"], d["grid"], d["lba"]
    main = (e.gate.eq("q_only") & e.q.eq(Q_CUTOFF) & e.min_length.eq(5) & e.max_length.eq(50))
    a = e[main & e.arm.eq("authentic")].iloc[0]
    s = e[main & e.arm.eq("all_position") & e.seed.eq(WORKED_SEED)].iloc[0]
    row = g[g.arm.eq(TOP_FRACTION_ARM) & g["sample"].eq(WORKED_ACQUISITION)]
    assert len(row) == 1
    row = row.iloc[0]
    assert int(row.s3_proteins) == int(a.selected_proteins)  # same 6,671,447 candidates
    q_real = int(lba.loc[lba.arm.eq("REAL"), "distinct_queries"].sum())
    assert q_real == int(row.s3_seeds)
    ctrl = e[main & ~e.arm.eq("authentic")]
    ctrl_all = e[~e.arm.eq("authentic")].groupby(["arm", "seed"]).n_peptides.max()
    res = int(d["stage1"]["proteins_searched"])
    return dict(queries=int(row.s3_seeds), candidates=int(a.selected_proteins),
                candidates_million=float(a.selected_proteins / 1e6), fasta=int(a.parsimony_proteins),
                peptides=int(a.n_peptides), reservoir_over_fasta=float(res / a.parsimony_proteins),
                scr_queries=int(lba.loc[lba.arm.eq(f"SIM_{WORKED_SEED}"), "distinct_queries"].sum()),
                scr_candidates=int(s.selected_proteins), scr_fasta=int(s.parsimony_proteins),
                scr_peptides=int(s.n_peptides),
                n_control_workflows=int(len(ctrl)), n_controls_zero=int(ctrl.n_peptides.eq(0).sum()),
                n_controls_zero_all_gates=int(ctrl_all.eq(0).sum()),
                controls_per_design=ctrl.groupby("arm").size().to_dict())


def stat_c_scramble(d):
    """Scrambled-query match fractions by length over five seeds (Fig. 1c; Results line 59).

    Statistical unit: distinct query within one seed of one acquisition. n per seed and
    length is ``distinct_queries`` in ``length_by_arm.tsv``. Range over five seeds
    (minimum and maximum); not a confidence interval.

    Returns
    -------
    dict
        ``len9_min``, ``len9_max``, ``len10_min``, ``len10_max`` (percent), recomputed from
        ``length_by_arm.tsv`` and asserted equal to ``five_seed_summary.tsv``;
        ``n_ge13`` (query-seed evaluations at >= 13 residues) and ``matches_ge13``;
        ``n_seeds``; per-seed denominators at 9 and 10.
    """
    lba, fss = d["lba"], d["fss"]
    scr = lba[~lba.arm.eq("REAL")]
    assert sorted(scr.seed.astype(int).unique()) == sorted(SCRAMBLE_SEEDS)
    out = {}
    for L in (9, 10):
        x = scr[scr.length.eq(L)].fraction_matched
        f = fss[fss.length.eq(L)].iloc[0]
        assert x.min() == f.scrambled_min and x.max() == f.scrambled_max
        out[f"len{L}_min"], out[f"len{L}_max"] = float(100 * x.min()), float(100 * x.max())
        out[f"len{L}_n_per_seed"] = {int(k): int(v) for k, v in scr[scr.length.eq(L)].set_index("seed").distinct_queries.items()}
    ge = scr[scr.length.ge(SCRAMBLED_MIN_LENGTH_NO_MATCH)]
    out.update(n_ge13=int(ge.distinct_queries.sum()), matches_ge13=int(ge.matching_queries.sum()),
               n_seeds=int(scr.seed.nunique()))
    return out


def stat_d_support(d):
    """Per-cohort medians of single-query support and FL corroboration (Fig. 1d; Results line 63).

    Statistical unit: acquisition; n = 15, 58, 306, 347. Median over acquisitions of the
    per-acquisition fraction of selected proteins (population ``all_selected``).

    Returns
    -------
    dict
        ``one_top`` and ``fl_corr``: cohort -> median percent; ``*_min``/``*_max`` over cohorts.
    """
    sup = d["sup"][d["sup"].population.eq(SUPPORT_POPULATION)]
    assert sup.groupby("cohort").size().to_dict() == {"CAMPI": 15, "HSTOOL": 58, "PREDICT": 306, "CAPSCAN": 347}
    one = {k: float(v) for k, v in (100 * sup.groupby("cohort").one_top_prediction_fraction.median()).items()}
    fl = {k: float(v) for k, v in (100 * sup.groupby("cohort").with_FL_corroborated_prediction_fraction.median()).items()}
    return dict(one_top=one, fl_corr=fl, one_top_min=min(one.values()), one_top_max=max(one.values()),
                fl_corr_min=min(fl.values()), fl_corr_max=max(fl.values()))


def stat_d_not_predicted(d):
    """Accepted peptides absent from all archived AlphaNovo predictions (Results line 63).

    Statistical unit: peptide (distinct I/L-equivalent sequence), one acquisition;
    n = 29,663 accepted at peptide q <= 0.01. Counts from the A09 record of
    ``WORKED_ACQUISITION``; the fraction is computed here.

    Returns
    -------
    dict
        ``accepted``, ``not_predicted``, ``percent``.
    """
    w = d["worked"]
    assert w["sample"] == WORKED_ACQUISITION and w["status"] == "PASS"
    return dict(accepted=int(w["accepted_peptides"]), not_predicted=int(w["accepted_not_predicted_any"]),
                percent=100 * w["accepted_not_predicted_any"] / w["accepted_peptides"])


def stat_e_entrapment(d):
    """Donor-selected entrapment, median over 58 acquisitions (Fig. 1e; Results line 65).

    Statistical unit: acquisition (target); n = 58. Each acquisition value is the mean of
    three entrapment seeds (``ENTRAPMENT_SEEDS``). The seed means are recomputed here from
    ``donor_recount.tsv`` and asserted equal to ``donor_seed_means.tsv``. Median and
    quartiles use numpy (``np.median``, ``np.percentile`` with linear interpolation), as in
    the plot code.

    Returns
    -------
    dict
        For L in (5, 9): ``eL_median``, ``eL_q25``, ``eL_q75`` (percent), ``eL_n``.
        ``seeds``: sorted distinct seeds; ``max_abs_diff_seed_means``.
    """
    rc, ent = d["recount"], d["ent"]
    sel = rc[rc.tier.eq(ENTRAPMENT_TIER) & rc.q_threshold.eq(Q_CUTOFF)]
    seeds = sorted(int(x) for x in sel.seed.unique())
    assert tuple(seeds) == ENTRAPMENT_SEEDS
    assert (sel.groupby(["target", "min_length"]).size() == 3).all()
    means = sel.groupby(["target", "donor", "min_length"]).paired_estimate.mean().reset_index()
    pub = ent[ent.tier.eq(ENTRAPMENT_TIER) & ent.q_threshold.eq(Q_CUTOFF)]
    m = means.merge(pub[["target", "donor", "min_length", "paired_estimate"]], on=["target", "donor", "min_length"],
                    suffixes=("_re", "_pub"), validate="one_to_one")
    diff = float((m.paired_estimate_re - m.paired_estimate_pub).abs().max())
    assert len(m) == 116 and diff < 1e-15
    out = dict(seeds=seeds, max_abs_diff_seed_means=diff)
    for L in ENTRAPMENT_LENGTHS:
        v = 100 * pub.loc[pub.min_length.eq(L), "paired_estimate"].to_numpy()
        out.update({f"e{L}_n": int(len(v)), f"e{L}_median": float(np.median(v)),
                    f"e{L}_q25": float(np.percentile(v, 25)), f"e{L}_q75": float(np.percentile(v, 75))})
    return out


WORDS = {0: "no", 5: "five", 10: "ten"}


def as_word(n):
    """Return the English word the manuscript uses for a small count (0 -> 'no')."""
    return WORDS[int(n)]


def numbers(d):
    """Every value in NUMBERS.tsv, keyed by the row id used in that file's ``code`` column.

    Returns
    -------
    dict
        row id -> recomputed value (full precision, as computed).
    """
    b, a, w, c, s, nd, e = (stat_b_sizes(d), stat_a_reservoir(d), stat_a_worked(d), stat_c_scramble(d),
                            stat_d_support(d), stat_d_not_predicted(d), stat_e_entrapment(d))
    return {
        "L51_acquisitions": a["hstool_acquisitions"],
        "L51_reservoir_million": a["reservoir_million"],
        "L51_human": b["human_entries"],
        "L51_ratio": a["reservoir_over_human"],
        "L53_acquisitions": a["hstool_acquisitions"],
        "L53_pooled_queries": a["pooled_queries"],
        "L53_study_supported_million": a["study_supported_million"],
        "L53_queries": w["queries"],
        "L53_candidates_million": w["candidates_million"],
        "L53_fasta": w["fasta"],
        "L53_fold": w["reservoir_over_fasta"],
        "L53_peptides": w["peptides"],
        "L55_total": b["total"],
        "L55_human": b["human_entries"],
        "L55_campi_median": b["median"]["CAMPI"],
        "L55_hstool_median": b["median"]["HSTOOL"],
        "L59_len9_min": c["len9_min"],
        "L59_len9_max": c["len9_max"],
        "L59_len10_min": c["len10_min"],
        "L59_len10_max": c["len10_max"],
        "L59_no_match": as_word(c["matches_ge13"]),
        "L59_n_ge13": c["n_ge13"],
        "L59_scr_fasta": w["scr_fasta"],
        "L59_scr_peptides": as_word(w["scr_peptides"]),
        "L59_fasta": w["fasta"],
        "L59_peptides": w["peptides"],
        "L61_controls": as_word(w["n_controls_zero"]),
        "L61_per_design": as_word(sorted(set(w["controls_per_design"].values()))[0]),
        "L61_controls_no": as_word(w["n_control_workflows"] - w["n_controls_zero"]),
        "L61_peptides": w["peptides"],
        "L63_one_top_min": s["one_top_min"],
        "L63_one_top_max": s["one_top_max"],
        "L63_fl_min": s["fl_corr_min"],
        "L63_fl_max": s["fl_corr_max"],
        "L63_not_predicted": nd["not_predicted"],
        "L63_accepted": nd["accepted"],
        "L63_percent": nd["percent"],
        "L65_e5": e["e5_median"],
        "L65_e9": e["e9_median"],
        "L71_human": b["human_entries"],
        "L71_acq_a": a["hstool_acquisitions"],
        "L71_total": b["total"],
        "L71_seeds": as_word(c["n_seeds"]),
        "L71_acq_e": e["e5_n"],
    }
