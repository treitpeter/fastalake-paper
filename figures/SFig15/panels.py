"""Supplementary Figure 15: measured computational resources.

``figure_v30`` draws the figure as placed from manuscript v30 on: the same values in the house
style (port of the v30 restyle script ``build_sf15_style.py`` and ``sf_common.py``).

Four heatmaps (panels a-d) of seven workflow stages by four cohorts, drawn together in
one matplotlib figure by the original producer: a, median task execution minutes; b,
median usable sampled memory (GiB); c, sum of allocated CPU-hours; d, usable memory
samples / tasks. Shading is log10(1 + value); printed values are untransformed.

Source data
-----------
source_data/SFig15/compute_stage_summary.tsv
    One row per cohort x stage (28 rows); the table the figure draws.
source_data/SFig15/compute_task_records.tsv
    One row per Slurm task (2,905 rows); the stage summary is recomputed from it by
    ``stage_summary`` and must equal the shipped summary.
source_data/SFig15/campi_alphanovo_inference_sacct.tsv
    Slurm accounting of the CAMPI AlphaNovo inference arrays (text line 299 only).
source_data/Fig1/acquisition_sizes_726.tsv
    Reservoir sizes per cohort (text line 193 only); owned by Fig1, not duplicated.

Producer
--------
``draw_figure`` is the figure block of the producer ``build_revision.py`` (lines 13,
34-36, 51, 59, 229-237), as copied into ``build_sf15.py`` with the summary read from a
table and the x tick labels derived from it. ``redact_title`` is ``redact_sf15.py`` (the
logic of ``prepare_supplement_figures.py``). ``stage_summary`` is the summary loop of
the producer's Computing block (lines 224-228). Logic and values are unchanged; paths are
arguments and literal values are named constants below.
"""

import collections
import hashlib
import re
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# ---- choices (named constants) ------------------------------------------------------
SUMMARY_SHA256 = "37c1d94d821d5ba9defb5ff5c8495905ce93209a136952dec8b600b74816498f"
RECORDS_SHA256 = "fb248054124262d7deff9a9e7653b85d6c4f055c08081cc8b1c29c15ff458c2a"
SACCT_SHA256 = "d007807d2819bc8b6213c9ee05aebd6d02428beb9ab640197a3b75d2eb0d6029"
SIZES_SHA256 = "fd9ade3bb4ebd2ee6c129f69819b9697cf20abd7c10c77a12644bef25396e386"   # Fig1 table

STAGES = ['02_evidence', '05_persample', '06_parsimony', '07_sage', '08_aggregate',
          '09_entrap_shuffle', '10_fdp_shuffle']            # producer l.229; Methods v28d line 903
STAGE_NAMES = ['Shared lookup', 'Acquisition lookup', 'Parsimony', 'Sage', 'Aggregation',
               'Shuffle search', 'Shuffle diagnostic']      # producer l.230 (y tick labels)
COHORTS = ['CAMPI', 'HSTOOL', 'PREDICT', 'CAPSCAN']         # producer l.230 (x order)
PER_ACQUISITION_STAGES = {'05_persample', '06_parsimony', '07_sage', '09_entrap_shuffle'}  # producer l.203
N_STAGE = '07_sage'               # stage whose task count is the acquisition n (x tick labels)
PANEL_COLUMNS = [('wall_median_min', 'Median execution time (min)'),
                 ('memory_median_GiB', 'Median sampled memory (GiB)'),
                 ('allocated_CPU_hours', 'Allocated CPU-hours (sum)'),
                 ('n_memory_sampled', 'Tasks with usable memory samples')]  # producer l.231
SHADE_TEXT_THRESHOLD = .62        # producer l.235: white text above this shade
CMAP = 'Blues'                    # producer l.232
FIGSIZE_IN = (7.1, 6.9)           # producer l.230
RC = {'font.family': 'DejaVu Sans', 'font.size': 8, 'axes.titlesize': 9, 'axes.labelsize': 8,
      'xtick.labelsize': 7, 'ytick.labelsize': 7, 'pdf.fonttype': 42, 'svg.fonttype': 'none',
      'axes.spines.top': False, 'axes.spines.right': False, 'savefig.dpi': 180}   # producer l.13
SUPTITLE = 'Computational resources across manuscript cohorts'   # removed again by redact_title
RASTER_DPI = 300                  # prepare_supplement_figures.py raster step
SELECTION_STAGES = ('05_persample', '06_parsimony')   # "database selection" in v28d line 193 (see OPEN_ITEMS)
LOOKUP_STAGE = '02_evidence'      # "study-level lookup" (v28d line 193), "pooled extraction" (line 299)
LAKE_STAGE = 'lake'               # CAMPI reservoir construction (v28d line 299)
CAMPI_INFERENCE_ARRAY = '5450700' # CAMPI AlphaNovo inference array, 15 tasks (v28d line 299)

_DEFAULT_RC = {k: v for k, v in mpl.rcParamsDefault.items() if k != "backend"}


def check_sha256(path, expected):
    """Raise if the sha256 of ``path`` differs from ``expected``.

    Parameters
    ----------
    path : pathlib.Path
        File to check.
    expected : str
        Hex sha256.

    Returns
    -------
    str
        The sha256.
    """
    s = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    if s != expected:
        raise ValueError(f"{Path(path).name}: sha256 {s} differs from expected {expected}")
    return s


def load_summary(path):
    """Read the stage summary (28 rows, cohort x stage); floats parsed round-trip exact."""
    return pd.read_csv(path, sep='\t', float_precision='round_trip')


def load_records(path):
    """Read the task records (one row per Slurm task); floats parsed round-trip exact."""
    return pd.read_csv(path, sep='\t', dtype={'jobid': str}, float_precision='round_trip')


def load_sacct(path):
    """Read the CAMPI inference sacct export (job, batch and extern step rows)."""
    return pd.read_csv(path, sep='\t', dtype={'JobID': str})


def stage_summary(jobs):
    """Recompute the stage summary from task records (producer Computing block, l.224-228).

    Statistical unit: Slurm task (one task per acquisition for per-acquisition stages, one
    task per cohort for shared stages). Descriptive only.

    Parameters
    ----------
    jobs : pandas.DataFrame
        Output of ``load_records``.

    Returns
    -------
    pandas.DataFrame
        Same columns and row order as ``compute_stage_summary.tsv``.
    """
    stages = STAGES
    summary = []
    warnings.filterwarnings('ignore', 'Mean of empty slice', RuntimeWarning)  # median of the empty CAMPI shuffle stages
    for co in COHORTS:
        for st in stages:
            g = jobs[jobs.cohort.eq(co) & jobs.stage.eq(st)]
            summary.append(dict(cohort=co, stage=st, n_tasks=len(g),
                                n_memory_sampled=int(g.memory_GiB.notna().sum()),
                                wall_median_min=g.wall_minutes.median(),
                                memory_median_GiB=g.memory_GiB.median(),
                                allocated_CPU_hours=g.allocated_CPU_hours.sum() if len(g) else np.nan,
                                population_scope=';'.join(g.population_scope.unique())))
    return pd.DataFrame(summary)


def draw_figure(s, out_pdf):
    """Draw panels a-d exactly as the original producer (with its in-image title).

    Parameters
    ----------
    s : pandas.DataFrame
        Output of ``load_summary``.
    out_pdf : pathlib.Path
        PDF to write (producer name ``COMPUTE_RESOURCE_MAP.pdf``).

    Returns
    -------
    list of str
        The x tick labels (cohort and acquisition n).
    """
    out_pdf = Path(out_pdf)
    with mpl.rc_context(_DEFAULT_RC):
        plt.rcParams.update(RC)

        def title(ax, l, t):
            ax.set_title(l + '  ' + t, loc='left', pad=9)

        def save(fig, title):
            fig.suptitle(title, fontsize=11, fontweight='bold')
            fig.tight_layout(rect=(0, 0, 1, .96), pad=1.6)
            fig.savefig(out_pdf)
            plt.close(fig)

        # ---- moved from build_sf15.py (figure block of build_revision.py), unchanged logic ----
        stages = STAGES
        n = {co: int(s[s.cohort.eq(co) & s.stage.eq(N_STAGE)].n_tasks.iloc[0]) for co in COHORTS}
        xt = [f'CAMPI\n{n["CAMPI"]}', f'Stool\n{n["HSTOOL"]}', f'PREDICT\n{n["PREDICT"]}', f'CapScan\n{n["CAPSCAN"]}']
        fig, axes = plt.subplots(2, 2, figsize=FIGSIZE_IN)
        axes = axes.ravel()
        coh = COHORTS
        namesst = STAGE_NAMES
        for i, (col, lab) in enumerate(PANEL_COLUMNS):
            mat = s.pivot(index='stage', columns='cohort', values=col).reindex(index=stages, columns=coh)
            vals = mat.values.astype(float)
            im = axes[i].imshow(np.ma.masked_invalid(np.log10(1 + vals)), cmap=CMAP, aspect='auto')
            for y, st in enumerate(stages):
                for x, co in enumerate(coh):
                    row = s[s.cohort.eq(co) & s.stage.eq(st)].iloc[0]
                    v = vals[y, x]
                    label = ('NA' if not np.isfinite(v) or row.n_tasks == 0 else f'{v:.2g}') if i < 3 else (
                        f'{row.n_memory_sampled}/{row.n_tasks}' if row.n_tasks else 'NA')
                    shade = np.log10(1 + v) / np.nanmax(np.log10(1 + vals)) if np.isfinite(v) else 0
                    axes[i].text(x, y, label, ha='center', va='center', fontsize=8,
                                 color='white' if shade > SHADE_TEXT_THRESHOLD and label != 'NA' else 'black')
            axes[i].set_xticks(range(4), xt)
            axes[i].set_yticks(range(7), namesst)
            title(axes[i], 'abcd'[i], lab)
        save(fig, SUPTITLE)
    return xt


def redact_title(in_pdf, out_pdf, out_png):
    """Remove the in-image title, check numeric labels are unchanged, save PDF and PNG.

    Logic of ``redact_sf15.py`` (``prepare_supplement_figures.py``): the title is found by
    its text, covered by a white redaction, the multiset of numeric strings on the page is
    asserted unchanged, the PDF is saved with ``garbage=4, deflate=True`` and rendered at
    ``RASTER_DPI``.

    Parameters
    ----------
    in_pdf : pathlib.Path
        Output of ``draw_figure``.
    out_pdf, out_png : pathlib.Path
        Composite PDF and PNG to write.

    Returns
    -------
    dict
        Title bounding box (points) and page rectangle.
    """
    import fitz
    doc = fitz.open(in_pdf)
    p = doc[0]
    before = p.get_text()
    ttext = SUPTITLE
    spans = [sp for b in p.get_text('dict')['blocks'] for l in b.get('lines', []) for sp in l['spans']
             if sp['text'].strip() == ttext]
    assert len(spans) == 1
    bbox = spans[0]['bbox']
    assert bbox[1] < 20
    p.add_redact_annot(fitz.Rect(bbox), fill=(1, 1, 1))
    p.apply_redactions(images=0, graphics=0)
    assert ttext not in p.get_text()
    before_nums = collections.Counter(re.findall(r'\d[\d,.]*', before.replace(ttext, '')))
    after_nums = collections.Counter(re.findall(r'\d[\d,.]*', p.get_text()))
    assert before_nums == after_nums
    doc.save(out_pdf, garbage=4, deflate=True)
    p = fitz.open(out_pdf)[0]
    p.get_pixmap(dpi=RASTER_DPI, alpha=False).save(out_png)
    return dict(title_bbox=[round(v, 3) for v in bbox], rect=list(p.rect))


def panel_labels(s, letter):
    """Return the values printed in one panel as a stage x cohort table (display aid).

    Uses the same label rule as ``draw_figure``: panels a-c print ``f'{v:.2g}'`` (``NA``
    for no task), panel d prints ``usable memory samples/tasks``.

    Parameters
    ----------
    s : pandas.DataFrame
        Output of ``load_summary``.
    letter : str
        ``'a'``, ``'b'``, ``'c'`` or ``'d'``.

    Returns
    -------
    pandas.DataFrame
        Rows: stage display names; columns: cohorts.
    """
    i = 'abcd'.index(letter)
    col = PANEL_COLUMNS[i][0]
    out = {}
    for co in COHORTS:
        lab = []
        for st in STAGES:
            row = s[s.cohort.eq(co) & s.stage.eq(st)].iloc[0]
            v = float(row[col]) if pd.notna(row[col]) else np.nan
            if i < 3:
                lab.append('NA' if not np.isfinite(v) or row.n_tasks == 0 else f'{v:.2g}')
            else:
                lab.append(f'{row.n_memory_sampled}/{row.n_tasks}' if row.n_tasks else 'NA')
        out[co] = lab
    return pd.DataFrame(out, index=STAGE_NAMES)


def stat_acquisition_counts(s):
    """Acquisitions per cohort in the per-acquisition stages, and shared-stage scopes.

    Statistical unit: acquisition. n = tasks of ``N_STAGE`` per cohort; shared scopes are
    the acquisition counts written in ``population_scope`` of the shared stages.

    Parameters
    ----------
    s : pandas.DataFrame
        Output of ``load_summary``.

    Returns
    -------
    dict
        ``n`` (cohort -> acquisitions; also checked equal across all per-acquisition
        stages) and ``shared_scope`` (cohort -> acquisitions in the shared stage scope).
    """
    n = {}
    for co in COHORTS:
        per = s[s.cohort.eq(co) & s.stage.isin(PER_ACQUISITION_STAGES) & s.n_tasks.gt(0)]
        vals = set(per.n_tasks)
        assert len(vals) == 1, (co, vals)
        n[co] = int(s[s.cohort.eq(co) & s.stage.eq(N_STAGE)].n_tasks.iloc[0])
    shared = {}
    for co in COHORTS:
        sc = s[s.cohort.eq(co) & ~s.stage.isin(PER_ACQUISITION_STAGES)].population_scope.dropna()
        m = {int(x) for t in sc for x in re.findall(r'shared stage as run: (\d+) acquisitions', t)}
        if m:
            assert len(m) == 1, (co, m)
            shared[co] = m.pop()
    return dict(n=n, shared_scope=shared)


def stat_stage_medians(s):
    """Median wall minutes per cohort and stage, as a nested dict.

    Statistical unit: Slurm task; n = ``n_tasks`` of each cell. Descriptive only.

    Parameters
    ----------
    s : pandas.DataFrame
        Output of ``load_summary``.

    Returns
    -------
    dict
        ``{cohort: {stage: median minutes}}``.
    """
    return {co: {st: float(s[s.cohort.eq(co) & s.stage.eq(st)].wall_median_min.iloc[0]) for st in STAGES}
            for co in COHORTS}


def stat_lookup_range(s):
    """Range over the four cohorts of the study-level lookup wall time (``LOOKUP_STAGE``).

    Statistical unit: cohort (one shared task each); n = 4.

    Returns
    -------
    dict
        ``min``, ``max`` (minutes) and the cohorts at which they occur.
    """
    g = s[s.stage.eq(LOOKUP_STAGE)].set_index('cohort').wall_median_min
    return dict(min=float(g.min()), max=float(g.max()), argmin=g.idxmin(), argmax=g.idxmax())


def stat_selection_per_acquisition(jobs):
    """Median per-acquisition database-selection time, per cohort.

    Selection time of one acquisition = wall minutes of its ``SELECTION_STAGES`` tasks
    summed. Tasks are matched by pseudonymous acquisition id, or for CAMPI (no id in the
    records) by Slurm array index, which is the same acquisition in both arrays.
    Statistical unit: acquisition; n = 15, 58, 306, 347.

    Parameters
    ----------
    jobs : pandas.DataFrame
        Output of ``load_records``.

    Returns
    -------
    dict
        ``per_acquisition_median`` (cohort -> minutes), ``n`` (cohort -> acquisitions),
        ``sum_of_stage_medians`` (cohort -> minutes; the alternative reading) and the
        ranges of both over cohorts.
    """
    g = jobs[jobs.stage.isin(SELECTION_STAGES) & jobs.cohort.isin(COHORTS)].copy()
    g['key'] = g.acquisition.where(g.acquisition.notna(), g.jobid.str.split('_').str[1])
    med, n, alt = {}, {}, {}
    for co in COHORTS:
        x = g[g.cohort.eq(co)]
        per = x.groupby('key').agg(total=('wall_minutes', 'sum'), k=('stage', 'nunique'))
        assert (per.k == len(SELECTION_STAGES)).all(), co
        med[co] = float(per.total.median())
        n[co] = int(len(per))
        alt[co] = float(sum(x[x.stage.eq(st)].wall_minutes.median() for st in SELECTION_STAGES))
    return dict(per_acquisition_median=med, n=n, sum_of_stage_medians=alt,
                range_per_acquisition=(min(med.values()), max(med.values())),
                range_sum_of_medians=(min(alt.values()), max(alt.values())))


def stat_sage_max(s):
    """Largest cohort median Sage wall time (minutes); the text says less than 1 min.

    Statistical unit: acquisition (Sage task); n = 15, 58, 306, 347.
    """
    return float(s[s.stage.eq(N_STAGE)].wall_median_min.max())


def stat_campi_lake(jobs):
    """CAMPI reservoir construction job: wall minutes and program threads.

    Statistical unit: one Slurm job (the CAMPI reservoir build); n = 1.

    Returns
    -------
    dict
        ``minutes``, ``program_threads``, ``allocated_logical_cpus``, ``jobid``.
    """
    r = jobs[jobs.cohort.eq('CAMPI') & jobs.stage.eq(LAKE_STAGE)]
    assert len(r) == 1
    r = r.iloc[0]
    return dict(minutes=float(r.elapsed_s) / 60, program_threads=int(r.program_threads),
                allocated_logical_cpus=int(r.allocated_logical_cpus), jobid=str(r.jobid))


def stat_campi_inference(sacct):
    """Median AlphaNovo wall minutes per CAMPI acquisition (array ``CAMPI_INFERENCE_ARRAY``).

    Uses the job rows (not the ``.batch``/``.extern`` steps) of the array. Statistical
    unit: acquisition; n = 15.

    Returns
    -------
    dict
        ``median_minutes``, ``n``, ``min_minutes``, ``max_minutes``.
    """
    j = sacct[sacct.JobID.str.match(rf'^{CAMPI_INFERENCE_ARRAY}_\d+$')]
    m = j.ElapsedRaw.astype(float) / 60
    return dict(median_minutes=float(m.median()), n=int(len(m)), min_minutes=float(m.min()),
                max_minutes=float(m.max()))


def stat_reservoir_range(sizes):
    """Smallest and largest cohort reservoir (millions of sequences), text line 193.

    Statistical unit: sequence; one reservoir per cohort (4 cohorts).

    Parameters
    ----------
    sizes : pandas.DataFrame
        ``source_data/Fig1/acquisition_sizes_726.tsv``.

    Returns
    -------
    dict
        ``min_millions``, ``max_millions``.
    """
    r = sizes.groupby('cohort').reservoir_entries.agg(['nunique', 'first'])
    assert (r['nunique'] == 1).all()
    return dict(min_millions=float(r['first'].min()) / 1e6, max_millions=float(r['first'].max()) / 1e6)


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


# v30 house-style figure (LAB_ROOT/.../final/v30/figures_style/SFig15/build_sf15_style.py)
V30_SIZE_MM = (160, 155.5)       # placed size in the v29b docx
V30_CELL_PT = 7                  # cell text (producer: 8 pt)
V30_TICK_LEN = 2.6
V30_H_PAD, V30_W_PAD, V30_SPACE = 4 / 72, 2 / 72, 0.06   # constrained-layout pads and spacing


def _v30_lum(rgb):
    """Relative luminance (WCAG) of an sRGB triple in 0..1."""
    f = lambda c: c / 12.92 if c <= 0.03928 else ((c + 0.055) / 1.055) ** 2.4
    r, g, b = (f(c) for c in rgb)
    return 0.2126 * r + 0.7152 * g + 0.0722 * b


def _v35_plain(t):
    """v35: a '.2g' cell value printed in scientific notation (e.g. '1.4e+02') as a plain number ('140')."""
    return f'{float(t):.0f}' if 'e+' in t else t


def figure_v30(s, out_stem, ref_pdf):
    """Draw Supplementary Figure 15 in the v30 house style and prove its values unchanged.

    Port of ``build_sf15_style.py``: the producer (``draw_figure``) is run with savefig
    intercepted (reference fingerprint); the same table is redrawn with fl_style (Söhne,
    160 x 155.5 mm, white-to-blue scale on log10(1 + value), cell text 7 pt in white or ink by
    higher contrast, panel titles removed: their content is the legend's a-d sentences); every
    heatmap array, colour limit and printed value is compared exactly with the reference.

    Parameters
    ----------
    s : pandas.DataFrame
        Stage summary (``load_summary``).
    out_stem : pathlib.Path
        Output path without suffix; ``.pdf`` and 300 dpi ``.png`` are written.
    ref_pdf : pathlib.Path
        Path passed to the intercepted producer (nothing is written there).

    Returns
    -------
    dict
        ``panels`` (comparison report lines), ``clip_ok``. Statistical unit: Slurm task;
        descriptive, no test.
    """
    from matplotlib.colors import LinearSegmentedColormap, to_rgb
    ref = v30_capture(lambda: draw_figure(s, ref_pdf))
    assert len(ref) == 1
    ref = ref[0]
    with mpl.rc_context():
        fl = _v30_fl(); INK = fl.INK
        stages = STAGES
        n = {co: int(s[s.cohort.eq(co) & s.stage.eq(N_STAGE)].n_tasks.iloc[0]) for co in COHORTS}
        xt = [f'CAMPI\n{n["CAMPI"]}', f'Healthy stool\n{n["HSTOOL"]}', f'PREDICT\n{n["PREDICT"]}', f'CapScan\n{n["CAPSCAN"]}']   # v35: 'Stool' -> 'Healthy stool'
        for r in ref:   # the producer's tick strings and cell texts, mapped to the v35 forms for the exact comparison
            r['texts'] = sorted((_v35_plain(t), x) for t, x in r['texts'])
            r['xticklabels'] = [t.replace('Stool\n', 'Healthy stool\n') for t in r['xticklabels']]
        coh = COHORTS
        namesst = STAGE_NAMES
        # white -> Okabe-Ito blue (fl.seq_cmap without its navy end, so ink and white text both stay legible)
        cmap = LinearSegmentedColormap.from_list('fl_seq_blue', ['#ffffff', fl.BLUE])
        fig, axes = plt.subplots(2, 2, figsize=fl.mm(*V30_SIZE_MM), layout='constrained')
        fig.get_layout_engine().set(w_pad=V30_W_PAD, h_pad=V30_H_PAD, wspace=V30_SPACE, hspace=V30_SPACE)
        axes = axes.ravel()
        for i, (col, lab) in enumerate(PANEL_COLUMNS):
            ax = axes[i]; fl.style(ax)
            mat = s.pivot(index='stage', columns='cohort', values=col).reindex(index=stages, columns=coh); vals = mat.values.astype(float)
            im = ax.imshow(np.ma.masked_invalid(np.log10(1 + vals)), cmap=cmap, aspect='auto')
            lo, hi = im.get_clim()
            for y, st in enumerate(stages):
                for x, co in enumerate(coh):
                    row = s[s.cohort.eq(co) & s.stage.eq(st)].iloc[0]; v = vals[y, x]
                    label = ('NA' if not np.isfinite(v) or row.n_tasks == 0 else _v35_plain(f'{v:.2g}')) if i < 3 else (f'{row.n_memory_sampled}/{row.n_tasks}' if row.n_tasks else 'NA')
                    if np.isfinite(v):
                        bg = to_rgb(cmap((np.log10(1 + v) - lo) / (hi - lo) if hi > lo else 0.0))
                    else:
                        bg = (1, 1, 1)
                    ax.text(x, y, label, ha='center', va='center', fontsize=V30_CELL_PT,
                            color='white' if (label != 'NA' and 1.05 / (_v30_lum(bg) + .05) > (_v30_lum(bg) + .05) / (_v30_lum(to_rgb(INK)) + .05)) else INK)   # higher contrast
            ax.set_xticks(range(4), xt); ax.set_yticks(range(7), namesst)
            ax.tick_params(length=V30_TICK_LEN)
            for sp in ax.spines.values():
                sp.set_visible(False)
            fl.panel(ax, 'abcd'[i])
        fig.canvas.draw()
        rep = v30_compare(ref, v30_fingerprint(fig))
        ok = _v30_finish(fl, fig, out_stem, letters=list(axes))
    assert ok
    return dict(panels=rep, clip_ok=ok)
