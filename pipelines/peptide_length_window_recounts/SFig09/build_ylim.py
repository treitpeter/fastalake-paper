"""Supplementary Fig. 9 at 9-50 with one change to the producer figure: the panel-a gain axis lower
limit is lowered from 0 to the smallest value that keeps every marker fully inside the axes
(one specimen has gain -0.03% at 9-50). figure_ab is repeated here call for call; only that
ylim is added after panel_a. Control: with the 5-50 data the lower limit stays 0 and the
render equals expected/SFig09_v28f.png."""
import json, os, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image
REPO = Path((os.environ.get('FASTALAKE_PAPER_REPO') or str(Path(__file__).resolve().parents[3])))
UNIT = REPO / 'figures/SFig09'; sys.path.insert(0, str(UNIT))
import panels as P
HERE = Path(__file__).resolve().parent; FD = os.environ['FASTALAKE_FONT_DIR']
t = P.load(REPO / 'source_data/SFig09')
acq9 = pd.read_csv(HERE.parent / 'capscan_acquisitions_9_50.tsv', sep='\t')
MARKER_S = 11  # panel_a scatter size (pt^2)

def lowest_ok(ax, fig, spec):
    """Smallest |ymin| step (0.01 % units) so the lowest marker's full disc lies above the x axis."""
    r_px = np.sqrt(MARKER_S) / 2 * fig.dpi / 72
    ymin = 0.0
    lo = spec.gain_percent.min()
    while True:
        ax.set_ylim(ymin, spec.gain_percent.max() * 1.2)
        y_pt = ax.transData.transform((1, lo))[1]; y_ax = ax.transAxes.transform((0, 0))[1]
        if y_pt - r_px >= y_ax - 1e-9:
            return ymin
        ymin -= 0.01

def figure(spec, summary):
    plt.rcdefaults()
    family = P._font_family(FD)
    rc = {'font.size': 8, 'axes.titlesize': 9, 'pdf.fonttype': 42, 'svg.fonttype': 'none'}
    if family: rc['font.family'] = family
    plt.rcParams.update(rc)
    fig, axes = plt.subplots(1, 2, figsize=P.FIGSIZE, gridspec_kw={'width_ratios': [1.5, 1]})
    P.panel_a(axes[0], spec)
    ymin = lowest_ok(axes[0], fig, spec) if spec.gain_percent.min() < 0 else 0.0
    axes[0].set_ylim(ymin, spec.gain_percent.max() * 1.2)
    P.panel_b(axes[1], summary)
    for ax in axes: ax.spines[['top', 'right']].set_visible(False)
    fig.tight_layout()
    return fig, ymin

def cmp(png):
    a = np.asarray(Image.open(UNIT / 'expected/SFig09_v28f.png').convert('RGB')).astype(int)
    b = np.asarray(Image.open(png).convert('RGB')).astype(int)
    return dict(n_diff_px=int((np.abs(a - b).max(axis=2) > 0).sum()))

rep = {}
for L, a in ((5, t['acq282']), (9, acq9)):
    spec, caps, summary = P.derive_tables(a)
    fig, ymin = figure(spec, summary); stem = HERE / f'sfig09_ylim_{L}'; P.save(fig, stem); plt.close(fig)
    rep[L] = dict(ymin=ymin, vs_v28f=cmp(f'{stem}.png'), spec_gain_min=float(spec.gain_percent.min()))
assert rep[5]['vs_v28f']['n_diff_px'] == 0 and rep[5]['ymin'] == 0
# against the previous 9-50 render: only panel a may differ
a = np.asarray(Image.open(HERE / 'sfig09_9.png').convert('RGB')).astype(int)
b = np.asarray(Image.open(HERE / 'sfig09_ylim_9.png').convert('RGB')).astype(int)
d = np.abs(a - b).max(axis=2) > 0; c = np.where(d.any(0))[0]
rep['ylim9_vs_prev9'] = dict(n_diff_px=int(d.sum()), cols=[int(c.min()), int(c.max())] if len(c) else None, width=int(a.shape[1]))
(HERE / 'stats_ylim.json').write_text(json.dumps(rep, indent=1)); print(json.dumps(rep, indent=1))
