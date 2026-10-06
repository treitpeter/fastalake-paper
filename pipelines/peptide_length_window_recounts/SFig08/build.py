"""Supplementary Fig. 8 at 9-50 residues: panel a redrawn from the 9-50 recount; b-d unchanged.

Uses the repo panels.py (read-only import). Steps:
1. control: draw a-d from the 5-50 recount (must equal the shipped table) with the producer code,
   redact the title, compare with expected/SFig08_v28f.png (must be 0 differing pixels).
2. 9-50: same code, panel a from peptides_9_50. Axes positions are fixed to those of the control
   render (tight_layout would otherwise be free to move b-d), so only panel a can change.
3. Pixel comparison against the v28f PNG, and panel-a statistics at 5-50 and 9-50.
"""
import os  # release: paths come from environment variables
import json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from PIL import Image

REPO = Path((os.environ.get('FASTALAKE_PAPER_REPO') or str(Path(__file__).resolve().parents[3])))
UNIT = REPO / 'figures/SFig08'
sys.path.insert(0, str(UNIT))
import panels as P  # noqa: E402

HERE = Path(__file__).resolve().parent
OUT = HERE.parent
rc5 = pd.read_csv(OUT / 'curation_counts_5_9.tsv', sep='\t')
_, abl = P.read_tables(REPO / 'source_data/SFig08')
shipped, _ = P.read_tables(REPO / 'source_data/SFig08')


def counts_at(L):
    c = rc5[['sample', 'arm', f'peptides_{L}_50']].rename(columns={f'peptides_{L}_50': 'peptides'})
    # same row order as the shipped table (groupby order is alphabetical anyway)
    return shipped[['sample', 'arm']].merge(c, on=['sample', 'arm'], how='left')


c5, c9 = counts_at(5), counts_at(9)
assert (c5.peptides.values == shipped.peptides.values).all()
P.check_tables(c9, abl)


def draw(counts, stem, positions=None):
    rc = {k: v for k, v in mpl.rcParamsDefault.items() if k not in ('backend', 'backend_fallback')}
    rc.update(P.PRODUCER_RC)
    with mpl.rc_context(rc):
        fig, axes = plt.subplots(2, 2, figsize=P.FIGSIZE); axes = axes.ravel()
        P.panel_a(axes[0], counts)
        for f, ax in zip([P.panel_b, P.panel_c, P.panel_d], axes[1:]):
            f(ax, abl)
        fig.suptitle(P.FIGURE_TITLE, fontsize=11, fontweight='bold')
        fig.tight_layout(rect=(0, 0, 1, .96), pad=1.6)
        pos = [ax.get_position().bounds for ax in axes]
        if positions is not None:
            for ax, p in zip(axes, positions):
                ax.set_position(p)
        fig.savefig(stem.with_suffix('.pdf'))
        plt.close(fig)
    png = P.redact_title(stem.with_suffix('.pdf'), stem.parent / (stem.name + '_redacted.pdf'))
    return pos, png


def cmp(png):
    a = np.asarray(Image.open(UNIT / 'expected/SFig08_v28f.png').convert('RGB')).astype(int)
    b = np.asarray(Image.open(png).convert('RGB')).astype(int)
    assert a.shape == b.shape
    d = np.abs(a - b).max(axis=2) > 0
    r, c = np.where(d.any(1))[0], np.where(d.any(0))[0]
    return dict(n_diff_px=int(d.sum()), rows=[int(r.min()), int(r.max())] if len(r) else None,
                cols=[int(c.min()), int(c.max())] if len(c) else None, shape=list(a.shape[:2]))


pos5, png5 = draw(c5, HERE / 'control_5')
rep = {'control_5_vs_v28f': cmp(png5)}
assert rep['control_5_vs_v28f']['n_diff_px'] == 0
pos9_free, png9f = draw(c9, HERE / 'len9_free')
rep['len9_free_layout_positions_equal_control'] = bool(np.allclose(pos9_free, pos5, atol=0, rtol=0))
pos9, png9 = draw(c9, HERE / 'len9', positions=pos5)
rep['len9_vs_v28f'] = cmp(png9)
# panel a box (fig fraction) -> pixels, from the control axes tight bbox
H_px, W_px = rep['len9_vs_v28f']['shape']
rep['panel_a_axes_fraction'] = pos5[0]
rep['panel_b_axes_fraction'] = pos5[1]
rep['panel_c_axes_fraction'] = pos5[2]
for L, c in ((5, c5), (9, c9)):
    rep[f'stats_{L}'] = {
        'curated_over_MG': P.stat_a_paired_gain(c, 'MG_own_curated', 'MG'),
        'FLpublic_over_curated': P.stat_a_paired_gain(c, 'FL_public', 'MG_own_curated'),
        'FLpublic_over_MG': P.stat_a_paired_gain(c, 'FL_public', 'MG'),
        'FLcohortMG_over_MG': P.stat_a_paired_gain(c, 'FL_cohortMG', 'MG'),
        'arm_medians': c.groupby('arm').peptides.median().to_dict()}
(HERE / 'stats.json').write_text(json.dumps(rep, indent=2, default=float))
print(json.dumps(rep, indent=1, default=float))
