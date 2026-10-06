"""Rebuild Supplementary Figure 3 with panels c, d at 9-50 residues (len950, 2026-10-06).

Usage (Slurm only): python build_sfig03_len950.py control|len950
Same steps as figures/SFig03/ASSEMBLY.md: (1) draw a-d with panels_len950.build_figure,
(2) slant the 16 taxon labels (panels_len950.slant), (3) rasterise. The control
(PLOTTED_MIN_LENGTH_CD = 5) must reproduce expected/SFig03_v28f.png at 400 dpi.
The len950 PNG is rendered at 300 dpi as requested (the v28d image was 400 dpi).
"""
import json, os, sys, hashlib
from pathlib import Path
import numpy as np
HERE = Path(__file__).resolve().parent
REPO = Path((os.environ.get('FASTALAKE_PAPER_REPO') or str(Path(__file__).resolve().parents[3])))
os.environ.setdefault('FASTALAKE_FONT_DIR', os.environ.get('FASTALAKE_FONT_DIR', ''))
sys.path.insert(0, str(HERE))
import panels_len950 as P
MODE = sys.argv[1]; assert MODE in ('control', 'len950')
P.PLOTTED_MIN_LENGTH_CD = 5 if MODE == 'control' else 9
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
OUT = HERE / 'work' / MODE; OUT.mkdir(parents=True, exist_ok=True)
assert P.init()
SD = REPO / 'source_data/SFig03'
t = P.load(SD)
rc = P.stat_recount_membership(SD, t)
assert rc['equal_support'] and rc['equal_acq_counts'] and rc['equal_recovery'] and rc['lengths_in_window']
fig, axes, plotted = P.build_figure(t)
raw = OUT / 'SFig03_raw.pdf'; fig.savefig(raw)
import matplotlib.pyplot as plt; plt.close(fig)
stem = 'SupplementaryFigure3_control' if MODE == 'control' else 'SupplementaryFigure3_len950'
final = (OUT if MODE == 'control' else HERE) / f'{stem}.pdf'
n_sl = P.slant(raw, final, set(P.short_names(list(t['support'].member.unique()))))
assert n_sl == P.N_SLANTED_LABELS
mt, k = P.min_text_pt(final); assert abs(k - 1) < 1e-6 and mt >= P.MIN_TEXT_PT
log = dict(mode=MODE, plotted_min_length_cd=P.PLOTTED_MIN_LENGTH_CD, recount=dict(equal_support=True, equal_acq_counts=True, equal_recovery=True),
           pdf=str(final), pdf_sha256=sha(final), smallest_text_pt=mt)
med = P.stat_panel_medians(plotted); log['panel_medians'] = med.to_dict('records')
log['line89'] = P.stat_line89_msfragger(t, lo=P.PLOTTED_MIN_LENGTH_CD)
plotted.to_csv(OUT / 'plotted_values.tsv', sep='\t', index=False)
if MODE == 'control':
    png = OUT / f'{stem}_400dpi.png'; P.render_png(final, png, 400)
    from PIL import Image
    a = np.asarray(Image.open(png).convert('RGB')); b = np.asarray(Image.open(REPO / 'figures/SFig03/expected/SFig03_v28f.png').convert('RGB'))
    log['control_vs_v28f_400dpi'] = dict(shape_equal=a.shape == b.shape, n_diff_px=int(np.any(a != b, axis=2).sum()) if a.shape == b.shape else None)
else:
    png = HERE / f'{stem}.png'; log['png_px'] = P.render_png(final, png, 300); log['png'] = str(png); log['png_sha256'] = sha(png)
(OUT / 'BUILD_LOG.json').write_text(json.dumps(log, indent=2, default=str) + '\n')
print(json.dumps(log, indent=1, default=str))
