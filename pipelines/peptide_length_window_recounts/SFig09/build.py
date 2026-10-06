"""Supplementary Fig. 9 at 9-50 residues: same producer code (repo panels.py, read-only), with the
FastaLake and MG peptide counts of each of the 282 acquisitions replaced by their 9-50 counts from
a08_benchmark282_old.tsv (the production-arm recount whose 5-50 columns equal the plotted counts).
Control: the 5-50 render must equal expected/SFig09_v28f.png pixel for pixel."""
import json, os, sys
from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image
REPO = Path((os.environ.get('FASTALAKE_PAPER_REPO') or str(Path(__file__).resolve().parents[3])))
UNIT = REPO / 'figures/SFig09'
sys.path.insert(0, str(UNIT))
import panels as P  # noqa: E402
HERE = Path(__file__).resolve().parent
FD = os.environ['FASTALAKE_FONT_DIR']
t = P.load(REPO / 'source_data/SFig09')
acq = t['acq282']
chk = P.check_counts_against_recount(t); assert chk == dict(rows_joined=282, unmatched=0, fl_disagree=0, mg_disagree=0, region_disagree=0), chk
m = acq.merge(t['a08_old'][['sample', 'n_fl_5', 'n_mg_5', 'n_fl_9', 'n_mg_9']], on='sample', validate='1:1')
assert len(m) == 282 and (m.n_fl_5 == m.FL_peptides).all() and (m.n_mg_5 == m.MG_peptides).all()
acq9 = m.assign(FL_peptides=m.n_fl_9, MG_peptides=m.n_mg_9).drop(columns=['n_fl_5', 'n_mg_5', 'n_fl_9', 'n_mg_9'])
acq9['gain_percent'] = (100 * (acq9.FL_peptides / acq9.MG_peptides - 1)).round(3)
acq9.to_csv(HERE.parent / 'capscan_acquisitions_9_50.tsv', sep='\t', index=False)
def cmp(png):
    a = np.asarray(Image.open(UNIT / 'expected/SFig09_v28f.png').convert('RGB')).astype(int)
    b = np.asarray(Image.open(png).convert('RGB')).astype(int)
    if a.shape != b.shape: return dict(shape_old=list(a.shape), shape_new=list(b.shape))
    d = np.abs(a - b).max(axis=2) > 0
    return dict(n_diff_px=int(d.sum()), shape=list(a.shape[:2]))
rep = {}
for L, a in ((5, acq), (9, acq9)):
    spec, caps, summary = P.derive_tables(a)
    fig = P.figure_ab(spec, summary, font_dir=FD); stem = HERE / f'sfig09_{L}'; P.save(fig, stem)
    import matplotlib.pyplot as plt; plt.close(fig)
    spec.to_csv(HERE / f'per_specimen_{L}.tsv', sep='\t', index=False); summary.to_csv(HERE / f'tertiles_{L}.tsv', sep='\t', index=False)
    rep[L] = dict(vs_v28f=cmp(f'{stem}.png'), tertiles=summary.to_dict('records'),
                  n_specimens=int(len(spec)), spec_gain_max=float(spec.gain_percent.max()), spec_gain_min=float(spec.gain_percent.min()),
                  pooled_yield=P.stat_129_pooled_yield(a))
assert rep[5]['vs_v28f'].get('n_diff_px') == 0, rep[5]
(HERE / 'stats.json').write_text(json.dumps(rep, indent=2, default=str))
print(json.dumps(rep, indent=1, default=str))
