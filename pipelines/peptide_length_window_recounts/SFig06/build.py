"""Supplementary Fig. 6 at 9-50 residues: rebuild b,c and d,e with only the length filter changed.

Run twice by build.sbatch: LEN950_MIN=5 (control, must reproduce the v28f composite pixel for pixel)
and LEN950_MIN=9. Then writes stats and the pixel comparison.
"""
import json, os, sys
from pathlib import Path
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import panels_len950 as P  # noqa: E402

OUT = HERE.parent
EXP = P.UNIT / 'expected'
L = P.MIN_LENGTH
d = HERE / f'render_{L}_{P.ROSTER}'
d.mkdir(exist_ok=True)
P.panel_abc(d); P.panel_de(d); res = P.compose(d)

old = np.asarray(Image.open(EXP / 'SFig06_v28f.png').convert('RGB')).astype(int)
new = np.asarray(Image.open(res['png']).convert('RGB')).astype(int)
assert old.shape == new.shape, (old.shape, new.shape)
diff = np.abs(old - new).max(axis=2) > 0
rows = np.where(diff.any(axis=1))[0]
cols = np.where(diff.any(axis=0))[0]
px_per_pt = new.shape[1] / 669.6
rep = dict(min_length=L, roster=P.ROSTER, png_px=list(new.shape[:2][::-1]), n_diff_px=int(diff.sum()),
           diff_rows_px=[int(rows.min()), int(rows.max())] if len(rows) else None,
           diff_cols_px=[int(cols.min()), int(cols.max())] if len(cols) else None,
           diff_top_pt=float(rows.min() / px_per_pt) if len(rows) else None)
rep['split_710'] = P.stat_split_text(P.load_split(710))
rep['split_726'] = P.stat_split_text(P.load_split(726))
rep['bench'] = P.stat_bench(P.load_benchmark())
(HERE / f'stats_{L}_{P.ROSTER}.json').write_text(json.dumps(rep, indent=2, default=float))
print(json.dumps(rep, indent=1, default=float))
