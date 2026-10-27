#!/usr/bin/env python3
"""Execute a unit notebook's code cells in-process (cwd = unit dir, as tools/run_notebooks.py does) and harvest the
recomputed numbers in NUMBERS.tsv row order, without the cells' "equals the committed files" assertions stopping the
harvest (a 295 run is expected to differ from the committed 282 values). A cell that raises is recorded and execution
continues with the next cell.

usage: extract.py <repo copy> <unit path, e.g. figures/Fig4>  ->  <repo>/HARVEST_<unit>.tsv
columns: NUMBERS.tsv columns + recomputed_now + cell_errors
"""
import csv, json, os, sys, traceback, io, contextlib
from pathlib import Path
repo, unit = Path(sys.argv[1]).resolve(), sys.argv[2]
U = repo / unit; os.chdir(U); sys.path[:0] = [str(repo), str(U)]
import matplotlib; matplotlib.use('Agg')
nb = json.load(open(next(U.glob('*.ipynb'))))
ns = {'__name__': '__main__', 'display': lambda *a, **k: None}
errors = []
for i, c in enumerate(nb['cells']):
    if c['cell_type'] != 'code':
        continue
    src = ''.join(c['source'])
    src = '\n'.join(l for l in src.splitlines() if not l.lstrip().startswith(('%', '!')))
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            exec(compile(src, f'cell{i}', 'exec'), ns)
    except BaseException as e:   # noqa: BLE001 (SystemExit from value checks included)
        errors.append(f'cell{i}: {type(e).__name__}: {str(e)[:300]}')
num = list(csv.DictReader(open(U / 'NUMBERS.tsv', encoding='utf-8'), delimiter='\t'))
if U.name == 'SFig06':
    num = [r for r in num if r['recomputed_value'] not in ('', 'NA')]   # rows the notebook recomputes
name = U.name
(repo / f'HARVEST_{name}.errors.txt').write_text('\n'.join(errors) + '\n')
vals = None
if name == 'Fig4':
    vals = [v for _, v in ns['rows']]
elif name == 'SFig06':
    vals = [v for _, v in ns['numbers']]
elif name == 'SFig09':
    vals = [v for _, _, v in ns['rows']]
elif name == 'supplementary_methods_6_7':
    vals = list(ns['rows']['recomputed_value'])
elif name == 'abstract_discussion':
    st = ns['stats']; vals = [v for _, v in st.numbers(ns['gain'], ns['sizes'], ns['times'], ns['n_members'])]
elif name == 'methods':
    R = ns['R']; vals = [R[r['printed_value']] for r in num]
fmt = lambda v: (repr(v.item() if hasattr(v, 'item') else v) if isinstance(getattr(v, 'item', lambda: v)(), float) else str(v.item() if hasattr(v, 'item') else v))
if name in ('Fig4', 'abstract_discussion', 'SFig06'):
    fmt = lambda v: repr(v.item() if hasattr(v, 'item') else v) if not isinstance(v, str) else v
if len(vals) != len(num):
    errors.append(f'harvest length {len(vals)} != NUMBERS rows {len(num)}; aligned by order up to the shorter')
out = repo / f"HARVEST_{name}.tsv"
with open(out, 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=list(num[0]) + ['recomputed_now'], delimiter='\t', lineterminator='\n')
    w.writeheader()
    for r, v in zip(num, vals):
        w.writerow(dict(r, recomputed_now=fmt(v)))
(repo / f"HARVEST_{name}.errors.txt").write_text('\n'.join(errors) + '\n')
print(name, 'rows', min(len(vals), len(num)), 'errors', len(errors))
