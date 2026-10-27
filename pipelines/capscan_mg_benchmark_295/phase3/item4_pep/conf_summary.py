#!/usr/bin/env python3
"""Item 4, same-spectrum PEP retention ("84.6% in CapScan", para 326): the ordinary_summary.tsv step of the frozen
Review_v16 Confidence collector (c1_copy/collect_v2.py, lines t['gain_pct'] .. s.to_csv), statements copied verbatim.
mode control : input = Review_v16/Confidence/ordinary_all_pairs.tsv (578 pairs); output must be byte-identical to
               Review_v16/Confidence/ordinary_summary.tsv (sha256 34bf8961...).
mode n295    : input = the same 578 x 36 rows + mg295/phase2/per_acq13/confidence_rows.tsv (13 x 36 rows, the frozen
               worker_v2 reader, control 1311 36/36 PASS); derived columns recomputed with the collector's statements.
mode negctrl : n295 input with the 13 new CapScan pep01 5-50 removed_fl values raised by 5000; must NOT reproduce the n295 summary.
"""
import sys, hashlib
from pathlib import Path
import pandas as pd
H = Path(__file__).resolve().parent
CONF = Path('LAB_ROOT/review/v15_integration/ready_integration_2026-09-28/source_additions/Review_v16/Confidence')
ROWS13 = H.parents[1] / 'phase2/per_acq13/confidence_rows.tsv'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
mode = sys.argv[1]
t = pd.read_csv(CONF / 'ordinary_all_pairs.tsv', sep='\t', dtype={'subject_id': str, 'specimen_id': str})
assert sha(CONF / 'ordinary_all_pairs.tsv').startswith('47c6ab99563bf570')
t = t.drop(columns=['gain_pct', 'fl_retained', 'mg_retained'])
if mode in ('n295', 'negctrl'):
    n = pd.read_csv(ROWS13, sep='\t')
    n = n[n['index'] != 335]  # 1311 is the phase-2 control row, already in the 578
    assert n['index'].nunique() == 13 and len(n) == 13 * 36 and set(n['index']) == set(range(578, 591))
    t = pd.concat([t, n], ignore_index=True)
    if mode == 'negctrl':
        k = t.index[(t['index'] >= 578) & (t.gate == 'pep01') & (t.q == 0.01) & (t.min_length == 5) & (t.max_length == 50)]
        assert len(k) == 13
        t.loc[k, 'removed_fl'] += 5000  # planted: 13 rows' FL retention lowered
exp = {'control': 578, 'n295': 591, 'negctrl': 591}[mode]
assert t['index'].nunique() == exp and len(t) == exp * 36
keys = ['index', 'gate', 'q', 'min_length', 'max_length']; assert not t.duplicated(keys).any()
assert ((t.n_fl == t.shared + t.fl_only) & (t.n_mg == t.shared + t.mg_only)).all()
# ---- verbatim from collect_v2.py ----
t['gain_pct']=(100*(t.n_fl/t.n_mg-1)).where(t.n_mg>0)
t['fl_retained']=t.n_fl/(t.n_fl+t.removed_fl).replace(0,float('nan'))
t['mg_retained']=t.n_mg/(t.n_mg+t.removed_mg).replace(0,float('nan'))
g=['cohort','gate','q','min_length','max_length']
s=t.groupby(g).agg(n=('index','size'),defined_ratios=('gain_pct','count'),median_gain_pct=('gain_pct','median'),median_fl_retained=('fl_retained','median'),median_mg_retained=('mg_retained','median'),median_fl_only=('fl_only','median'),median_mg_only=('mg_only','median')).reset_index()
s['zero_mg']=s.n-s.defined_ratios
# ---- end verbatim ----
out = H / mode; out.mkdir(exist_ok=True); s.to_csv(out / 'ordinary_summary.tsv', sep='\t', index=False)
t.to_csv(out / 'ordinary_all_pairs_recomputed.tsv', sep='\t', index=False)
print(mode, 'summary sha256', sha(out / 'ordinary_summary.tsv'))
if mode == 'control':
    print('CONTROL', 'IDENTICAL' if sha(out / 'ordinary_summary.tsv') == sha(CONF / 'ordinary_summary.tsv') else 'DIFFERENT')
if mode == 'negctrl':
    print('NEGCTRL', 'DETECTED' if sha(out / 'ordinary_summary.tsv') != sha(H / 'n295/ordinary_summary.tsv') else 'NOT DETECTED')
print(s[(s.gate == 'pep01') & (s.q == 0.01) & (s.min_length == 5) & (s.max_length == 50)].to_string(index=False))
