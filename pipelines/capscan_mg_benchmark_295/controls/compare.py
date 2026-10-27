#!/usr/bin/env python3
"""Control (existing 282 pair 20210513_EXPL0_PeTr_SA_T760_1311, roster i=335) and pilot completeness.
Each mg295 step re-run on the control must reproduce the archived accepted-peptide COUNT exactly (sets are
reported too). Archives: (b) Review_v16 Confidence frozen MG sets; (d) iteration3/disjoint MG_ODD;
(c) public347 cross-fit pooled FL/FULL/MG (same 347 caches); (e) second_engine_all MG receipt.
usage: compare.py [pilot sample ...]  -> controls/CONTROL.json; exit 1 if any control count differs."""
import json, sys, csv
from pathlib import Path
M = Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295'); C = '20210513_EXPL0_PeTr_SA_T760_1311'
I3 = Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration3/disjoint'); PCF = Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/public347/crossfit_capscan347'); S = Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all')
def lines(p): return set(Path(p).read_text().split())
def cmp(mine, arch):
    a, b = lines(mine), lines(arch)
    return dict(rerun=len(a), archived=len(b), count_equal=len(a) == len(b), set_equal=a == b, only_rerun=len(a - b), only_archived=len(b - a))
res = {}
try:
    res['b_mg_same_spectrum'] = json.loads((M / f'mg_full/out/{C}/COUNTS.json').read_text())
    res['b_mg_same_spectrum']['count_equal'] = res['b_mg_same_spectrum']['control_status'] == 'PASS'
except Exception as e: res['b_mg_same_spectrum'] = dict(error=repr(e), count_equal=False)
try:
    d = M / f'disjoint/outputs/CAPSCAN/{C}'
    j = json.loads((d / 'COMPLETE.json').read_text())
    res['d_mg_odd'] = dict(odd_mzml_identical=j['odd_mzml_identical'], **cmp(d / 'MG_ODD/ACCEPTED.txt', I3 / f'outputs/CAPSCAN/{C}/MG_ODD/ACCEPTED.txt'))
    res['d_mg_full_vs_frozen'] = {k: v for k, v in j.items() if k.startswith(('mg_full_', 'frozen_only'))}
except Exception as e: res['d_mg_odd'] = dict(error=repr(e), count_equal=False)
try:
    o = M / f'crossfit/outputs/CAPSCAN/{C}'
    for arm in ['FL', 'FULL', 'MG']:
        r = cmp(o / f'{arm}_POOLED_ACCEPTED.txt', PCF / f'outputs/CAPSCAN/{C}/{arm}_POOLED_ACCEPTED.txt')
        res[f'c_crossfit_{arm}_pooled'] = r
        res[f'c_crossfit_{arm}_union'] = cmp(o / f'{arm}_UNION_ACCEPTED.txt', PCF / f'outputs/CAPSCAN/{C}/{arm}_UNION_ACCEPTED.txt')
        res[f'c_crossfit_{arm}_union']['diagnostic_only'] = 'union is known to vary by 1-4 peptides run to run (public347 control); not gated'
except Exception as e: res['c_crossfit'] = dict(error=repr(e), count_equal=False)
try:
    sid = f'CAPSCAN__{C}__MG'
    for L in ['5', '9']:
        res[f'e_second_engine_MG_{L}_50'] = cmp(M / f'second_engine/results/{sid}/accepted_{L}_50.txt', S / f'results/{sid}/accepted_{L}_50.txt')
except Exception as e: res['e_second_engine'] = dict(error=repr(e), count_equal=False)
gated = {k: v for k, v in res.items() if not k.endswith('_union') and k != 'd_mg_full_vs_frozen'}
res['CONTROL_STATUS'] = 'PASS' if all(v.get('count_equal') for v in gated.values()) else 'FAIL'
pil = {}
for s in sys.argv[1:]:
    i = next(r['i'] for r in csv.DictReader(open(M / 'disjoint/roster_mg295.tsv'), delimiter='\t') if r['sample'] == s)
    f = dict(db=M / f'db/_db/{s}.metag_hashed_plus_human.fasta', b_counts=M / f'mg_full/out/{s}/COUNTS.json',
             d_complete=M / f'disjoint/outputs/CAPSCAN/{s}/COMPLETE_odd.json', c_complete=M / f'crossfit/outputs/CAPSCAN/{s}/COMPLETE.json',
             e_receipt=M / f'second_engine/results/CAPSCAN__{s}__MG/RECEIPT.json')
    p = {k: v.exists() for k, v in f.items()}
    if p['b_counts']: p['b'] = json.loads(f['b_counts'].read_text())
    if p['d_complete']: p['d'] = json.loads(f['d_complete'].read_text()).get('MG_ODD')
    if p['c_complete']:
        cc = json.loads(f['c_complete'].read_text()); p['c'] = {a: cc[f'{a}_pooled'] for a in ['FL', 'FULL', 'MG']}; p['c']['folds'] = len(cc['folds'])
    if p['e_receipt']: p['e'] = json.loads(f['e_receipt'].read_text())['peptides_q01']
    p['complete'] = all(v for k, v in p.items() if k in f)
    pil[s] = p
res['PILOT'] = pil
json.dump(res, open(M / 'controls/CONTROL.json', 'w'), indent=1, default=str); print(json.dumps(res, indent=1, default=str))
sys.exit(0 if res['CONTROL_STATUS'] == 'PASS' and all(p['complete'] for p in pil.values()) else 1)
