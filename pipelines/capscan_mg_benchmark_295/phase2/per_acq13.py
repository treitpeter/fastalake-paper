#!/usr/bin/env python3
"""Phase 2, step 1: per-acquisition same-spectrum counts for the 13 new pairs (+ control 1311), using the
FROZEN Review_v16 Confidence reader (lineage/worker_v2.py read()/selection(), imported unchanged), i.e. the
code that produced confidence_ordinary_all_pairs.tsv and the expected FL/MG sets of the 282.

Arms: FL_public = public347/A08_347/results/<s>/search/results.sage.tsv (Sage, LFQ; razor FASTA)
      MG        = mg295/mg_full/out/<s>/results.sage.tsv
      FL_cohortMG = canonical noclust_2026-08-07 CAPSCAN 5_sage razor_hash-acc/<s>/results.sage.tsv (SFig06 rows)
Control: 1311 must reproduce its production confidence rows (all 36) and same_spectrum rows exactly.
Writes phase2/per_acq13/{confidence_rows.tsv, pairs.tsv, CONTROL.json}.
"""
import csv, importlib.util, json, sys
from pathlib import Path
M = Path(__file__).resolve().parents[1]
FINAL = M.parent
LAB = Path('LAB_ROOT')
CL = Path('LAB_ROOT/review/v15_integration/ready_integration_2026-09-28/source_additions/Review_v16/Confidence/lineage')
spec = importlib.util.spec_from_file_location('w2', CL / 'worker_v2.py'); W = importlib.util.module_from_spec(spec); spec.loader.exec_module(W)
REPO = M / 'phase2/repo'
OUT = M / 'phase2/per_acq13'; OUT.mkdir(parents=True, exist_ok=True)
A08 = FINAL / 'public347/A08_347/results'
NOCL = LAB / 'canonical/publication_data/arms/noclust_2026-08-07/cohorts/CAPSCAN/5_sage/razor_hash-acc'
acq = list(csv.DictReader(open(M / 'inputs/ACQ13.tsv'), delimiter='\t'))
CONTROL = '20210513_EXPL0_PeTr_SA_T760_1311'
todo = [(CONTROL, 'colon', '3', 335)] + [(a['sample'], a['region'], a['participant'], 578 + k) for k, a in enumerate(acq)]

conf_rows, pair_rows = [], []
for s, region, part, idx in todo:
    paths = dict(FL_public=A08 / s / 'search/results.sage.tsv', MG=M / f'mg_full/out/{s}/results.sage.tsv',
                 FL_cohortMG=NOCL / s / 'results.sage.tsv')
    arms = {a: W.read(p)[0] for a, p in paths.items()}
    for q in [.01, .005, .001]:
        for lo, hi in [(5, 50), (5, 8), (9, 50)]:
            for gate in arms['MG']:
                f = W.selection(arms['FL_public'][gate], q, lo, hi); g = W.selection(arms['MG'][gate], q, lo, hi)
                fb = W.selection(arms['FL_public']['q_only'], q, lo, hi); gb = W.selection(arms['MG']['q_only'], q, lo, hi)
                conf_rows.append(dict(index=idx, cohort='CAPSCAN', sample=s, gate=gate, q=q, min_length=lo, max_length=hi, n_fl=len(f), n_mg=len(g),
                                      shared=len(f & g), fl_only=len(f - g), mg_only=len(g - f), removed_fl=len(fb - f), removed_mg=len(gb - g), group=region))
    for test, ref in [('FL_public', 'MG'), ('FL_cohortMG', 'MG')]:
        for lo in (5, 9):
            t = W.selection(arms[test]['q_only'], .01, lo, 50); r = W.selection(arms[ref]['q_only'], .01, lo, 50)
            pair_rows.append(dict(sample=s, group=region, subject_id=part, test_arm=test, reference_arm=ref, min_length=lo,
                                  n_test=len(t), n_reference=len(r), overlap=len(t & r)))
    print('done', s, flush=True)

for name, rows in [('confidence_rows.tsv', conf_rows), ('pairs.tsv', pair_rows)]:
    with open(OUT / name, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)

# control: production confidence rows and SFig06 rows of 1311
prod = [r for r in csv.DictReader(open(REPO / 'source_data/supplementary_methods_6_7/confidence_ordinary_all_pairs.tsv'), delimiter='\t') if r['sample'] == CONTROL]
mine = [r for r in conf_rows if r['sample'] == CONTROL]
key = lambda r: (r['gate'], float(r['q']), int(r['min_length']), int(r['max_length']))
pm = {key(r): r for r in prod}
cmp = [all(int(pm[key(r)][c]) == r[c] for c in ['n_fl', 'n_mg', 'shared', 'fl_only', 'mg_only', 'removed_fl', 'removed_mg']) for r in mine]
sf = [r for r in csv.DictReader(open(REPO / 'source_data/SFig06/same_spectrum_benchmark_1_capscan.tsv'), delimiter='\t') if r['sample'] == CONTROL and r['reference_arm'] == 'MG']
sm = {(r['test_arm'], int(r['min_length'])): r for r in pair_rows if r['sample'] == CONTROL}
cmp2 = [all(int(r[c]) == sm[(r['test_arm'], int(r['min_length']))][c] for c in ['n_test', 'n_reference', 'overlap']) for r in sf]
res = dict(confidence_rows_compared=len(cmp), confidence_rows_equal=sum(cmp), sfig06_rows_compared=len(cmp2), sfig06_rows_equal=sum(cmp2))
res['status'] = 'PASS' if len(cmp) == 36 and all(cmp) and len(cmp2) == 4 and all(cmp2) else 'FAIL'
json.dump(res, open(OUT / 'CONTROL.json', 'w'), indent=1); print(res)
sys.exit(0 if res['status'] == 'PASS' else 1)
