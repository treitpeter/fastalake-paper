"""M09 fusion control on the CURRENT broad-reservoir CAMPI fusion (campi_fusion_rerun, 375,333,654 entries; manuscript P37).
Wave 2 (Codex) reran the historical 233M series instead; this adds the matching control. Usage: run_fusion.py <0..29>.
Task i: seed = [20260907,20260908][i//15], sample = S01..S15 order of the rerun. Only change: min_len 5 -> 7 (asserted).
Scoring: the frozen run_experiment.read_scores/estimates (imported, unchanged). Also re-scores the ORIGINAL 5-50 results
at floors 5/7/9 and asserts that floors 5 and 9 reproduce the rerun's saved summary exactly."""
import sys, json, copy, subprocess, hashlib, importlib.util, csv, os, time
from pathlib import Path
FR = Path('LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/campi_fusion_rerun')
spec = importlib.util.spec_from_file_location('fx', FR/'run_experiment.py'); fx = importlib.util.module_from_spec(spec); spec.loader.exec_module(fx)
SAGE = FR/'bin/sage'; SAGE_SHA = 'd3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625'
OUTROOT = Path('LAB_ROOT/canonical/publication_data/arms/M09_minlen7_2026-09-24/CAMPI_FUSION_375M')
HERE = Path(__file__).resolve().parent
sha = fx.sha
i = int(sys.argv[1]); seed = fx.SEEDS[i // 15]
samples = sorted(p.name for p in (FR/'runs_v2/CAMPI_MOCK'/str(seed)/'samples').iterdir() if p.is_dir())
assert len(samples) == 15; sample = samples[i % 15]
base = FR/'runs_v2/CAMPI_MOCK'/str(seed)/'samples'/sample
done = json.loads((base/'COMPLETE.json').read_text())
assert sha(base/'5_sage/results.sage.tsv') == done['results_sha256']
assert sha(SAGE) == SAGE_SHA
src = json.loads((base/'search_input/sage.json').read_text())
out = OUTROOT/str(seed)/sample; assert not (out/'results.sage.tsv').exists(); out.mkdir(parents=True, exist_ok=True)
cfg = copy.deepcopy(src); cfg['database']['enzyme']['min_len'] = 7; cfg['output_directory'] = str(out)
a, b = copy.deepcopy(src), copy.deepcopy(cfg)
for d in (a, b): d.pop('output_directory', None); d['database']['enzyme'].pop('min_len')
assert a == b
json.dump(cfg, open(out/'config.json', 'w'), indent=2)
t0 = time.time(); p = subprocess.run([str(SAGE), str(out/'config.json'), '--write-pin'], stdout=open(out/'sage.stdout','w'), stderr=open(out/'sage.stderr','w'))
assert p.returncode == 0, 'sage failed'
res = json.loads((out/'results.json').read_text()); assert res['database']['enzyme']['min_len'] == 7
pairs = {r['target']: r['entrapment'] for r in csv.DictReader(open(base/'search_input/peptide_pairs.tsv'), delimiter='\t')}
bs, bc = fx.read_scores(base/'5_sage/results.sage.tsv'); ns, nc = fx.read_scores(out/'results.sage.tsv')
assert all(len(q) >= 7 for q, v in ns.items() if v <= 0.05), 'accepted peptide shorter than 7'
Q = [.001, .005, .01, .02, .05]
rows = []
for q in Q:
    for lo in (5, 7, 9): rows.append(dict(search='baseline', **fx.estimates(bs, bc, pairs, q, lo)))
    for lo in (7, 9): rows.append(dict(search='fresh7', **fx.estimates(ns, nc, pairs, q, lo)))
# reproduction of the rerun's saved per-acquisition summary (floors 5 and 9)
saved = {(r['q_threshold'], r['min_residues']): r for r in done['summary']}
for r in rows:
    if r['search'] == 'baseline' and (r['q_threshold'], r['min_residues']) in saved:
        s = saved[(r['q_threshold'], r['min_residues'])]
        for k in ('n_canonical', 'paired_numerator', 'n_unpaired_accepted'): assert r[k] == s[k], (k, r[k], s[k])
with open(HERE/'receipts'/f'summary_{i:02d}.tsv', 'w', newline='') as f:
    w = csv.DictWriter(f, ['seed', 'sample'] + list(rows[0]), delimiter='\t'); w.writeheader()
    for r in rows: w.writerow(dict(seed=seed, sample=sample, **r))
json.dump(dict(status='PASS', index=i, seed=seed, sample=sample, elapsed_s=round(time.time()-t0,1), baseline_reproduced=True,
               results_sha256=sha(out/'results.sage.tsv'), source_results_sha256=done['results_sha256'], sage_sha256=SAGE_SHA,
               code_sha256=sha(__file__), interpretation='Fusion-linked diagnostic; not FDR certification.'),
          open(HERE/'receipts'/f'receipt_{i:02d}.json', 'w'), indent=1)
print('PASS', i, seed, sample)
