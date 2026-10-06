#!/usr/bin/env python3
"""Second collector (runs after collect_all.py and the participant-15 tasks): CapScan
FastaLake public-sequence arm at 347 acquisitions, reported separately from the paired
282 FL/MG comparison in RESULTS.md (which this script does not touch).

Per acquisition: MSFragger+Percolator peptides (RECEIPT.json, q<=0.01, canonical) and the
Sage count from the task's own source results.sage.tsv (same FASTA and mzML), counted with
the rule collect_all.py uses for healthy stool (rank 1, label 1, peptide_q<=0.01, canonical).
Writes RESULTS_CAPSCAN_FL347.md and results_capscan_fl347.tsv."""
import csv, json, math, pathlib, statistics, sys
RUN = pathlib.Path(sys.argv[1]).resolve()
sys.path.insert(0, str(RUN / 'frozen'))
from peptide import canonical  # noqa: E402
m = json.loads((RUN / 'inputs/manifest.json').read_text())
tasks = [(i, t) for i, t in enumerate(m['tasks']) if t['cohort'] == 'CAPSCAN' and t['arm'] == 'FL_public']
mg = {t['sample'] for t in m['tasks'] if t['cohort'] == 'CAPSCAN' and t['arm'] == 'MG'}
first_new = m.get('appended', {}).get('first_index', len(m['tasks']))
rows = []
for i, t in tasks:
    r = RUN / 'results' / t['id'] / 'RECEIPT.json'
    x = json.loads(r.read_text()) if r.exists() else None
    cfg = json.loads(pathlib.Path(t['source_config']['path']).read_text())
    tsv = [p for p in cfg['output_paths'] if p.endswith('results.sage.tsv')][0]
    sets = {5: set(), 9: set()}
    for q in csv.DictReader(open(tsv), delimiter='\t'):
        if q['rank'] != '1' or q['label'] != '1' or float(q['peptide_q']) > .01:
            continue
        p = canonical(q['peptide'])
        for lo in (5, 9):
            if lo <= len(p) <= 50:
                sets[lo].add(p)
    n5, n9 = (x['peptides_q01']['len5_50'], x['peptides_q01']['len9_50']) if x else (None, None)
    rows.append(dict(index=i, sample=t['sample'], group='participant15_added' if i >= first_new else 'roster331',
                     matched_MG=t['sample'] in mg, status='PASS' if x else 'MISSING_OR_FAILED',
                     msf_5=n5, msf_9=n9, sage_5=len(sets[5]), sage_9=len(sets[9]),
                     ratio_9=n9 / len(sets[9]) if x and sets[9] else None, ratio_5=n5 / len(sets[5]) if x and sets[5] else None))
with open(RUN / 'results_capscan_fl347.tsv', 'w') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)


def med(v):
    """Median of the finite values, NaN if none."""
    v = [a for a in v if a is not None and math.isfinite(a)]
    return statistics.median(v) if v else float('nan')


fmt = lambda v: 'nan' if v is None or (isinstance(v, float) and math.isnan(v)) else (f'{v:.3f}' if isinstance(v, float) else str(v))
L = ['# Second search engine: CapScan FastaLake public-sequence arm at 347 acquisitions', '',
     'FastaLake only; reported separately from the paired FL/MG comparison (282 acquisitions with a matched metagenome, RESULTS.md), which is unchanged. '
     'MSFragger 4.1 + Percolator 3.06.5, rank-1 target peptides at q<=0.01, modifications removed, I/L collapsed. '
     'Sage counts: the same FASTA and mzML, from each task\'s source results.sage.tsv (rank 1, target, peptide_q<=0.01, canonical). '
     'The 331 FASTAs are the A08 production databases; the 16 participant-15 FASTAs come from final/public347/A08_347.', '',
     '| set | PASS / planned | median MSFragger peptides 9-50 | median 5-50 | median Sage peptides 9-50 | median MSFragger/Sage 9-50 | median MSFragger/Sage 5-50 |',
     '|---|---:|---:|---:|---:|---:|---:|']
for label, sel in [('all 347', rows), ('331 roster', [r for r in rows if r['group'] == 'roster331']),
                   ('16 participant-15 added', [r for r in rows if r['group'] == 'participant15_added']),
                   ('282 with matched MG (paired set)', [r for r in rows if r['matched_MG']])]:
    ok = [r for r in sel if r['status'] == 'PASS']
    L.append(f"| {label} | {len(ok)} / {len(sel)} | {fmt(med([r['msf_9'] for r in ok]))} | {fmt(med([r['msf_5'] for r in ok]))} | "
             f"{fmt(med([r['sage_9'] for r in sel]))} | {fmt(med([r['ratio_9'] for r in ok]))} | {fmt(med([r['ratio_5'] for r in ok]))} |")
miss = [r['sample'] for r in rows if r['status'] != 'PASS']
L += ['', f'Missing or failed: {len(miss)}' + (': ' + ', '.join(miss) if miss else ''), '', 'File: results_capscan_fl347.tsv.']
(RUN / 'RESULTS_CAPSCAN_FL347.md').write_text('\n'.join(L) + '\n')
print('\n'.join(L))
