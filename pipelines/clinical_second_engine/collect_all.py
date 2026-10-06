#!/usr/bin/env python3
"""Collect second-engine receipts and write RESULTS.md (runs afterany on all arrays).

Headline: per cohort, median paired FastaLake/MG peptide-count ratio (MSFragger 4.1 +
Percolator 3.06.5, peptide q<=0.01, canonical I/L-collapsed sequences) at 9-50
residues, with 5-50 alongside. Healthy stool: FL MSFragger/Sage count ratio.
Sage comparison numbers come from the main-benchmark Source Data (Fig. 4 benchmark
tables) and, for healthy stool, from the benchmark Sage results.sage.tsv counted
with the A10 rule (rank 1, label 1, peptide_q<=0.01, canonical)."""
import csv, json, math, os, pathlib, statistics, sys
RUN = pathlib.Path(sys.argv[1]).resolve()
sys.path.insert(0, str(RUN / 'frozen'))
from peptide import canonical  # noqa: E402
SD = pathlib.Path(os.environ['FASTALAKE_DATA_DIR']) / 'review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/figure_sd/source_data/package/FastaLake_Source_Data'
A07 = pathlib.Path(os.environ['FASTALAKE_DATA_DIR']) / 'MM_review_analysis/A07_4_2_second_engine/results/within_MSFragger_comparisons.tsv'
m = json.loads((RUN / 'inputs/manifest.json').read_text())


def med(x):
    """Median of the finite values, NaN if none."""
    x = [v for v in x if v is not None and math.isfinite(v)]
    return statistics.median(x) if x else float('nan')


def rows(p):
    """Read a tab-separated file as a list of dicts."""
    with open(p) as f:
        return list(csv.DictReader(f, delimiter='\t'))


def write(p, rs):
    """Write a list of dicts as a tab-separated file."""
    if not rs:
        p.write_text('');return
    with open(p, 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(rs[0]), delimiter='\t', lineterminator='\n');w.writeheader();w.writerows(rs)


per = [];done = {}
for i, t in enumerate(m['tasks']):
    r = RUN / 'results' / t['id'] / 'RECEIPT.json'
    if r.exists():
        x = json.loads(r.read_text());n5, n9 = x['peptides_q01']['len5_50'], x['peptides_q01']['len9_50']
        done[t['cohort'], t['arm'], t['sample']] = (n5, n9);status = 'PASS'
        res = x['resources']['msfragger']
    else:
        n5 = n9 = None;status = 'MISSING_OR_FAILED';res = {}
    per.append(dict(index=i, id=t['id'], cohort=t['cohort'], arm=t['arm'], sample=t['sample'], status=status, peptides_5_50=n5, peptides_9_50=n9,
                    msfragger_max_rss_kb=res.get('max_rss_kb'), msfragger_elapsed=res.get('elapsed')))
write(RUN / 'results_per_task.tsv', per)

bench = rows(SD / 'Figure4/benchmark_0_per_acquisition.tsv') + rows(SD / 'Figure4/benchmark_1_per_acquisition.tsv')
sage = {(r['cohort'], r['test_arm'], r['reference_arm'], int(r['min_length']), r['sample']): r for r in bench}
pairs = [('PREDICT', 'FL_public', 'MG'), ('CAPSCAN', 'FL_public', 'MG'), ('CAPSCAN', 'FL_cohortMG', 'MG')]
paired = [];summary = []
for cohort, fl, ref in pairs:
    planned = sorted({t['sample'] for t in m['tasks'] if t['cohort'] == cohort and t['arm'] == ref})
    ok = [s for s in planned if (cohort, fl, s) in done and (cohort, ref, s) in done]
    for s in ok:
        a, b = done[cohort, fl, s], done[cohort, ref, s]
        sg9, sg5 = sage.get((cohort, fl, ref, 9, s)), sage.get((cohort, fl, ref, 5, s))
        paired.append(dict(cohort=cohort, fl_arm=fl, ref_arm=ref, sample=s, msf_FL_9=a[1], msf_MG_9=b[1], msf_ratio_9=a[1] / b[1] if b[1] else None,
                           msf_FL_5=a[0], msf_MG_5=b[0], msf_ratio_5=a[0] / b[0] if b[0] else None,
                           sage_ratio_9=float(sg9['count_ratio']) if sg9 and sg9['count_ratio'] not in ('', 'nan') else None,
                           sage_ratio_5=float(sg5['count_ratio']) if sg5 and sg5['count_ratio'] not in ('', 'nan') else None))
    pr = [p for p in paired if p['cohort'] == cohort and p['fl_arm'] == fl]
    allsage = lambda lo: med([float(r['count_ratio']) for k, r in sage.items() if k[:4] == (cohort, fl, ref, lo) and r['count_ratio'] not in ('', 'nan')])
    nsage = len([k for k in sage if k[:4] == (cohort, fl, ref, 9)])
    summary.append(dict(cohort=cohort, comparison=fl + '/' + ref, planned=len(planned), n_paired=len(pr),
                        msf_median_ratio_9_50=med([p['msf_ratio_9'] for p in pr]), msf_median_ratio_5_50=med([p['msf_ratio_5'] for p in pr]),
                        msf_FL_more_9_50=sum(1 for p in pr if p['msf_ratio_9'] and p['msf_ratio_9'] > 1),
                        sage_median_ratio_9_50_same_acqs=med([p['sage_ratio_9'] for p in pr]), sage_median_ratio_5_50_same_acqs=med([p['sage_ratio_5'] for p in pr]),
                        sage_median_ratio_9_50_benchmark=allsage(9), sage_median_ratio_5_50_benchmark=allsage(5), sage_benchmark_n=nsage))
write(RUN / 'results_paired.tsv', paired)

# Healthy stool: FL MSFragger / Sage (benchmark Sage search of the same FASTA and mzML)
hs = []
for t in m['tasks']:
    if t['cohort'] != 'HSTOOL':continue
    cfg = json.loads(pathlib.Path(t['source_config']['path']).read_text())
    tsv = [p for p in cfg['output_paths'] if p.endswith('results.sage.tsv')][0]
    sets = {5: set(), 9: set()}
    for r in csv.DictReader(open(tsv), delimiter='\t'):
        if r['rank'] != '1' or r['label'] != '1' or float(r['peptide_q']) > .01:continue
        p = canonical(r['peptide'])
        for lo in (5, 9):
            if lo <= len(p) <= 50:sets[lo].add(p)
    d = done.get(('HSTOOL', 'FL', t['sample']))
    hs.append(dict(sample=t['sample'], sage_5=len(sets[5]), sage_9=len(sets[9]), msf_5=d[0] if d else None, msf_9=d[1] if d else None,
                   ratio_9=(d[1] / len(sets[9])) if d and sets[9] else None, ratio_5=(d[0] / len(sets[5])) if d and sets[5] else None))
write(RUN / 'results_hstool.tsv', hs)
hsok = [h for h in hs if h['msf_9'] is not None]

campi = rows(A07)
cam = {(r['comparator'], r['min_length']): med([float(x['count_ratio']) for x in campi if x['comparator'] == r['comparator'] and x['min_length'] == r['min_length']]) for r in campi}
f = lambda v: '' if v is None or (isinstance(v, float) and math.isnan(v)) else f'{v:.3f}'
L = ['# Second search engine (MSFragger 4.1 + Percolator 3.06.5) across cohorts', '',
     'Pipeline: frozen CAMPI A07_4_2_second_engine code, unchanged; inputs are the main-benchmark FASTAs and mzML (manifest.tsv). Peptides: rank-1 Percolator target peptides, q<=0.01, modifications removed, I/L collapsed. Ratios are per acquisition, then the median.', '',
     '## Paired FastaLake / MG peptide-count ratio', '',
     '| cohort | comparison | n paired / planned | MSFragger median ratio 9-50 (headline) | MSFragger 5-50 | FL more (9-50) | Sage median ratio 9-50, same acquisitions | Sage 5-50, same acquisitions | Sage 9-50, full benchmark | Sage 5-50, full benchmark (n) |',
     '|---|---|---:|---:|---:|---:|---:|---:|---:|---:|']
for s in summary:
    L.append(f"| {s['cohort']} | {s['comparison']} | {s['n_paired']} / {s['planned']} | {f(s['msf_median_ratio_9_50'])} | {f(s['msf_median_ratio_5_50'])} | {s['msf_FL_more_9_50']}/{s['n_paired']} | {f(s['sage_median_ratio_9_50_same_acqs'])} | {f(s['sage_median_ratio_5_50_same_acqs'])} | {f(s['sage_median_ratio_9_50_benchmark'])} | {f(s['sage_median_ratio_5_50_benchmark'])} ({s['sage_benchmark_n']}) |")
L += ['', 'Fig. 4b/c uses FL_public. CapScan FL_cohortMG (noclust_2026-08-07, assembly-inclusive reservoir) is the arm named in the request; it is reported separately and is not the Fig. 4b arm.', '',
      '## Healthy stool (FastaLake only): MSFragger / Sage peptide-count ratio', '',
      '| n / planned | median MSFragger/Sage 9-50 | median 5-50 | median MSFragger peptides 9-50 | median Sage peptides 9-50 |', '|---:|---:|---:|---:|---:|',
      f"| {len(hsok)} / {len(hs)} | {f(med([h['ratio_9'] for h in hsok]))} | {f(med([h['ratio_5'] for h in hsok]))} | {med([h['msf_9'] for h in hsok]) if hsok else ''} | {med([h['sage_9'] for h in hs]) if hs else ''} |",
      '', 'Sage counts for healthy stool were recounted from the benchmark results.sage.tsv of the same FASTA and mzML (rank 1, target, peptide_q<=0.01, canonical).', '',
      '## CAMPI reference (archived A07, 15 acquisitions)', '',
      f"FL/REF median ratio 5-50 {f(cam.get(('REF','5')))}, 9-50 {f(cam.get(('REF','9')))}; FL/MG 5-50 {f(cam.get(('MG','5')))}, 9-50 {f(cam.get(('MG','9')))}.", '',
      '## Task status', '']
from collections import Counter
for (c, a), n in sorted(Counter((p['cohort'], p['arm']) for p in per if p['status'] == 'PASS').items()):
    tot = sum(1 for p in per if p['cohort'] == c and p['arm'] == a)
    L.append(f'- {c} {a}: {n}/{tot} PASS')
miss = [p['id'] for p in per if p['status'] != 'PASS']
L.append(f'- missing or failed: {len(miss)}' + (' (listed in results_per_task.tsv)' if miss else ''))
L += ['', 'Files: results_per_task.tsv, results_paired.tsv, results_hstool.tsv, results/<task>/RECEIPT.json.']
(RUN / 'RESULTS.md').write_text('\n'.join(L) + '\n')
print('\n'.join(L))
