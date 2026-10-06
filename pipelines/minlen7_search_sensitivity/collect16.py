"""Supp. Fig. 4c 'CapScan FL' (FL_public) row, 331 vs 347, and the search count.
read() is a verbatim copy of M09 wave2/code/collect_peptides.py read() (that script runs on import, so it cannot be imported).
Gate: recomputing the first three production CapScan FL_public rows from their archived searches must equal
per_acquisition_wide_claude.tsv exactly (exit 1 otherwise). --selftest runs only the gate.
"""
from pathlib import Path
import sys, json, csv, statistics as st
F = "LAB_ROOT"
sys.path.insert(0, F + '/MM_review_analysis/A10_4_3_4_5_identity_recovery/frozen')
from peptide import canonical
csv.field_size_limit(1 << 30)
RUN = Path(__file__).resolve().parent.parent
M09 = Path(F + "/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7")
def read(p):
 sets={5:set(),7:set(),9:set()};proteins=set();hist={};n=0
 with Path(p).open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   n+=1
   if r['rank']!='1' or r['label']!='1' or float(r['peptide_q'])>.01:continue
   seq=canonical(r['peptide'])
   if not 5<=len(seq)<=50:continue
   for low in sets:
    if len(seq)>=low:sets[low].add(seq)
   if float(r['protein_q'])<=.01:proteins.update(x for x in r['proteins'].split(';') if x and not x.startswith('rev_'))
 for seq in sets[5]:hist[len(seq)]=hist.get(len(seq),0)+1
 return sets,dict(psm_rows=n,native_identified_accessions=len(proteins),canonical_length_hist=hist)
def row(sample, base_tsv, fresh_tsv):
    b, _ = read(base_tsv); n, _ = read(fresh_tsv)
    r = dict(cohort='CAPSCAN', arm='FL_public', sample=sample, baseline5=float(len(b[5])), fresh7=float(len(n[7])), fresh7_nested9=float(len(n[9])),
             nested7=float(len(b[7])), nested9=float(len(b[9])), gained_fresh_nested7=len(n[7] - b[7]), lost_fresh_nested7=len(b[7] - n[7]))
    r['fresh7_vs_base5'] = 100 * (r['fresh7'] / r['baseline5'] - 1); r['fresh7_vs_nested7'] = 100 * (r['fresh7'] / r['nested7'] - 1)
    r['nested7_vs_base5'] = 100 * (r['nested7'] / r['baseline5'] - 1)
    return r
W = [r for r in csv.DictReader(open(M09 / 'wave2/peptide_summary/per_acquisition_wide_claude.tsv'), delimiter='\t') if r['cohort'] == 'CAPSCAN' and r['arm'] == 'FL_public']
assert len(W) == 331
tasks = {r['sample']: r for r in csv.DictReader(open(M09 / 'tasks.tsv'), delimiter='\t') if r['cohort'] == 'CAPSCAN' and r['arm'] == 'FL_public'}
ok = True
for w in W[:3]:
    t = tasks[w['sample']]; r = row(w['sample'], t['baseline_dir'] + '/results.sage.tsv', t['out_dir'] + '/results.sage.tsv')
    for k in ['baseline5', 'fresh7', 'fresh7_nested9', 'nested7', 'nested9', 'gained_fresh_nested7', 'lost_fresh_nested7']:
        ok &= float(w[k]) == float(r[k])
    for k in ['fresh7_vs_base5', 'nested7_vs_base5']:
        ok &= abs(float(w[k]) - r[k]) < 1e-9
print('GATE', 'PASS' if ok else 'FAIL')
if not ok: sys.exit(1)
if '--selftest' in sys.argv: sys.exit(0)
new = []
for t in csv.DictReader(open(RUN / 'tasks.tsv'), delimiter='\t'):
    rec = json.load(open(RUN / 'receipts' / ('receipt_%04d.json' % int(t['index'])))); assert rec['pass']
    new.append(row(t['sample'], t['baseline_dir'] + '/results.sage.tsv', t['out_dir'] + '/results.sage.tsv'))
assert len(new) == 16
with open(RUN / 'per_acquisition_new16.tsv', 'w') as f:
    w = csv.DictWriter(f, fieldnames=list(new[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(new)
old = dict(n=331, filtered_to_7=st.median(float(r['nested7_vs_base5']) for r in W), new_search_7=st.median(float(r['fresh7_vs_base5']) for r in W))
allr = [float(r['nested7_vs_base5']) for r in W] + [r['nested7_vs_base5'] for r in new]
alln = [float(r['fresh7_vs_base5']) for r in W] + [r['fresh7_vs_base5'] for r in new]
res = dict(supp_fig4c_capscan_FL_public=dict(old=old, new=dict(n=347, filtered_to_7=st.median(allr), new_search_7=st.median(alln))),
           search_count=dict(old=1990, after_reruns_347_A2=2006, with_FL_public_16=2022))
json.dump(res, open(RUN / 'RESULTS.json', 'w'), indent=1); print(json.dumps(res))
