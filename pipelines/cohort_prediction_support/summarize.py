"""Acquisition support rows from A09 protein_support.tsv, using the record logic of
MM_review_analysis/A09_4_1_prediction_support/frozen/collect.py::support() (loop body copied unchanged).
Usage: summarize.py repro | new"""
import sys, json, csv, io, pathlib, hashlib
import pandas as pd
OUT = pathlib.Path(__file__).resolve().parents[1]; R = OUT / 'a09_run'
A9 = pathlib.Path('LAB_ROOT/MM_review_analysis/A09_4_1_prediction_support')
ORIG = A9 / 'results/acquisition_support_summary.tsv'
def sha(p):return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def completion(p,names=None):
 c=json.loads((p/'COMPLETE.json').read_text());assert c['status']=='PASS'
 for x in c['outputs']:
  f=pathlib.Path(x['path'])
  if names is None or f.name in names:assert sha(f)==x['sha256'],f
 return c
m = json.loads((R / 'inputs/manifest.json').read_text())
mode = sys.argv[1]
tasks = m['tasks'][:3] if mode == 'repro' else m['tasks'][3:]
rows=[];hist=[]
for t in tasks:
  # ---- begin: verbatim from collect.py support() ----
  p=R/'results'/t['cohort']/t['sample'];c=completion(p,{'protein_support.tsv'});x=pd.read_csv(p/'protein_support.tsv',sep='\t');assert len(x)==c['populations']['all_selected']['n'];assert not x.protein.duplicated().any();assert (x.production_prediction_sequences>0).all() and (x.parsimony_prediction_sequences>0).all()
  for pop,xx in [('all_selected',x),('identified_accessions',x[x.identified_accession==1]),('native_protein_q01',x[x.native_protein_q01==1])]:
   record=dict(cohort=t['cohort'],sample=t['sample'],population=pop,n=len(xx),one_top_prediction=int((xx.production_prediction_sequences==1).sum()),one_parsimony_prediction=int((xx.parsimony_prediction_sequences==1).sum()),with_FL_corroborated_prediction=int((xx.predictions_corroborated_by_FL_peptides>0).sum()))
   if t['cohort']=='CAMPI':record.update(with_REF_correct_top=int((xx.same_spectrum_REF_correct_top>0).sum()),with_REF_correct_all=int((xx.same_spectrum_REF_correct_all>0).sum()))
   assert record['n']==c['populations'][pop]['n']
   for key in c['populations'][pop]:assert record[key]==c['populations'][pop][key],(t['sample'],pop,key)
   for key in ['production_prediction_sequences','parsimony_prediction_sequences','production_prediction_spectra','parsimony_prediction_spectra','unique_among_selected_proteins','shared_among_selected_proteins']:
    for n,v in xx[key].value_counts().items():hist.append(dict(cohort=t['cohort'],sample=t['sample'],population=pop,support_definition=key,support_count=int(n),protein_acquisition_entries=int(v)))
    record[key+'_median']=float(xx[key].median()) if len(xx) else None
   for key in list(record):
    if key.startswith(('one_','with_')):record[key+'_fraction']=record[key]/len(xx) if len(xx) else None
   rows.append(record)
  # ---- end verbatim ----
header = open(ORIG).readline().rstrip('\n').split('\t')
d = pd.DataFrame(rows).reindex(columns=header)   # the 710 table has the CAMPI-only REF columns; empty for CapScan
buf = io.StringIO(); d.to_csv(buf, sep='\t', index=False); lines = buf.getvalue().splitlines()[1:]
if mode == 'repro':
    orig = {tuple(l.split('\t')[1:3]): l for l in open(ORIG).read().splitlines()[1:]}
    res = []
    for l in lines:
        k = tuple(l.split('\t')[1:3]); ok = orig[k] == l; res.append(dict(sample=k[0], population=k[1], identical=ok))
        print('REPRO', k, 'identical' if ok else 'DIFF\n new ' + l + '\n old ' + orig[k])
        assert ok
    (OUT / 'checks/support_repro.json').write_text(json.dumps(res, indent=1) + '\n'); print('SUPPORT REPRO PASS')
else:
    assert len(lines) == 48
    (OUT / 'work/acquisition_support_summary_new16.tsv').write_text('\t'.join(header) + '\n' + '\n'.join(lines) + '\n')
    print('SUPPORT NEW16 WRITTEN')
