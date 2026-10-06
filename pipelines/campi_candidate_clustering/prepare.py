from pathlib import Path
import json,hashlib,copy,subprocess,sys
B=Path(__file__).resolve().parent;ROOT=B.parents[3]
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for block in iter(lambda:f.read(8<<20),b''):h.update(block)
 return h.hexdigest()
def rec(p,expected=None):
 p=Path(p);st=p.stat();r={'path':str(p),'bytes':st.st_size,'mtime_ns':st.st_mtime_ns}
 if expected:r['sha256']=expected
 elif st.st_size<150_000_000:r['sha256']=sha(p)
 return r
h4p=ROOT/'Hartout_edits/H04_cluster_retention/attempt_v2/manifest.json';h8p=ROOT/'Hartout_edits/H08_cluster_searches/attempt_v2/manifest.json';h4=json.loads(h4p.read_text());h8=json.loads(h8p.read_text());m=copy.deepcopy(h4);m['new_thresholds']=['0.95','0.90','0.75','0.50'];m['pilot_indices']=[0,5];m['pilot_samples']=[m['tasks'][i]['sample'] for i in m['pilot_indices']];m['sage']=h8['sage'];m['s02_diagnostic']=str(ROOT/'MM_review_analysis/H08_S02_recovery_2026-09-14/diagnostic')
for r in h4['binaries']+h8['binaries']:
 assert sha(r['path'])==r['sha256'];st=Path(r['path']).stat();assert st.st_size==r['bytes'] and st.st_mtime_ns==r['mtime_ns']
m['binaries']=list({r['path']:r for r in h4['binaries']+h8['binaries']}.values())
for t,u in zip(m['tasks'],h8['tasks']):
 assert t['sample']==u['sample'];prior4=Path(u['H04']);prior8=ROOT/'Hartout_edits/H08_cluster_searches/attempt_v2/results'/t['sample']
 if t['sample']=='S02':prior8=ROOT/'MM_review_analysis/H08_S02_recovery_2026-09-14/results/S02'
 c4=json.loads((prior4/'COMPLETE.json').read_text());c8=json.loads((prior8/'COMPLETE.json').read_text());assert c4['status']==c8['status']=='PASS';assert c4['baseline_sequence_set_reproduced'];assert c8['baseline_exact'] or (t['sample']=='S02' and c8['baseline_drift_documented'])
 checks={r['path']:r['sha256'] for r in json.loads((prior4/'SOURCE_HASHES.json').read_text())};checks.update({r['path']:r['sha256'] for r in c4['outputs']})
 for r in t['sources']:
  st=Path(r['path']).stat();assert st.st_size==r['bytes'] and st.st_mtime_ns==r['mtime_ns'];r['sha256']=checks[r['path']]
 t.update(prior_h04=str(prior4),prior_search=str(prior8),config=u['config'],archived_search=u['archived_search'],reference_peptides=u['reference_peptides'],members=u['members'])
 if t['sample'] in m['pilot_samples']:
  for rel in ['cluster_099/representatives.fasta','cluster_099_parsimony/'+t['sample']+'_razor.fasta']:
   p=prior4/rel;t['sources'].append(rec(p,checks[str(p)]))
 t['search_sources']=[]
 for r in u['sources']:
  st=Path(r['path']).stat();assert st.st_size==r['bytes'] and st.st_mtime_ns==r['mtime_ns'];t['search_sources'].append(rec(r['path'],r.get('sha256')))
 t['search_sources'] += [rec(prior4/'COMPLETE.json'),rec(prior8/'COMPLETE.json')]
 cfg=json.loads(Path(t['config']).read_text());assert cfg['version']=='0.14.6';assert cfg.get('report_psms',1)==1,cfg.get('report_psms')
folder=ROOT/'reviews/supplementary_completion_2026-09-13/assignment_exports/CAMPI_FL_REF_MG';inputs=json.loads((folder/'SEARCH_INPUTS.json').read_text());assert len(inputs)==15
m['annotations']=str(folder/'candidate_annotations.tsv.gz');m['annotation_inputs']=[rec(folder/n) for n in ['candidate_annotations.tsv.gz','SEARCH_INPUTS.json','COMPLETE.json']];m['function_fastas']=[rec(r['searched_fasta'],r['fasta_sha256']) for r in inputs]
m['lineage']=[rec(p) for p in [h4p,h8p,ROOT/'MM_review_analysis/H08_S02_recovery_2026-09-14/independent_full/full/COMPLETE.json',Path(m['s02_diagnostic'])/'COMPLETE.json',B/'PROTOCOL.md']]
m['frozen_code']=[rec(p) for p in sorted((B/'frozen').glob('*.py'))]
m['scope']='New95/90/75/50 candidate clustering and searches, all15 CAMPI acquisitions; existing97/99 and fixed baselines retained. Functional evidence uses fixed sequence/KO dictionary; no whole-workflow FDR claim.'
(B/'manifest.json').write_text(json.dumps(m,indent=2)+'\n')
py='python3'
for p in (B/'frozen').glob('*.py'):compile(p.read_text(),str(p),'exec')
for label,lines in {'pilot':[f'{py} {B}/frozen/search.py pilot',f'{py} {B}/frozen/collect.py pilot'],'full':[f'{py} {B}/frozen/search.py "$SLURM_ARRAY_TASK_ID"'],'collect':[f'{py} {B}/frozen/collect.py']}.items():
 (B/(label+'.sh')).write_text('#!/bin/bash\nset -euo pipefail\nexport OMP_NUM_THREADS=2 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 RAYON_NUM_THREADS=2\n'+'\n'.join(lines)+'\n')
print(json.dumps({'tasks':len(m['tasks']),'new_thresholds':m['new_thresholds'],'pilot_samples':m['pilot_samples'],'manifest_sha256':sha(B/'manifest.json'),'syntax':'PASS'},indent=2))
