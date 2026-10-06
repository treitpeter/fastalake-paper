#!/usr/bin/env python3
"""Collect full PREDICT rosters with three explicitly repaired acquisition inputs.

Original results and frozen workers remain untouched. Repaired cohort outputs
are separate, and every reused acquisition retains its original receipt.
"""
import copy,hashlib,json,subprocess,sys
from pathlib import Path
MM=Path(__file__).resolve().parents[1];REC=MM/'recovery_2026-09-15';prefix=sys.argv[1]
source=next(MM.glob(prefix+'_*'));run=REC/'cohort_recovery'/prefix;run.mkdir(parents=True,exist_ok=True)
(run/'inputs').mkdir(exist_ok=True);(run/'results').mkdir(exist_ok=True);(run/'results_revision2').mkdir(exist_ok=True);(run/'results_revision4').mkdir(exist_ok=True)
if not (run/'frozen').exists():(run/'frozen').symlink_to(source/'frozen',target_is_directory=True)
def meta(p):
 p=Path(p).resolve();st=p.stat();h=hashlib.sha256(p.read_bytes()).hexdigest();return dict(path=str(p),sha256=h,bytes=st.st_size,mtime_ns=st.st_mtime_ns)
def save(p,x):p.write_text(json.dumps(x,indent=2)+'\n')
rawm=json.loads((REC/'raw_recovery/manifest_revision2.json').read_text());fix={}
for t in rawm['tasks']:
 p=REC/'task_recovery'/t['sample']/'COMPLETE.json';c=json.loads(p.read_text());assert c['status']=='PASS';fix[t['sample']]=c
manifest='manifest_revision4.json' if prefix=='A28' else 'manifest_revision2.json' if prefix in ['A25','A30'] else 'manifest.json'
m=copy.deepcopy(json.loads((source/'inputs'/manifest).read_text()))
for t in m['tasks']:
 sample=t['sample'];c=fix.get(sample)
 if prefix=='A24' and t['cohort']=='PREDICT':
  p=REC/'cohort_recovery/A20/results/PREDICT'/sample/'COMPLETE.json';assert p.exists();t['routes_receipt']=str(p.resolve())
 if c is None:continue
 if prefix=='A20':
  replacement=json.loads((REC/'task_recovery'/sample/'A20/inputs/manifest.json').read_text())['tasks'][0]
  t.update(replacement)
 elif prefix=='A30':t['receipt']=c['receipts']['A29'][0]['path']
 elif prefix in ['A28','A31']:
  t['receipts']={'available':c['receipts']['A27'][0]['path'],'unspecific':c['receipts']['A27'][1]['path']}
 elif prefix=='A25':
  s=json.loads(Path(c['matched_searches']['path']).read_text());tr=json.loads(Path(s['receipts']['tryptic']['path']).read_text());oo={Path(x['path']).name:x for x in tr['outputs']}
  t.update(receipt=c['receipts']['A27'][1]['path'],unspecific_receipt=s['receipts']['unspecific']['path'],baseline_peptides=oo['peptide_q01.txt'],raw=s['raw'])
m['recovery_lineage']=dict(original_manifest=meta(source/'inputs'/manifest),repaired_acquisitions=[meta(REC/'task_recovery'/s/'COMPLETE.json') for s in fix],scope='Full original PREDICT roster; five repaired input paths, all other frozen acquisitions reused')
save(run/'inputs'/manifest,m)
if prefix=='A20':
 folder=run/'results/PREDICT';folder.mkdir(exist_ok=True)
 for t in m['tasks']:
  if t['cohort']!='PREDICT':continue
  target=Path(fix[t['sample']]['receipts']['A20'][0]['path']).parent if t['sample'] in fix else source/'results/PREDICT'/t['sample']
  assert json.loads((target/'COMPLETE.json').read_text())['status']=='PASS'
  p=folder/t['sample']
  if not p.exists():p.symlink_to(target,target_is_directory=True)
 cm=copy.deepcopy(json.loads((source/'inputs/collector_manifest.json').read_text()));cm['source_manifest']=meta(run/'inputs/manifest.json');save(run/'inputs/collector_manifest.json',cm)
 worker,index='collector.py',3
elif prefix=='A24':worker,index='search.py',3
elif prefix=='A28':worker,index='collect_revision4.py',4
elif prefix=='A30':worker,index='collect_revision2.py',2
elif prefix=='A25':worker,index='search_revision2.py',2
elif prefix=='A31':worker,index='search.py',2
else:raise ValueError(prefix)
subprocess.run([sys.executable,'-B',str(source/'frozen'/worker),str(run),str(index)],check=True)
print('REPAIRED FULL PREDICT COLLECTION PASS',prefix,flush=True)
