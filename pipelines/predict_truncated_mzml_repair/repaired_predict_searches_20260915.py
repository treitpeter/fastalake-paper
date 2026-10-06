#!/usr/bin/env python3
"""Separate matched tryptic/unspecific reruns on a validated RAW reconversion."""
import copy,hashlib,json,os,re,subprocess,sys
os.environ['RAYON_NUM_THREADS']='8'
from pathlib import Path
MM=Path(__file__).resolve().parents[1];REC=MM/'recovery_2026-09-15'
rm=json.loads((REC/'raw_recovery/manifest.json').read_text());rt=rm['tasks'][int(sys.argv[1])]
source=MM/'A13_4_7_predict_unspecific';sm=json.loads((source/'inputs/manifest.json').read_text());task=sm['tasks'][rt['A13_index']]
gate=REC/'raw_recovery'/rt['sample']/'COMPLETE.json';g=json.loads(gate.read_text());assert g['status']=='PASS' and g['sample']==rt['sample']
base=REC/'search_recovery'/rt['sample'];base.mkdir(parents=True,exist_ok=True)
assert not (base/'COMPLETE.json').exists()
def meta(p):
 p=Path(p).resolve();s=p.stat();h=hashlib.sha256()
 with p.open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns,sha256=h.hexdigest())
raw=meta(g['new_mzml']);assert raw['sha256']==g['new_mzml_sha256']
receipts={}
for arm in ['tryptic','unspecific']:
 run=base/arm
 existing=run/'results'/(rt['sample']+'__'+arm)/'COMPLETE.json'
 if existing.exists():
  c=json.loads(existing.read_text());assert c['status']=='PASS' and c['input_raw']['sha256']==raw['sha256']
  assert not re.search(r'\bERROR\b|unhandled XML|panicked at',(existing.parent/'sage.err').read_text())
  for x in c['outputs']:
   actual=meta(x['path']);assert all(actual[k]==x[k] for k in ['sha256','bytes','mtime_ns'])
  receipts[arm]=meta(existing);continue
 run.mkdir(exist_ok=True);(run/'inputs').mkdir(exist_ok=True);(run/'results').mkdir(exist_ok=True)
 if not (run/'frozen').exists():(run/'frozen').symlink_to(source/'frozen',target_is_directory=True)
 m=copy.deepcopy(sm);t=copy.deepcopy(task);t['arm']=arm;t['id']=rt['sample']+'__'+arm;t['raw']=raw;t['require_nontryptic']=arm=='unspecific'
 cfg=json.loads(Path(task['config']['path']).read_text());cfg['mzml_paths']=[raw['path']];cfg['output_directory']=str(run/'results'/t['id']/'search')
 assert 'output_directory' in json.loads(Path(task['config']['path']).read_text())
 if arm=='tryptic':
  original=json.loads(Path(task['source_config']['path']).read_text());cfg['database']['enzyme']=original['database']['enzyme']
  assert cfg['database']['enzyme']['cleave_at']=='KR' and not cfg['database']['enzyme'].get('semi_enzymatic')
 cfgpath=run/'inputs/config.json';cfgpath.write_text(json.dumps(cfg,indent=2)+'\n');t['config']=meta(cfgpath)
 m['tasks']=[t];m['recovery_lineage']=dict(source_manifest=meta(source/'inputs/manifest.json'),source_reconversion=meta(gate),historical_receipt=rt['original_receipt'],scope='Historical comparison preserved; primary tryptic/unspecific comparison uses these matched repaired inputs')
 (run/'inputs/manifest.json').write_text(json.dumps(m,indent=2)+'\n')
 subprocess.run([sys.executable,'-B',str(source/'frozen/search.py'),str(run),'0'],check=True)
 p=run/'results'/t['id']/'COMPLETE.json';c=json.loads(p.read_text());assert c['status']=='PASS'
 log=p.parent/'sage.err';assert not re.search(r'\bERROR\b|unhandled XML|panicked at',log.read_text()),'Repaired source still triggers parser errors'
 receipts[arm]=meta(p)
summary=dict(worker=meta(__file__),RAYON_NUM_THREADS=8,status='PASS',sample=rt['sample'],source_reconversion=meta(gate),raw=raw,receipts=receipts,job_id=os.environ.get('SLURM_JOB_ID'),scientific_review='PENDING_MATCHED_REPAIRED_INPUT_REVIEW')
(base/'COMPLETE.json').write_text(json.dumps(summary,indent=2)+'\n')
print('MATCHED REPAIRED SEARCHES PASS',rt['sample'],flush=True)
