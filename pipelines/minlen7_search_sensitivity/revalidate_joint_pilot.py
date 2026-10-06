from pathlib import Path
import json, csv, re, time, os, shutil
from run_search import sha,streamsha,clean
from joint_lfq_audit import audit_joint
O=Path(__file__).resolve().parent.parent
p=O/'receipts/pilot_/088.json'
old=json.loads(p.read_text())
assert old['status']=='FAIL' and old['error']=='''NameError("name 'Path' is not defined")'''
archive=p.with_name('088.failed_missing_import.json')
assert not archive.exists()
shutil.copy2(p,archive)
t=json.loads((O/'tasks.json').read_text())[88]
dest=Path(old['output'])
assert sha(t['source'])==t['source_sha'] and sha(t['config'])==t['config_sha'] and sha(t['binary'])==t['binary_sha']
cfg=json.loads((dest/'config.json').read_text())
expected=clean(json.loads(Path(t['source']).read_text()))
expected['database']['enzyme']['min_len']=7
expected['mzml_paths']=expected['mzml_paths'][:2]
assert clean(cfg)==expected and clean(json.loads((dest/'results.json').read_text()))==expected
assert streamsha(cfg['database']['fasta'])==t['fasta_sha']
assert all(streamsha(k)==v for k,v in old['raw_hashes'].items())
audit=audit_joint(dest,cfg)
accepted=0; files=set()
with (dest/'results.sage.tsv').open() as f:
 for r in csv.DictReader(f,delimiter='	'):
  if r['rank']=='1' and r['label']=='1' and float(r['peptide_q'])<=.01:
   seq=re.sub(r'\[[^\]]*\]|[^A-Z]','',r['peptide'])
   assert 7<=len(seq)<=50
   accepted+=1;files.add(r['filename'])
assert accepted>0 and len(files)==2 and (dest/'results.sage.pin').stat().st_size>0
rec=dict(old)
rec.update(status='PASS',revalidated_utc_epoch=time.time(),revalidation_job=os.getenv('SLURM_JOB_ID'),prior_failed_receipt=str(archive),prior_failed_receipt_sha256=sha(archive),revalidation='Existing successful Sage outputs; fixed missing Path import in audit only; search not rerun',lfq_audit=audit,accepted_psms=accepted,n_acquisitions=2,resolved_settings_match=True,outputs={n:streamsha(dest/n) for n in ['results.json','results.sage.tsv','lfq.tsv','lfq_target_decoy_audit.tsv']})
rec.pop('error');rec.pop('traceback')
tmp=p.with_suffix('.tmp');tmp.write_text(json.dumps(rec,indent=2)+chr(10));tmp.replace(p)
print(json.dumps(rec,indent=2))
