from pathlib import Path
import json,sys,os,subprocess,hashlib,csv,re,time,traceback,copy
O=Path(__file__).resolve().parent.parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def streamsha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def clean(d):
 d=copy.deepcopy(d)
 for k in ['version','output_paths','output_directory']:d.pop(k,None)
 return d
def main():
 idx=int(sys.argv[1]);pilot='--pilot' in sys.argv
 t=json.loads((O/'tasks.json').read_text())[idx]
 dest=O/'pilot'/str(idx) if pilot else Path(t['output'])
 dest.mkdir(parents=True,exist_ok=False)
 receipt=O/'receipts'/('pilot_' if pilot else '')/('%03d.json'%idx);receipt.parent.mkdir(exist_ok=True)
 rec=dict(index=idx,kind=t['kind'],pilot=pilot,status='RUNNING',output=str(dest),job=os.getenv('SLURM_JOB_ID'),started=time.time())
 try:
  assert sha(t['source'])==t['source_sha'];assert sha(t['config'])==t['config_sha'];assert sha(t['binary'])==t['binary_sha']
  src=json.loads(Path(t['source']).read_text());cfg=json.loads(Path(t['config']).read_text())
  assert streamsha(cfg['database']['fasta'])==t['fasta_sha']
  expected=clean(src);expected['database']['enzyme']['min_len']=7
  assert clean(cfg)==expected
  if pilot and t['kind']=='joint':cfg['mzml_paths']=cfg['mzml_paths'][:2];expected['mzml_paths']=expected['mzml_paths'][:2]
  cfg['output_directory']=str(dest);(dest/'config.json').write_text(json.dumps(cfg,indent=2)+'\n')
  rec['raw_hashes']={p:streamsha(p) for p in cfg['mzml_paths']}
  env=os.environ.copy();env['RAYON_NUM_THREADS']=os.getenv('SLURM_CPUS_PER_TASK','8')
  cmd=[t['binary'],str(dest/'config.json'),'--write-pin']
  if t['kind']=='joint':
   env['MM_LFQ_RANK_AUDIT']=str(dest/'ranked_lfq_scores.tsv')
   cmd+=['--batch-size','4','--disable-telemetry-i-dont-want-to-improve-sage']
  rec['command']=cmd
  with (dest/'sage.stdout').open('w') as out,(dest/'sage.stderr').open('w') as err:
   run=subprocess.run(cmd,stdout=out,stderr=err,env=env)
  assert run.returncode==0,('sage_exit',run.returncode)
  result=json.loads((dest/'results.json').read_text());assert clean(result)==expected,'Resolved settings differ'
  accepted=0;files=set()
  with (dest/'results.sage.tsv').open() as f:
   for row in csv.DictReader(f,delimiter='\t'):
    if row['rank']=='1' and row['label']=='1' and float(row['peptide_q'])<=.01:
     seq=re.sub(r'\[[^\]]*\]|[^A-Z]','',row['peptide']);assert 7<=len(seq)<=50
     accepted+=1;files.add(row.get('filename',''))
  assert accepted>0 and (dest/'results.sage.pin').stat().st_size>0
  assert (dest/'lfq.tsv').exists()==bool(cfg['quant']['lfq'])
  if t['kind']=='joint':
   from joint_lfq_audit import audit_joint
   rec['lfq_audit']=audit_joint(dest,cfg)
   assert len(files)==len(cfg['mzml_paths']),(len(files),len(cfg['mzml_paths']))
  rec.update(status='PASS',accepted_psms=accepted,n_acquisitions=len(cfg['mzml_paths']),resolved_settings_match=True,outputs={p.name:streamsha(p) for p in dest.iterdir() if p.name in ['results.json','results.sage.tsv','lfq.tsv','lfq_target_decoy_audit.tsv']})
 except Exception as e:
  rec.update(status='FAIL',error=repr(e),traceback=traceback.format_exc());raise
 finally:
  rec['finished']=time.time();tmp=receipt.with_suffix('.tmp');tmp.write_text(json.dumps(rec,indent=2)+'\n');tmp.replace(receipt)
  print(json.dumps({k:v for k,v in rec.items() if k not in ['raw_hashes','outputs','traceback']}),flush=True)
if __name__=='__main__':main()
