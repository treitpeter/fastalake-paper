from pathlib import Path
import json,copy,hashlib,functools
O=Path(__file__).resolve().parent.parent;ROOT=Path.cwd().resolve()
sha=functools.lru_cache(None)(lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest())
for n in ['configs','receipts','logs','pilot','outputs']:(O/n).mkdir(exist_ok=True)
assert not (O/'tasks.json').exists(),'Existing prepared wave2 must not be replaced'
sage=ROOT/'MM_review_analysis/A13_4_7_predict_unspecific/bin/sage'
patched=ROOT/'MM_review_analysis/A12_7_quantification/source/target/release/sage'
assert sha(sage)=='d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625'
assert sha(patched)=='21241c22a7e886789e3abfb3f129c299ec06080b8ccb504d54906e17689a3420'
a21=json.loads((ROOT/'MM_review_analysis/A21_7_hstool_quant_comparison/inputs/manifest_revision2.json').read_text())
sources=[]
for arm in a21['arms']:
 if arm['name']=='union_individual':
  for r in arm['records']:sources.append(dict(kind='common',cohort='HSTOOL',sample=r['id'],source=r['config']['path'],source_sha=r['config']['sha256'],binary=str(sage),extra={}))
 if arm['name']=='joint_corrected':
  r=arm['records'][0];joint=dict(kind='joint',cohort='HSTOOL',sample='joint58',source=r['config']['path'],source_sha=r['config']['sha256'],binary=str(patched),extra={})
a11=json.loads((ROOT/'MM_review_analysis/A11_4_2_control_reconciliation/inputs/manifest.json').read_text())
for r in a11['fusion']:
 if r['cohort']=='CAMPI_MOCK':
  p=Path(r['directory'])/'5_sage/results.json'
  sources.append(dict(kind='fusion',cohort='CAMPI',sample=r['sample']+'_seed'+str(r['seed']),source=str(p),source_sha=sha(p),binary=str(sage),extra=r))
assert len(sources)==88
sources.append(joint);tasks=[]
for idx,r in enumerate(sources):
 assert sha(r['source'])==r['source_sha']
 src=json.loads(Path(r['source']).read_text());assert src['database']['enzyme']['min_len']==5
 out=O/'outputs'/r['sample'] if r['kind']=='joint' else O/'outputs'/r['kind']/r['sample']
 assert not out.exists()
 cfg=copy.deepcopy(src)
 for k in ['version','output_paths']:cfg.pop(k,None)
 cfg['database']['enzyme']['min_len']=7;cfg['output_directory']=str(out)
 path=O/'configs'/('%03d.json'%idx);path.write_text(json.dumps(cfg,indent=2)+'\n')
 task=dict(r,index=idx,config=str(path),config_sha=sha(path),output=str(out),binary_sha=sha(r['binary']),n_acquisitions=len(cfg['mzml_paths']),fasta_sha=sha(cfg['database']['fasta']),baseline_dir=str(Path(r['source']).parent))
 assert all(Path(x).is_file() for x in cfg['mzml_paths'])
 tasks.append(task)
assert len(tasks)==89 and tasks[-1]['n_acquisitions']==58
assert {t['sample'] for t in tasks if t['kind']=='common'}==set(a21['samples'])
(O/'tasks.json').write_text(json.dumps(tasks,indent=2)+'\n')
pilots=[max([t for t in tasks if t['kind']==kind],key=lambda t:Path(json.loads(Path(t['config']).read_text())['database']['fasta']).stat().st_size)['index'] for kind in ['common','fusion','joint']]
(O/'pilots.json').write_text(json.dumps(pilots)+'\n')
(O/'PREPARE.json').write_text(json.dumps(dict(status='PASS',tasks=89,common=58,fusion=30,joint=1,pilot_indices=pilots,scope='Search-only min_len 5 to 7. Joint engineering pilot uses first two original acquisitions; full uses all58. No new selection or shuffled database generation.'),indent=2)+'\n')
print((O/'PREPARE.json').read_text())
