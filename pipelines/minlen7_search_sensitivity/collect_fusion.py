from pathlib import Path
import json,sys,hashlib,csv,importlib.util
ROOT=Path.cwd().resolve();O=Path(__file__).resolve().parent.parent
src=ROOT/'MM_review_analysis/A11_4_2_control_reconciliation/frozen/audit_fusion.py'
spec=importlib.util.spec_from_file_location('audit_fusion',src);audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)
idx=int(sys.argv[1]);t=json.loads((O/'tasks.json').read_text())[idx];assert t['kind']=='fusion'
receipt=json.loads((O/'receipts'/('%03d.json'%idx)).read_text());assert receipt['status']=='PASS'
out=O/'fusion_collect'/('%03d'%idx);out.mkdir(parents=True,exist_ok=False)
base=Path(t['extra']['directory']);pairs,reasons,pairing=audit.reconstruct_pairs(base/'search_input/search.fasta')
saved={r['target']:r['entrapment'] for r in csv.DictReader((base/'search_input/peptide_pairs.tsv').open(),delimiter='\t')};assert saved==pairs
old,nold=audit.read_search(base/'5_sage/results.sage.tsv');new,nnew=audit.read_search(Path(t['output'])/'results.sage.tsv')
rows=[]
for label,psm in [('baseline',old),('fresh7',new)]:
 for q in audit.THRESHOLDS:
  for low,high in [(5,50),(7,50),(9,50),(5,6),(7,8)]:
   if label=='fresh7' and low<7:continue
   r=audit.decompose(psm,pairs,reasons,q,low,high);r.update(search=label,sample=t['sample'],seed=t['extra']['seed'],cohort='CAMPI');rows.append(r)
audit.write_tsv(out/'decomposition.tsv',rows)
audit.write_json(out/'COMPLETE.json',dict(status='PASS',index=idx,pairing=pairing,baseline_rows=nold,fresh_rows=nnew,frozen_estimator_sha256=audit.sha(src),inputs={str(p):audit.sha(p) for p in [base/'5_sage/results.sage.tsv',Path(t['output'])/'results.sage.tsv',base/'search_input/search.fasta',base/'search_input/peptide_pairs.tsv']},scope='Frozen artificial-fusion diagnostic, not production workflow FDR; all paired/unresolved/adverse components retained.'))
print('PASS fusion',idx)
