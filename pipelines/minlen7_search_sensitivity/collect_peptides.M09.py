from pathlib import Path
import sys,json,csv,gzip,hashlib,time,traceback
ROOT=Path.cwd().resolve();O=Path(__file__).resolve().parent.parent;R=O.parent
sys.path.insert(0,str(ROOT/'MM_review_analysis/A10_4_3_4_5_identity_recovery/frozen'))
from peptide import canonical
csv.field_size_limit(1<<30)
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
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
def run(t):
 idx=int(t['index']);out=O/'peptide_tasks'/('%04d'%idx)
 if (out/'COMPLETE.json').exists():
  c=json.loads((out/'COMPLETE.json').read_text());assert c['status']=='PASS';return
 out.mkdir(parents=True,exist_ok=False)
 receipt=R/'receipts'/('receipt_%04d.json'%idx)
 try:
  c=json.loads(receipt.read_text());assert c['pass'] is True
  assert sha(t['source_config'])==t['source_sha256']
  b=Path(t['baseline_dir'])/'results.sage.tsv';n=Path(t['out_dir'])/'results.sage.tsv'
  base,bqc=read(b);new,nqc=read(n);assert new[5]==new[7]
  pops={'baseline5':base[5],'nested7':base[7],'nested9':base[9],'fresh7':new[7],'fresh7_nested9':new[9]}
  for name,seqs in pops.items():
   with gzip.open(out/(name+'.txt.gz'),'wt') as h:h.write('\n'.join(sorted(seqs))+'\n')
  summary=dict(index=idx,cohort=t['cohort'],arm=t['arm'],sample=t['sample'],counts={k:len(v) for k,v in pops.items()},fresh_vs_nested7=dict(shared=len(new[7]&base[7]),gained=len(new[7]-base[7]),lost=len(base[7]-new[7])),baseline_qc=bqc,fresh_qc=nqc)
  (out/'SUMMARY.json').write_text(json.dumps(summary,indent=2)+'\n')
  proof=dict(status='PASS',index=idx,source_psms={str(b):sha(b),str(n):sha(n)},receipt_sha=sha(receipt),outputs={p.name:sha(p) for p in out.iterdir() if p.is_file()},gate='rank1 target peptide_q<=0.01; canonical modifications removed and I/L collapsed. Protein count additionally protein_q<=0.01.')
  (out/'COMPLETE.json').write_text(json.dumps(proof,indent=2)+'\n')
 except Exception as e:
  (out/'FAILED.json').write_text(json.dumps(dict(error=repr(e),traceback=traceback.format_exc()),indent=2)+'\n');raise
tasks=list(csv.DictReader((R/'tasks.tsv').open(),delimiter='\t'))
if sys.argv[1]=='pilot':ids=[int(x) for x in sys.argv[2:]]
else:
 block=int(sys.argv[1]);ids=list(range(block*20,min(len(tasks),(block+1)*20)))
for idx in ids:run(tasks[idx]);print('PASS collected',idx,flush=True)
