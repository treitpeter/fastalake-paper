from pathlib import Path
import csv,json,hashlib,sys,os,shutil
from lib import sha,meta,dump,table,command
H=Path(__file__).resolve().parent
M=json.loads((H/'search_manifest.json').read_text())
sys.path.insert(0,str(H.parent/'confidence'))
from peptide import canonical
def run(i):
 t=M['tasks'][i];out=H/'end_to_end'/str(i);out.mkdir(parents=True,exist_ok=False)
 for b in M['binaries'].values():assert sha(b['path'])==b['sha256']
 assert sha(t['predictions'])==t['predictions_sha256']
 if t['arm']=='terminal':
  assert json.loads((H/'COMPLETE.json').read_text())['status']=='PASS'
 pred=out/'predictions';pred.mkdir();(pred/'sample_predictions.csv').symlink_to(t['predictions'])
 selection=out/'selection';selection.mkdir()
 threads=int(os.environ.get('SLURM_CPUS_PER_TASK','8'))
 command([M['binaries']['extractor']['path'],'--peptides-dir',pred,'--database',t['candidate_union'],'--output',selection/'selected.fasta','--output-hits',selection/'hits.tsv','--top-percent','0.30','--min-length','9','--max-length','50','--normalize-il','true','-t',str(threads)],out/'extract')
 stats=json.loads((selection/'selected.stats.json').read_text());n=stats['total_output_proteins']
 if n==0:
  dump(out/'COMPLETE.json',dict(status='PASS',index=i,arm=t['arm'],seed=t['seed'],selected_proteins=0,search_status='NO_SELECTED_PROTEINS',accepted_counts=[],scope='no search possible; not a zero-error certificate'));return
 razor=out/'razor';razor.mkdir()
 command([M['binaries']['parsimony']['path'],'--strategy','razor','--tiebreak','hash-acc','--fasta',selection/'selected.fasta','--peptides',t['predictions'],'--output-dir',razor,'--sample','sample','--min-length','9','--max-length','50','--normalize-il','-t',str(threads)],out/'parsimony')
 fasta=razor/'sample_razor.fasta';assert fasta.exists()
 cfg=json.loads(Path(M['sage_config']).read_text());cfg['database']['fasta']=str(fasta);cfg['output_directory']=str(out/'sage')
 dump(out/'sage_config.json',cfg)
 command([M['binaries']['sage']['path'],out/'sage_config.json','--write-pin'],out/'search')
 raw=out/'sage/results.sage.tsv';sets={};bad=0
 with raw.open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['rank']!='1' or r['label']!='1':continue
   import math
   q=float(r['peptide_q']);pep=float(r['posterior_error']);assert math.isfinite(q) and math.isfinite(pep)
   p=canonical(r['peptide'])
   for gate in ['q_only','pep01']:
    if gate=='pep01' and pep>-2:continue
    for cutoff in [.01,.005,.001]:
     if q>cutoff:continue
     for lo,hi in [(5,50),(5,8),(9,50)]:
      if lo<=len(p)<=hi:sets.setdefault((gate,cutoff,lo,hi),set()).add(p)
 rows=[]
 for gate in ['q_only','pep01']:
  for q in [.01,.005,.001]:
   for lo,hi in [(5,50),(5,8),(9,50)]:
    seq=sets.get((gate,q,lo,hi),set());rows.append(dict(gate=gate,q=q,min_length=lo,max_length=hi,n_peptides=len(seq)))
 table(out/'accepted_counts.tsv',rows)
 table(out/'accepted_sequences.tsv',[dict(gate=g,q=q,min_length=lo,max_length=hi,peptide=p) for (g,q,lo,hi),pp in sets.items() for p in sorted(pp)],['gate','q','min_length','max_length','peptide'])
 nfinal=sum(1 for l in fasta.open() if l.startswith('>'))
 dump(out/'COMPLETE.json',dict(status='PASS',index=i,arm=t['arm'],seed=t['seed'],selected_proteins=n,parsimony_proteins=nfinal,search_status='COMPLETED',accepted_counts=rows,results_sha256=sha(raw),fasta_sha256=sha(fasta),code_sha256=sha(__file__),manifest_sha256=sha(H/'search_manifest.json'),job_id=os.getenv('SLURM_JOB_ID')))
 print('PASS',i,t['arm'],t['seed'],flush=True)
if __name__=='__main__':run(int(sys.argv[1]))

