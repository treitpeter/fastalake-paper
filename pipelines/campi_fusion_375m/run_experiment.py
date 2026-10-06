"""Frozen full-reservoir fusion-control experiment; diagnostic, not FDR certification."""
from pathlib import Path
import argparse,collections,csv,hashlib,json,math,os,re,shutil,subprocess,time
ROOT=Path(__file__).resolve().parent
SEP='X'*51
SEEDS=[20260907,20260908]
COHORTS=['CAMPI_MOCK','HSTOOL','PREDICT','CAPSCAN']

def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
def save(p,obj):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True);tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(obj,indent=2)+'\n');tmp.replace(p)
def fasta(p):
 head=None;ss=[]
 with Path(p).open() as f:
  for line in f:
   if line.startswith('>'):
    if head is not None:yield head,''.join(ss)
    head=line[1:].split()[0];ss=[]
   else:ss.append(line.strip())
 if head is not None:yield head,''.join(ss)
def table(p,rows,fields=None):
 rows=list(rows);fields=fields or list(rows[0])
 with Path(p).open('w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(rows)
def config():return json.loads((ROOT/'config.json').read_text())
def verify_code():
 j=json.loads((ROOT/'SOURCE_SHA256.json').read_text())
 for name,digest in j.items():
  if sha(ROOT/name)!=digest:raise RuntimeError('Frozen code changed: '+name)
def binary(key):
 b=config()['binaries'][key]
 if sha(b['path'])!=b['sha256']:raise RuntimeError('Frozen binary changed: '+key)
 return b['path']
def run(args,out):
 out=Path(out);out.mkdir(parents=True,exist_ok=False);start=time.time()
 record={'argv':list(map(str,args)),'start':start}
 save(out/'command.json',record)
 env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',RAYON_NUM_THREADS=str(int(os.environ.get('SLURM_CPUS_PER_TASK','8'))))
 with (out/'stdout.log').open('w') as so,(out/'stderr.log').open('w') as se:r=subprocess.run(record['argv'],stdout=so,stderr=se,env=env)
 record.update(returncode=r.returncode,seconds=time.time()-start);save(out/'command.json',record)
 if r.returncode:raise RuntimeError('Command failed: '+str(out))
def samples(cohort):
 c=config()['cohorts'][cohort]
 if sha(c['sample_list'])!=c['sample_list_sha256']:raise RuntimeError('Frozen sample list changed')
 with Path(c['sample_list']).open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
 assert len(rows)==c['n'];return rows
def extract_args(pool,db,out,threads):
 return [binary('extractor'),'--peptides-dir',pool,'--database',db,'--output',out/'selected.fasta','--output-hits',out/'hits.tsv','--top-percent','0.30','--min-length','9','--max-length','50','--normalize-il','true','-t',str(threads)]
def link_predictions(rows,pool):
 pool.mkdir(parents=True,exist_ok=False)
 for r in rows:
  p=Path(r['denovo_csv']);s=p.stat()
  assert (s.st_size,s.st_mtime_ns)==(int(r['denovo_bytes']),int(r['denovo_mtime_ns']))
  (pool/(r['sample_id']+'_predictions.csv')).symlink_to(p)
def audit_linked(path):
 n=0
 for acc,s in fasta(path):
  k=(len(s)-51)//2
  if len(s)!=2*k+51 or s[k:k+51]!=SEP:raise ValueError('Broken fusion group: '+acc)
  n+=1
 if not n:raise ValueError('Empty linked candidate set')
 return {'groups':n,'both_members_retained_in_every_group':True,'sha256':sha(path)}
def prepare(cohort, chosen_seed=None):
 c=config()['cohorts'][cohort];rows=samples(cohort);threads=int(os.environ.get('SLURM_CPUS_PER_TASK','16'))
 for seed in ([chosen_seed] if chosen_seed is not None else SEEDS):
  base=ROOT/'runs_v2'/cohort/str(seed);base.mkdir(parents=True,exist_ok=True)
  if (base/'PREPARED.json').exists():raise FileExistsError('Already prepared; use a new run')
  fused=base/'0_fused_reservoir.fasta';original=Path(c['reservoir']);before=original.stat()
  run([binary('builder'),original,fused,seed,c['reservoir_sequences']],base/'0_build')
  after=original.stat();assert (before.st_size,before.st_mtime_ns)==(after.st_size,after.st_mtime_ns)
  pool=base/'predictions';link_predictions(rows,pool)
  out=base/'1_evidence';run(extract_args(pool,fused,out,threads),out)
  audit=audit_linked(out/'selected.fasta')
  stats=json.loads((out/'selected.stats.json').read_text());assert stats['proteins_searched']==c['reservoir_sequences']
  save(base/'PREPARED.json',dict(cohort=cohort,seed=seed,original=str(original),original_bytes=before.st_size,original_mtime_ns=before.st_mtime_ns,evidence=audit,stats=stats))
  print(cohort,seed,'cohort preparation complete',flush=True)

def digest_intervals(seq,missed=2):
 cuts=[0]+[i+1 for i,c in enumerate(seq) if c in 'KR' and (i+1==len(seq) or seq[i+1]!='P')]
 if cuts[-1]!=len(seq):cuts.append(len(seq))
 for i,start in enumerate(cuts[:-1]):
  for end in cuts[i+1:i+missed+2]:
   if 5<=end-start<=50:yield start,end
   if start==0 and seq.startswith('M') and 5<=end-1<=50:yield 1,end

def split_and_pairs(linked,out):
 te=collections.defaultdict(set);et=collections.defaultdict(set);counts=collections.Counter();pairs_path=out/'peptide_pairs.tsv'
 with (out/'search.fasta').open('w') as f:
  for acc,s in fasta(linked):
   k=(len(s)-51)//2;assert len(s)==2*k+51 and s[k:k+51]==SEP
   t=s[:k];e=s[k+51:];tn=t.replace('I','L');assert len(t)==len(e)
   assert list(digest_intervals(t))==list(digest_intervals(e))
   f.write('>T_'+acc+'\n'+t+'\n>E_'+acc+'\n'+e+'\n');counts['protein_pairs']+=1
   for a,b in digest_intervals(t):
    p=tn[a:b];q=e[a:b]
    if not re.fullmatch('[ACDEFGHKLMNPQRSTVWY]+',p+q):continue
    assert sorted(p)==sorted(q)
    te[p].add(q);et[q].add(p)
 shared=set(te)&set(et);valid={t:next(iter(es)) for t,es in te.items() if len(es)==1 and t not in shared and next(iter(es)) not in shared and len(et[next(iter(es))])==1}
 table(pairs_path,({'target':t,'entrapment':e} for t,e in sorted(valid.items())),['target','entrapment'])
 counts.update(unique_target_peptides=len(te),unique_entrapment_peptides=len(et),shared_sequences=len(shared),valid_disjoint_pairs=len(valid),unpaired_target_peptides=len(te)-len(valid),unpaired_entrapment_peptides=len(et)-len(valid))
 save(out/'PAIRING_AUDIT.json',dict(counts,globally_bijective_over_full_reservoir=False,paired_subset_bijective=True,global_absence_from_specimen_established=False))
 return valid,counts

def canonical(p):
 p=re.sub(r'\[[^\]]*\]','',p).replace('I','L')
 if not re.fullmatch('[A-Z]+',p):raise ValueError('Unknown peptide syntax '+p)
 return p

def read_scores(path):
 scores={};classes=collections.defaultdict(set)
 with Path(path).open() as f:
  for r in csv.DictReader(f,delimiter='\t'):
   if r['label']!='1':continue
   p=canonical(r['peptide']);q=float(r['peptide_q'])
   if not math.isfinite(q):raise ValueError('Nonfinite q')
   scores[p]=min(q,scores.get(p,math.inf))
   for acc in re.split('[;,]',r['proteins']):
    if acc.startswith('T_'):classes[p].add('target')
    elif acc.startswith('E_'):classes[p].add('entrapment')
    elif acc.startswith('rev_'):continue
    else:raise ValueError('Unknown protein class '+acc)
 return scores,classes

def estimates(scores,classes,pairs,threshold,floor):
 accepted={p for p,q in scores.items() if q<=threshold and floor<=len(p)<=50}
 valid={t:e for t,e in pairs.items() if len(t)>=floor}
 paired=set(valid)|set(valid.values());unknown=accepted-paired
 ne=sum(e in accepted for e in valid.values());e_only=sum(e in accepted and t not in accepted for t,e in valid.items())
 e_better=sum(e in accepted and t in accepted and scores[e]<=scores[t] for t,e in valid.items())
 n=len(accepted);numerator=ne+e_only+2*e_better
 return dict(q_threshold=threshold,min_residues=floor,n_canonical=n,n_target_only=sum(classes[p]=={'target'} for p in accepted),n_entrapment_only=sum(classes[p]=={'entrapment'} for p in accepted),n_shared=sum(len(classes[p])>1 for p in accepted),n_unpaired_accepted=len(unknown),n_paired_entrapment=ne,paired_numerator=numerator,paired_upper_diagnostic_with_unpaired_worst_case=(numerator+len(unknown))/n if n else None,combined_upper_diagnostic_with_unpaired_worst_case=(2*ne+len(unknown))/n if n else None,paired_subset_fraction=(n-len(unknown))/n if n else None,whole_procedure_fdr_control_established=False)

def execute_sample(index, chosen_seed=None):
 tasks=config()['tasks'];task=tasks[index];cohort=task['cohort'];r=samples(cohort)[task['index']];threads=int(os.environ.get('SLURM_CPUS_PER_TASK','8'))
 for seed in ([chosen_seed] if chosen_seed is not None else SEEDS):
  base=ROOT/'runs_v2'/cohort/str(seed);assert (base/'PREPARED.json').is_file();out=base/'samples'/r['sample_id']
  if (out/'COMPLETE.json').exists():
   done=json.loads((out/'COMPLETE.json').read_text());assert sha(out/'5_sage/results.sage.tsv')==done['results_sha256'];continue
  out.mkdir(parents=True,exist_ok=False);pool=out/'predictions';link_predictions([r],pool)
  selection=out/'3_selection';run(extract_args(pool,base/'1_evidence/selected.fasta',selection,threads),selection)
  selection_audit=audit_linked(selection/'selected.fasta')
  razor=out/'4_razor';run([binary('parsimony'),'--strategy','razor','--tiebreak','hash-acc','--fasta',selection/'selected.fasta','--peptides',r['denovo_csv'],'--output-dir',razor,'--sample',r['sample_id'],'--min-length','9','--max-length','50','--normalize-il','-t',threads],razor)
  linked=razor/(r['sample_id']+'_razor.fasta');razor_audit=audit_linked(linked)
  bundle=out/'search_input';bundle.mkdir();pairs,pair_counts=split_and_pairs(linked,bundle)
  sage=out/'5_sage';cfg=json.loads(Path(r['sage_config']).read_text());cfg['database']['fasta']=str(bundle/'search.fasta');cfg['output_directory']=str(sage);cfg['mzml_paths']=[r['mzml']]
  save(bundle/'sage.json',cfg);run([binary('sage'),bundle/'sage.json','--write-pin'],sage)
  scores,classes=read_scores(sage/'results.sage.tsv')
  rows=[estimates(scores,classes,pairs,q,lo) for q in [.001,.005,.01,.02,.05] for lo in [5,9]]
  table(out/'summary.tsv',rows)
  table(out/'accepted_entrapment_queries.tsv',({'peptide':p,'peptide_q':scores[p],'assignment':';'.join(sorted(classes[p]))} for p in sorted(scores) if scores[p]<=.05 and 'entrapment' in classes[p]),['peptide','peptide_q','assignment'])
  save(out/'COMPLETE.json',dict(cohort=cohort,sample=r['sample_id'],seed=seed,selection=selection_audit,razor=razor_audit,pairing=dict(pair_counts),results_sha256=sha(sage/'results.sage.tsv'),summary=rows,interpretation='Fusion-linked full-route diagnostic. Equal-chance and full-reservoir unique pairing remain unproven; not a certificate of FDR control.'))
  print(cohort,r['sample_id'],seed,'complete',flush=True)

def collect():
 rows=[];missing=[]
 for task in config()['tasks']:
  cohort=task['cohort'];r=samples(cohort)[task['index']]
  for seed in SEEDS:
   out=ROOT/'runs_v2'/cohort/str(seed)/'samples'/r['sample_id']
   if not (out/'COMPLETE.json').is_file():missing.append([cohort,r['sample_id'],seed]);continue
   d=json.loads((out/'COMPLETE.json').read_text())
   for row in d['summary']:rows.append(dict(cohort=cohort,sample=r['sample_id'],seed=seed,**row))
 dest=ROOT/'results';dest.mkdir(exist_ok=True)
 if rows:table(dest/'per_acquisition.tsv',rows)
 save(dest/'STATUS.json',dict(expected_searches=2*len(config()['tasks']),completed_searches=len(rows)//10,missing=missing,whole_procedure_fdr_control_established=False))
 print('Completed searches',len(rows)//10,'Missing',len(missing),flush=True)
 if missing:raise RuntimeError('Incomplete experiment; partial tables retained')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('stage',choices=['prepare','sample','collect']);p.add_argument('--cohort');p.add_argument('--index',type=int);p.add_argument('--seed',type=int,choices=SEEDS);a=p.parse_args();verify_code()
 if a.stage=='prepare':prepare(a.cohort,a.seed)
 elif a.stage=='sample':execute_sample(a.index,a.seed)
 else:collect()
