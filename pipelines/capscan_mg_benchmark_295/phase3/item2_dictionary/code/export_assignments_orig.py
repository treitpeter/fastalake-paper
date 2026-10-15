"""Export frozen assignments and exact local sequence alternatives. No searches/grouping."""
from pathlib import Path
from collections import defaultdict, Counter
from contextlib import ExitStack
import argparse, csv, gzip, hashlib, io, json, re, shutil, sys, time
import ahocorasick
import pandas as pd

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
BASE=ROOT/'simplified_grouping_results/manuscript_cohorts_2026-09-10_2200'
OUT=HERE/'assignment_exports'
MOD=re.compile(r'\[[^\]]*\]')
def canon(p):return MOD.sub('',p).replace('I','L')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(1048576),b''):h.update(b)
 return h.hexdigest()
def rows(p):
 op=gzip.open if str(p).endswith('.gz') else open
 with op(p,'rt') as f:yield from csv.DictReader(f,delimiter='\t')
def fasta(p):
 header=None;buf=[]
 with Path(p).open() as f:
  for line in f:
   if line.startswith('>'):
    if header is not None:yield header,''.join(buf)
    header=line[1:].strip();buf=[]
   else:buf.append(line.strip().upper())
 if header is not None:yield header,''.join(buf)
def write_json(p,x):Path(p).write_text(json.dumps(x,indent=2)+'\n')
def tasks():
 m=json.loads((BASE/'manifest.json').read_text()); ans={}
 for t in m['tasks']:
  t=t.copy();t['direct']=str(BASE/'cohorts'/t['id']/'direct');t['population']='full_cohort';t['manifest']=str(BASE/'manifest.json')
  ans[t['id']]=t
 for cohort,folder in [('PREDICT','PREDICT_three_arm_2026-09-11_immediate'),('CAPSCAN','CAPSCAN_two_arm_282_2026-09-11_immediate')]:
  b=ROOT/'simplified_grouping_results'/folder;mm=json.loads((b/'manifest.json').read_text())
  for t in mm['tasks']:
   tt=ans[cohort+'_noclust'].copy();tt.update(t);tt.update(id=f"{cohort}_{t['acquisitions']}_{t['id']}",arm=t['id'],population='matched_comparison',direct=t['output'],manifest=str(b/'manifest.json'),samples_file=mm['samples_file'])
   ans[tt['id']]=tt
 return ans

class Annotation:
 def __init__(self,task):
  self.source={};self.byhash=defaultdict(set);self.species={};self.lineage={};self.reference=defaultdict(set)
  sources={h.split()[0]:hashlib.sha256(s.encode()).hexdigest() for h,s in fasta(task['annotation_fasta'])}
  for r in rows(task['annotations']):
   a=r['accession'];ko=tuple(sorted(set(re.findall(r'K\d{5}',r.get('KEGG_ko','')))))
   if a in sources:
    self.source[a]=(sources[a],ko);self.byhash[sources[a]].add(ko)
  del sources
  genomes=list(rows(task['genomes']));fallback={r['Species_rep']:r['Lineage'] for r in genomes if r['Lineage']}
  for r in genomes:
   lineage=r['Lineage'] or fallback.get(r['Species_rep'],'')
   self.lineage[r['Genome']]=lineage
  if task.get('campi_reference'):
   names={r['campi_os_string']:r['member'] for r in rows(task['campi_members'])}
   for h,s in fasta(task['campi_reference']):
    osname=re.split(r' (?:GN|PE|SV|OX)=',h.split('OS=',1)[-1])[0].strip() if 'OS=' in h else ''
    if osname in names:self.reference[hashlib.sha256(s.encode()).hexdigest()].add(names[osname])
 def get(self,acc,seq):
  h=hashlib.sha256(seq.encode()).hexdigest();original=self.source.get(acc)
  kos=original[1] if original and original[0]==h else tuple(sorted({k for x in self.byhash.get(h,()) for k in x}))
  origin='same_accession_exact_sequence' if original and original[0]==h else 'exact_sequence_projection' if h in self.byhash else 'unannotated'
  ref=self.reference.get(h,set());mg=re.search(r'MGYG\d+',acc)
  lineage=self.lineage.get(mg[0],'') if mg else ''
  ranks={k:tuple(x[3:] for x in lineage.split(';') if x.startswith(k+'__') and len(x)>3) for k in ['d','p','c','o','f','g','s']}
  basis='UHGG_lineage' if lineage else 'unclassified'
  if ref:
   ranks['s']=tuple(sorted(ref));ranks['g']=tuple(sorted({s.split()[0] for s in ref}));basis='CAMPI_reference_full_sequence'
  elif '_HUMAN' in acc:
   ranks['s']=('Homo sapiens',);ranks['g']=('Homo',);ranks['d']=('Eukaryota',);basis='human_reference_accession'
  return (h,origin,tuple(kos),ranks['s'],ranks['g'],basis,lineage)

def run(task):
 started=time.time();d=Path(task['direct']);out=OUT/task['id'];out.mkdir(parents=True,exist_ok=True)
 if (out/'COMPLETE.json').exists():
  print(task['id'],'already complete',flush=True);return
 write_json(out/'TASK.json',{k:v for k,v in task.items() if k not in ['search_inputs','static_inputs']})
 s=json.loads((d/'summary.json').read_text());names=s['acquisition_names'];assert len(names)==task['acquisitions']
 original_registry=json.loads((d/'input_manifest.json').read_text())
 registry={str(Path(x['path']).resolve()):x['sha256'] for x in original_registry}
 receipts=[]
 def checked(p):
  p=Path(p);h=sha(p);expected=registry.get(str(p.resolve()));assert expected==h,(p,expected,h)
  receipts.append(dict(path=str(p),bytes=p.stat().st_size,sha256=h));return h
 mapping={r['canonical_peptide']:r['study_group'] for r in rows(d/'peptide_to_study_group.tsv')}
 group_rows=list(rows(d/'study_groups.tsv'));counts=Counter(mapping.values())
 assert len(mapping)==s['pooled_canonical_peptides'] and len(group_rows)==s['inferred_study_groups']
 assert set(counts)=={r['study_group'] for r in group_rows}
 assert all(counts[r['study_group']]==int(r['assigned_peptides']) for r in group_rows)
 annotation=Annotation(task);print(task['id'],'annotation loaded',flush=True)
 gstats=defaultdict(lambda:[set(),set(),set(),False,False,0])
 candidate_cache={};source_lists=[];summaries=[];pooled=set();rep_sequences={};seen_candidate=set()
 with ExitStack() as stack:
  def writer(name,fields):
   f=stack.enter_context(gzip.open(out/name,'wt',compresslevel=1,newline=''))
   w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(fields);return w
  pw=None
  cw=writer('peptide_compatibility.tsv.gz',['sample','canonical_peptide','study_group','n_psms','n_reported_proteins','reported_proteins','n_sequence_compatible_proteins','compatible_proteins','species_candidates','genus_candidates','KO_candidates','n_KO_annotation_sets','n_unclassified_species_candidates','n_missing_KO_candidates','species_conflict','genus_conflict','KO_set_difference','KO_sets_disjoint','representative_in_local_fasta'])
  aw=writer('candidate_annotations.tsv.gz',['accession','sequence_sha256','function_annotation_origin','KEGG_ko','species_candidates','genus_candidates','taxonomy_basis','UHGG_lineage'])
  rw=writer('representative_sequences.tsv.gz',['study_group','sequence_sha256','sequence'])
  for index,name in enumerate(names,1):
   folder=Path(task['search_root'])/name;config=folder/task['config_name'];checked(config)
   cfg=json.loads(config.read_text());fp=Path(cfg['database']['fasta']);fhash=checked(fp)
   sp=folder/'results.sage.tsv';raw=sp.read_bytes();h=hashlib.sha256(raw).hexdigest();assert registry[str(sp.resolve())]==h
   receipts.append(dict(path=str(sp),bytes=len(raw),sha256=h));df=pd.read_csv(io.BytesIO(raw),sep='\t',dtype=str,keep_default_na=False);del raw
   q=pd.to_numeric(df.peptide_q,errors='coerce');accepted=df.loc[df.label.eq('1') & q.between(0,.01)].copy();del df,q
   accepted.insert(0,'sample',name);accepted['canonical_peptide']=accepted.peptide.map(canon);accepted['study_group']=accepted.canonical_peptide.map(mapping)
   assert accepted.study_group.notna().all()
   if pw is None:pw=writer('accepted_psms.tsv.gz',accepted.columns.tolist())
   pw.writerows(accepted.itertuples(index=False,name=None))
   reported=defaultdict(set);npsms=Counter(accepted.canonical_peptide);peps=sorted(npsms);pooled.update(peps)
   for p,field in zip(accepted.canonical_peptide,accepted.proteins):reported[p].update(a for a in field.split(';') if a and not a.startswith(('rev_','DECOY_')))
   automaton=ahocorasick.Automaton()
   for p in peps:automaton.add_word(p,p)
   if peps:automaton.make_automaton()
   compatible=defaultdict(list);local=set();npairs=0
   for header,seq in fasta(fp):
    a=header.split()[0]
    if a.startswith(('rev_','DECOY_')):continue
    assert a not in local,(name,a,'duplicate target accession');local.add(a)
    if a in counts:
     hh=hashlib.sha256(seq.encode()).hexdigest()
     if a not in rep_sequences:rep_sequences[a]=hh;rw.writerow([a,hh,seq])
     else:assert rep_sequences[a]==hh
    hits=set(p for _,p in automaton.iter(seq.replace('I','L'))) if peps else set()
    if not hits:continue
    if a not in candidate_cache:
     candidate_cache[a]=annotation.get(a,seq)
     ah,origin,kos,species,genus,basis,lineage=candidate_cache[a]
     aw.writerow([a,ah,origin,';'.join(kos),';'.join(species),';'.join(genus),basis,lineage])
    else:assert candidate_cache[a][0]==hashlib.sha256(seq.encode()).hexdigest(),('accession sequence changed',a)
    for p in hits:compatible[p].append(a)
   assert all(p in compatible for p in peps),(name,'accepted peptide absent from searched target FASTA')
   stats=Counter()
   for p in peps:
    alts=sorted(compatible[p]);aset=set(alts);assert reported[p] <= aset,(name,p,'reported protein not sequence compatible')
    sg=mapping[p];ss=set();gg=set();ks=set();miss_s=miss_k=0
    for a in alts:
     info=candidate_cache[a];ss.update(info[3]);gg.update(info[4])
     if info[2]:ks.add(info[2])
     else:miss_k+=1
     if not info[3]:miss_s+=1
    ko_union=sorted({k for x in ks for k in x});disjoint=len(ks)>1 and not set.intersection(*(set(x) for x in ks))
    sc=len(ss)>1;gc=len(gg)>1;kc=len(ks)>1
    cw.writerow([name,p,sg,npsms[p],len(reported[p]),';'.join(sorted(reported[p])),len(alts),';'.join(alts),';'.join(sorted(ss)),';'.join(sorted(gg)),';'.join(ko_union),len(ks),miss_s,miss_k,int(sc),int(gc),int(kc),int(disjoint),int(sg in local)])
    g=gstats[sg];g[0].update(ss);g[1].update(gg);g[2].update(ks);g[3]|=bool(miss_s);g[4]|=bool(miss_k);g[5]+=1
    stats.update(dict(peptides=1,compatibility_pairs=len(alts),species_conflicting_peptides=int(sc),genus_conflicting_peptides=int(gc),KO_different_set_peptides=int(kc),KO_disjoint_set_peptides=int(disjoint),species_fully_unclassified_peptides=int(not ss),species_partly_unclassified_peptides=int(bool(ss) and bool(miss_s)),KO_fully_unannotated_peptides=int(not ks),multiple_candidate_peptides=int(len(alts)>1),reported_subset_expanded_peptides=int(len(alts)>len(reported[p]))))
   summaries.append(dict(sample=name,accepted_psms=len(accepted),**stats))
   source_lists.append(dict(sample=name,search_config=str(config),search_config_sha256=receipts[-3]['sha256'],searched_fasta=str(fp),fasta_sha256=fhash,sage_results=str(sp),sage_sha256=h,raw_spectra=cfg.get('mzml_paths',[])))
   if index%10==0 or index==len(names):
    write_json(out/'PROGRESS.json',dict(completed_acquisitions=index,total_acquisitions=len(names),elapsed_seconds=time.time()-started))
    print(task['id'],index,'/',len(names),'peptides',len(peps),'elapsed',round(time.time()-started),flush=True)
  gw=writer('group_candidate_conflicts.tsv.gz',['study_group','assigned_peptides','observed_support_peptides','species_candidates','genus_candidates','KO_candidates','n_KO_annotation_sets','species_conflict','genus_conflict','KO_set_difference','KO_sets_disjoint','any_missing_species_annotation','any_missing_KO_annotation','peptide_acquisition_observations'])
  for r in group_rows:
   a=r['study_group'];ss,gg,ks,ms,mk,n=gstats[a]
   gw.writerow([a,r['assigned_peptides'],r['observed_support_peptides'],';'.join(sorted(ss)),';'.join(sorted(gg)),';'.join(sorted({k for x in ks for k in x})),len(ks),int(len(ss)>1),int(len(gg)>1),int(len(ks)>1),int(len(ks)>1 and not set.intersection(*(set(x) for x in ks))),int(ms),int(mk),n])
 assert pooled==set(mapping) and set(rep_sequences)==set(counts)
 for filename in ['peptide_to_study_group.tsv','study_groups.tsv','acquisition_summary.tsv','study_group_matrix.tsv','study_group_long.tsv']:
  src=d/filename;dest=out/(filename+'.gz')
  with src.open('rb') as fi,gzip.open(dest,'wb',compresslevel=1) as fo:shutil.copyfileobj(fi,fo,1048576)
 for filename in ['summary.json','input_manifest.json']:shutil.copy2(d/filename,out/('original_'+filename))
 rep=d.parent/'reaggregated/representative_annotations.tsv.gz'
 if rep.exists():shutil.copy2(rep,out/'original_representative_annotations.tsv.gz')
 pd.DataFrame(summaries).fillna(0).to_csv(out/'acquisition_conflict_summary.tsv',sep='\t',index=False)
 write_json(out/'SEARCH_INPUTS.json',source_lists);write_json(out/'INPUT_RECEIPTS.json',receipts)
 summary=dict(id=task['id'],cohort=task['cohort'],arm=task['arm'],population=task['population'],acquisitions=len(names),peptides=len(mapping),groups=len(counts),candidate_proteins_exported=len(candidate_cache),one_assigned_peptide_groups=sum(v==1 for v in counts.values()),one_observed_support_peptide_groups=sum(int(r['observed_support_peptides'])==1 for r in group_rows),species_conflicting_groups=sum(len(g[0])>1 for g in gstats.values()),genus_conflicting_groups=sum(len(g[1])>1 for g in gstats.values()),KO_differing_set_groups=sum(len(g[2])>1 for g in gstats.values()),all_source_hashes_match=True,all_reported_members_sequence_compatible=True,all_dictionary_peptides_observed=True,all_representatives_exported=True,elapsed_seconds=time.time()-started,scope='Exact I/L-equivalent substring compatibility in each acquisition searched target FASTA; not a digestion-filtered or biological-origin claim. Candidate group summaries union alternatives for assigned peptides. No new FDR calibration.')
 write_json(out/'COMPLETE.json',summary);print(json.dumps(summary),flush=True)

if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('task',nargs='?',default='list');args=parser.parse_args();alltasks=tasks()
 HERE.mkdir(exist_ok=True);write_json(HERE/'TASKS.json',{k:{a:b for a,b in v.items() if a not in ['search_inputs','static_inputs']} for k,v in alltasks.items()})
 if args.task=='list':print('\n'.join(alltasks))
 elif args.task=='all':
  for key in sorted(alltasks,key=lambda k:(not k.startswith('CAMPI'),k)):run(alltasks[key])
 else:run(alltasks[args.task])
