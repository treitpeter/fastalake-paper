"""Frozen SIHUMIx-only yield decomposition; no searches or q-value refitting."""
from pathlib import Path
from datetime import datetime, timezone
from collections import Counter
import hashlib, json, os, sys, re
import numpy as np
import pandas as pd
from scipy.special import gammaln
from lxml import etree as E

HERE=Path(__file__).resolve().parent
FL=HERE.parents[4]
J=HERE.parents[1]
ARMS=['FL','REF','MG','FLR','FLRM']
ALIAS={'FL':'FL','REF':'REF','MG':'MG','FL_R':'FLR','FL_R_M':'FLRM'}
TASKS=dict(FL='CAMPI_FL_UNIPROT',FLR='CAMPI_FL_REF',FLRM='CAMPI_FL_REF_MG',REF='CAMPI_REF',MG='CAMPI_MG')
SEED=2026093017
def load(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def meta(p,hashit=True):
 p=Path(p).resolve();s=p.stat();r=dict(path=str(p),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
 if hashit:r['sha256']=sha(p)
 return r
def verify(r):
 p=Path(r['path']);s=p.stat();assert (s.st_size,s.st_mtime_ns)==(r['bytes'],r['mtime_ns']),p
 if 'sha256' in r:assert sha(p)==r['sha256'],p
def table(p,rows):pd.DataFrame(rows).to_csv(p,sep='\t',index=False)
def canonical(s):return re.sub(r'\[[^\]]*\]','',s).replace('I','L')
def guard():
 m=load(HERE/'INPUTS.json')
 for r in m['small_inputs']+m['code']:verify(r)
 assert len(m['tasks'])==15 and all(t['sample'].startswith('S') for t in m['tasks'])
 return m

def expected_unique(freq,N,k):
 freq=np.asarray(freq,dtype=np.int64);assert 0<=k<=N and np.all((freq>0)&(freq<=N))
 if k==N:return float(len(freq))
 if k==0:return 0.
 # P(no occurrence) = C(N-f,k)/C(N,k); zero if insufficient non-hits.
 possible=N-freq>=k;no=np.zeros(len(freq),float);f=freq[possible]
 no[possible]=np.exp(gammaln(N-f+1)-gammaln(N-f-k+1)-gammaln(N+1)+gammaln(N-k+1))
 return float(np.sum(1-np.clip(no,0,1)))

class HashReader:
 def __init__(self,p):self.f=open(p,'rb');self.h=hashlib.sha256()
 def read(self,n=-1):
  b=self.f.read(n);self.h.update(b);return b
 def close(self):self.f.close()

def read_ms2(record):
 verify(record);reader=HashReader(record['path']);ids=set();times=[];count=0
 for _,e in E.iterparse(reader,events=('end',),tag='{*}spectrum',huge_tree=True):
  cv=e.find("{*}cvParam[@accession='MS:1000511']")
  if cv is not None and cv.get('value')=='2':
   ident=e.get('id');assert ident not in ids;ids.add(ident);count+=1
   rt=e.find(".//{*}cvParam[@accession='MS:1000016']")
   if rt is not None:times.append(float(rt.get('value'))/(60 if rt.get('unitAccession')=='UO:0000010' else 1))
  e.clear()
  while e.getprevious() is not None:del e.getparent()[0]
 reader.close();verify(record);assert count==len(ids)>0
 return ids,dict(ms2=count,rt_min=min(times),rt_max=max(times),mzml_sha256=reader.h.hexdigest())

def normalized_config(c):
 d={k:v for k,v in c.items() if k not in ['mzml_paths','output_paths','output_directory','quant']}
 d['database']=dict(d['database']);d['database'].pop('fasta',None);return d

def audit(t,m):
 sample=t['sample'];out=HERE/'results'/sample;out.mkdir(parents=True,exist_ok=False)
 ids,inventory=read_ms2(t['raw']);counts=[];sourcehashes=[];frames={};params={}
 figure=pd.read_csv(HERE/'inputs/figure_metrics.tsv',sep='\t').set_index(['sample','strategy'])
 for arm in ARMS:
  inp=t['arms'][arm];verify(inp['search']);verify(inp['config']);cfg=load(inp['config']['path'])
  assert len(cfg['mzml_paths'])==1 and Path(cfg['mzml_paths'][0]).samefile(t['raw']['path'])
  params[arm]=normalized_config(cfg)
  cols=['scannr','rank','label','peptide','charge','peptide_q','hyperscore','spectrum_q']
  a=pd.read_csv(inp['search']['path'],sep='\t',usecols=cols);assert a['rank'].eq(1).all()
  assert a.scannr.is_unique and set(a.scannr)<=ids
  a['canonical']=a.peptide.map(canonical);a['length']=a.canonical.str.len()
  a['accepted']=a.label.eq(1)&a.peptide_q.le(.01)&a.length.between(5,50)
  a=a.set_index('scannr');frames[arm]=a
  for lo in [5,9]:
   accepted=a[a.accepted&a.length.ge(lo)];freq=accepted.canonical.value_counts()
   assert expected_unique(freq.to_numpy(),len(ids),len(ids))==len(freq)
   freq.rename_axis('peptide').rename('accepted_psms').to_csv(out/f'{arm}_{lo}_frequencies.tsv',sep='\t')
   row=dict(sample=sample,arm=arm,min_length=lo,ms2=len(ids),reported_candidates=len(a),accepted_psms=len(accepted),peptides=len(freq),identified_fraction=len(accepted)/len(ids),distinct_per_accepted=len(freq)/len(accepted))
   counts.append(row)
   if lo==5:
    assert len(freq)==int(figure.loc[(sample,arm),'peptides'])
    assert len(accepted)==int(figure.loc[(sample,arm),'psms'])
  sourcehashes.extend([meta(inp['search']['path']),meta(inp['config']['path'])]);verify(inp['search'])
 # Search settings must agree to attribute paired differences to existing database arms.
 dump(out/'SEARCH_PARAMETERS.json',params)
 assert len({json.dumps(c,sort_keys=True) for c in params.values()})==1,'Search parameter mismatch'
 fates=[]
 for a,b in [('FL','REF'),('FL','MG'),('FL','FLR'),('FL','FLRM'),('REF','FLR'),('MG','FLRM')]:
  x=frames[a][frames[a].accepted];y=frames[b].reindex(x.index)
  same=y.peptide.eq(x.peptide)&y.charge.eq(x.charge)&y.label.eq(1)
  flags={'same_form_pass':same&y.accepted.eq(True),'same_form_threshold_loss':same&~y.accepted.eq(True),
         'changed_target_winner':y.label.eq(1)&~same,'decoy_winner':y.label.eq(-1),'no_reported_candidate':y.label.isna()}
  assert sum(v.astype(int) for v in flags.values()).eq(1).all()
  for name,v in flags.items():fates.append(dict(sample=sample,source_arm=a,destination_arm=b,outcome=name,spectra=int(v.sum()),denominator=len(x),same_form_raw_unchanged=int((v&same&y.hyperscore.eq(x.hyperscore)).sum())))
 table(out/'counts.tsv',counts);table(out/'same_spectrum_fates.tsv',fates)
 gap=pd.read_csv(t['gap']['path'],sep='\t');verify(t['gap'])
 ref=set(frames['REF'].loc[frames['REF'].accepted,'canonical']);assert set(gap.peptide)==ref and gap.peptide.is_unique
 for arm,col in [('FL','accepted_FL'),('FLR','accepted_FL_R'),('FLRM','accepted_FL_R_M')]:
  actual=set(frames[arm].loc[frames[arm].accepted,'canonical']);assert gap[col].astype(bool).equals(gap.peptide.isin(actual))
 rows=[]
 for stage,g in gap.groupby('stage'):
  rows.append(dict(sample=sample,stage=stage,reference_peptides=len(ref),peptides=len(g),recovered_add_ref=int(g.accepted_FL_R.sum()),recovered_add_ref_mg=int(g.accepted_FL_R_M.sum())))
 table(out/'availability_partition.tsv',rows)
 dump(out/'COMPLETE.json',dict(status='PASS',sample=sample,inventory=inventory,inputs=sourcehashes+[meta(t['gap']['path'])],outputs=[meta(p) for p in sorted(out.iterdir()) if p.is_file()],job_id=os.getenv('SLURM_JOB_ID')))
 print('PASS',sample,inventory['ms2'],flush=True)

def group_audit(m):
 figure=pd.read_csv(HERE/'inputs/figure_metrics.tsv',sep='\t').set_index(['sample','strategy']);unions=[]
 for arm,r in m['matrices'].items():
  verify(r);a=pd.read_csv(r['path'],sep='\t',index_col=0).fillna(0)
  assert a.index.is_unique and a.columns.is_unique and set(a.columns)=={t['sample'] for t in m['tasks']}
  assert np.isfinite(a.to_numpy()).all() and (a.to_numpy()>=0).all()
  for sample in a:
   assert int((a[sample]>0).sum())==int(figure.loc[(sample,arm),'quantified_groups'])
  cols=[c for c in a if c.startswith('S11_')];assert len(cols)==4
  unions.append(dict(sample='S11_union',strategy=arm,quantified_groups=int(a[cols].gt(0).any(axis=1).sum()),note='Union of fixed study-group rows within this arm; no new inference'))
 return unions

def prepare():
 inputs=HERE/'inputs';inputs.mkdir(exist_ok=False);(HERE/'results').mkdir()
 sources={'figure_metrics.tsv':J/'audits/CAMPI_workflow_variability.tsv','workflow_metadata.tsv':J/'CAMPI_acquisitions.tsv',
 'method_overview.xlsx':FL/'00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/methods_comparison/campi_method_overview.xlsx'}
 import shutil
 for n,p in sources.items():shutil.copy2(p,inputs/n)
 old=load(FL/'MM_review_analysis/A10_4_3_4_5_identity_recovery/inputs/manifest.json');tasks=[]
 for t in old['tasks']:
  if t['cohort']!='CAMPI':continue
  arms={ALIAS[a]:{k:meta(v['path'],k=='config') for k,v in x.items()} for a,x in t['arms'].items()}
  tasks.append(dict(sample=t['sample'],raw=meta(t['raw'],False),arms=arms,gap=meta(FL/'MM_review_analysis/A14_4_3_campi_gap/results'/t['sample']/'reference_gap.tsv')))
 matrices={a:meta(FL/'simplified_grouping_results/manuscript_cohorts_2026-09-10_2200/cohorts'/task/'reaggregated/simplified_study_group_matrix.tsv.gz') for a,task in TASKS.items()}
 m=dict(scope='SIHUMIx mock only',tasks=tasks,matrices=matrices,small_inputs=[meta(p) for p in inputs.iterdir()],source_provenance=[meta(p) for p in sources.values()],code=[meta(HERE/'analyze.py'),meta(HERE/'PROTOCOL.md'),meta(HERE/'run.sh')],pilot_samples=['S01','S11_Fraction4'])
 assert len(tasks)==15;dump(HERE/'INPUTS.json',m);group_audit(m)
 dump(HERE/'PREFLIGHT.json',dict(status='PASS',acquisitions=15,arms=5,quantified_group_counts_reproduced=75,scope=m['scope']))

def summarize():
 m=guard();groups=group_audit(m);counts=[];fates=[];gaps=[];inventories=[];frequencies={}
 for t in m['tasks']:
  p=HERE/'results'/t['sample'];r=load(p/'COMPLETE.json');assert r['status']=='PASS'
  for x in r['outputs']:verify(x)
  inventories.append(dict(sample=t['sample'],**r['inventory']))
  counts.append(pd.read_csv(p/'counts.tsv',sep='\t'));fates.append(pd.read_csv(p/'same_spectrum_fates.tsv',sep='\t'));gaps.append(pd.read_csv(p/'availability_partition.tsv',sep='\t'))
  for arm in ARMS:
   for lo in [5,9]:frequencies[t['sample'],arm,lo]=pd.read_csv(p/f'{arm}_{lo}_frequencies.tsv',sep='\t').set_index('peptide').accepted_psms
 out=HERE/'summary';out.mkdir(exist_ok=False)
 c=pd.concat(counts,ignore_index=True);assert len(c)==150
 c.to_csv(out/'primary_counts.tsv',sep='\t',index=False);pd.concat(fates).to_csv(out/'same_spectrum_fates.tsv',sep='\t',index=False);pd.concat(gaps).to_csv(out/'availability_partition.tsv',sep='\t',index=False);table(out/'spectral_inventory.tsv',inventories)
 k=int(c.ms2.min());kp=int(c.accepted_psms.min());rare=[]
 for r in c.to_dict('records'):
  f=frequencies[r['sample'],r['arm'],r['min_length']].to_numpy()
  rare.append(dict(**r,equal_ms2_draw=k,expected_peptides_equal_ms2=expected_unique(f,r['ms2'],k),equal_accepted_psm_draw=kp,expected_peptides_equal_accepted_psms=expected_unique(f,r['accepted_psms'],kp)))
 table(out/'depth_standardization.tsv',rare)
 figure=pd.read_csv(HERE/'inputs/figure_metrics.tsv',sep='\t');union=[]
 for arm in ARMS:
  sets=[set(frequencies[t['sample'],arm,5].index) for t in m['tasks'] if t['sample'].startswith('S11_')]
  union.append(dict(sample='S11_union',strategy=arm,peptides=len(set.union(*sets)),quantified_groups=next(x['quantified_groups'] for x in groups if x['strategy']==arm)))
 table(out/'fraction_workflow_union.tsv',union)
 clean=figure[~figure['sample'].str.startswith('S11_')];workflow=pd.concat([clean,pd.DataFrame(union)],ignore_index=True)
 spreads=[];decomp=[];contrasts=[]
 for scope,a in [('all15',figure),('unfractionated11',clean),('complete_workflows12',workflow)]:
  for metric in ['peptides','quantified_groups']:
   for arm,g in a.groupby('strategy'):
    y=g[metric].to_numpy();spreads.append(dict(scope=scope,metric=metric,arm=arm,n=len(y),minimum=float(y.min()),maximum=float(y.max()),max_min=float(y.max()/y.min()),cv=float(y.std(ddof=1)/y.mean()),sd_log=float(np.std(np.log(y),ddof=1))))
   pivot=a.pivot(index='sample',columns='strategy',values=metric).reindex(columns=ARMS);y=np.log(pivot.to_numpy());grand=y.mean();ac=y.mean(axis=1,keepdims=True);ar=y.mean(axis=0,keepdims=True)
   total=float(((y-grand)**2).sum());terms={'acquisition':float((len(ARMS)*(ac-grand)**2).sum()),'database_arm':float((len(y)*(ar-grand)**2).sum()),'interaction_residual':float(((y-ac-ar+grand)**2).sum())};assert np.isclose(sum(terms.values()),total)
   for term,ss in terms.items():decomp.append(dict(scope=scope,metric=metric,term=term,ss=ss,total_ss=total,fraction_ss=ss/total,interpretation='Descriptive balanced log-yield decomposition; not causal'))
   for base,dest in [('FL','REF'),('FL','MG'),('FL','FLR'),('FL','FLRM'),('REF','FLR'),('MG','FLRM')]:
    for sample in pivot.index:contrasts.append(dict(scope=scope,metric=metric,sample=sample,base=base,destination=dest,ratio=float(pivot.loc[sample,dest]/pivot.loc[sample,base]),difference=float(pivot.loc[sample,dest]-pivot.loc[sample,base])))
 table(out/'yield_spreads.tsv',spreads);table(out/'descriptive_variance_decomposition.tsv',decomp);table(out/'paired_database_effects.tsv',contrasts)
 within=[]
 for arm,a in figure.groupby('strategy'):
  a=a.set_index('sample')
  for metric in ['peptides','quantified_groups']:
   for label,low,high in [('Kleiner_140_to_260_load_covaries',['S01','S02'],['S03','S04']),('Kleiner_140_to_460_load_covaries',['S01','S02'],['S05','S06']),('CEA_60_to_120_fixed_reported_load',['S09'],['S10'])]:
    l=float(a.loc[low,metric].mean());h=float(a.loc[high,metric].mean());within.append(dict(arm=arm,metric=metric,contrast=label,low_mean=l,high_mean=h,ratio=h/l,n_low=len(low),n_high=len(high)))
 table(out/'within_lab_contrasts.tsv',within)
 dump(out/'COMPLETE.json',dict(status='COMPUTATION_COMPLETE_REQUIRES_REVIEW',scope='15 SIHUMIx acquisitions only; no complexity inference',primary_count_checks=75,group_count_checks=75,spectrum_inventory_checks=15,seed=SEED,postsearch_rarefaction=True,scientific_review='pending',inputs=[meta(HERE/'INPUTS.json')],outputs=[meta(p) for p in out.iterdir() if p.is_file()],job_id=os.getenv('SLURM_JOB_ID')))
 print('SUMMARY COMPLETE',flush=True)

if __name__=='__main__':
 action=sys.argv[1]
 if action=='prepare':prepare()
 elif action=='collect':summarize()
 else:
  m=guard();group_audit(m)
  if action=='pilot':chosen=[t for t in m['tasks'] if t['sample'] in m['pilot_samples']]
  elif action=='full':
   assert load(HERE/'PILOT_GATE.json')['status']=='PASS';chosen=[t for t in m['tasks'] if t['sample'] not in m['pilot_samples']]
  else:raise ValueError(action)
  for t in chosen:audit(t,m)
  if action=='pilot':dump(HERE/'PILOT_GATE.json',dict(status='PASS',samples=m['pilot_samples'],receipts=[meta(HERE/'results'/t['sample']/'COMPLETE.json') for t in chosen],job_id=os.getenv('SLURM_JOB_ID')))
