"""Reproduce the frozen KO precision estimator, then compare revised annotations."""
from pathlib import Path
from itertools import combinations
from datetime import datetime,timezone
import hashlib,json,os,sys
import numpy as np
import pandas as pd
from scipy import sparse

H=Path(__file__).resolve().parent
R=next(p for p in H.parents if p.name=='FastaLake')
A33=R/'MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19'
OLD=A33/'results/HSTOOL__simplified/source_data/Figure3'
NEW=H/'profiles/HSTOOL_noclust'
TIERS=['TechnicalReplicates','WorkflowReplicates','OnePiece','Intra','Interindividual']
KEYS=['design','i','j'];METRICS=['r_per_pair','cv_per_pair_pct','n_features_per_pair']
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def load(p):return json.loads(Path(p).read_text())
def dump(p,x):Path(p).write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')
def meta(p):return dict(path=str(p),sha256=sha(p),bytes=Path(p).stat().st_size)
def write(frame,p):frame.to_csv(p,sep='\t',index=False)
def measure(v):
 return dict(r_per_pair=float(np.corrcoef(np.log2(v).T)[0,1]) if len(v)>2 else np.nan,
             cv_per_pair_pct=float(np.median(v.std(axis=1,ddof=0)/v.mean(axis=1))*100) if len(v) else np.nan,
             n_features_per_pair=len(v))
def valid(frame,samples):
 assert frame.index.is_unique and frame.columns.is_unique
 assert set(frame.columns)==set(samples)
 assert np.isfinite(frame.to_numpy()).all() and (frame.to_numpy()>=0).all()
 assert 'FunctionallyUnassigned' not in frame.index
 return frame.reindex(columns=samples)
def pairs(frame,pairlist):
 rows=[]
 for design,a,b in pairlist:
  v=frame[[a,b]].to_numpy();v=v[(v>0).all(axis=1)]
  rows.append(dict(design=design,i=a,j=b,**measure(v)))
 return pd.DataFrame(rows)
def summarize(frame,groupcols):
 return frame.groupby(groupcols,sort=False)[METRICS].median().reset_index()

def main(mode):
 assert mode in ['baseline','refit']
 im=load(H/'HSTOOL_PRECISION_INPUTS.json');cm=load(H/'HSTOOL_PRECISION_CODE.json')
 for record in im['files']+cm['files']:assert sha(record['path'])==record['sha256'],record['path']
 selection=load(OLD.parents[1]/'RUN_SELECTION.json')
 assert selection['selection']=='HSTOOL__simplified'
 for item in selection['resolved_inputs']:
  if 'HSTOOL_KO_local_candidate_' in item.get('path',''):assert sha(item['path'])==item['sha256']
 design=pd.read_csv(OLD/'design.tsv',sep='\t',dtype=str).fillna('')
 assert len(design)==58 and design.sample_id.is_unique
 assert set(design.group)==set(TIERS),sorted(set(design.group))
 samples=design.sample_id.tolist();groups={t:design.loc[design.group.eq(t),'sample_id'].tolist() for t in TIERS}
 assert [len(groups[t]) for t in TIERS]==[12,15,10,11,10]
 pairlist=[(t,a,b) for t in TIERS for a,b in combinations(groups[t],2)];assert len(pairlist)==316
 old={c:valid(pd.read_csv(Path(selection['matrices_dir'])/f'HSTOOL_KO_local_candidate_{c}.tsv.gz',sep='\t',index_col=0).drop(index='FunctionallyUnassigned',errors='ignore'),samples) for c in ['all','host','microbial']}
 out=H/('hstool_precision_'+mode);assert not out.exists();out.mkdir()
 if mode=='baseline':
  actual=pd.concat([pairs(frame,pairlist).assign(cls=c) for c,frame in old.items()],ignore_index=True)
  ref=pd.read_csv(OLD/'pyramid_pairwise_local_candidate.tsv',sep='\t');ref=ref[ref.space.eq('KO')]
  keys=['cls']+KEYS;actual=actual.sort_values(keys).reset_index(drop=True);ref=ref.sort_values(keys).reset_index(drop=True)
  assert len(actual)==len(ref)==948 and actual[keys].equals(ref[keys])
  assert np.array_equal(actual.n_features_per_pair,ref.n_features_per_pair)
  maxdiff={}
  for metric in METRICS[:2]:
   a=actual[metric].to_numpy();b=ref[metric].to_numpy()
   assert np.allclose(a,b,rtol=1e-10,atol=1e-10,equal_nan=True),metric
   maxdiff[metric]=float(np.nanmax(np.abs(a-b)))
  write(actual,out/'pairwise.tsv');write(summarize(actual,['cls','design']),out/'summary.tsv')
  checks=dict(reproduced_pair_rows=948,max_absolute_error=maxdiff)
 else:
  gate=load(H/'hstool_precision_baseline/COMPLETE.json')
  assert gate['status']=='PASS' and gate['code_manifest_sha256']==sha(H/'HSTOOL_PRECISION_CODE.json')
  assert gate['input_manifest_sha256']==sha(H/'HSTOOL_PRECISION_INPUTS.json')
  assert load(H/'HSTOOL_BASELINE_REVIEW.json')['status']=='PASS'
  assert load(NEW/'COMPLETE.json')['status']=='PASS'
  new={}
  for c in ['all','other_unresolved']:
   axes=load(NEW/f'KO_{c}_axes.json');m=sparse.load_npz(NEW/f'KO_{c}.npz')
   assert m.shape==(len(axes['samples']),len(axes['features']))
   new[c]=valid(pd.DataFrame(m.toarray().T,index=axes['features'],columns=axes['samples']),samples)
  raw=pd.concat([pairs(frame,pairlist).assign(profile='new_'+c) for c,frame in new.items()],ignore_index=True)
  raw=pd.concat([raw,*[pairs(old[c],pairlist).assign(profile='old_'+c) for c in ['all','microbial']]],ignore_index=True)
  assert len(raw)==1264
  matched=[];attrition=[]
  for oc,nc in [('all','all'),('microbial','other_unresolved')]:
   comparison='old_'+oc+'__new_'+nc
   union=old[oc].index.union(new[nc].index)
   a=old[oc].reindex(union,fill_value=0);b=new[nc].reindex(union,fill_value=0)
   for tier,i,j in pairlist:
    x=a[[i,j]].to_numpy();y=b[[i,j]].to_numpy()
    am=(x>0).all(axis=1);bm=(y>0).all(axis=1);common=am&bm
    nold=int(am.sum());nnew=int(bm.sum());ncommon=int(common.sum())
    attrition.append(dict(comparison=comparison,design=tier,i=i,j=j,old_eligible=nold,new_eligible=nnew,common=ncommon,lost=int((am&~bm).sum()),gained=int((bm&~am).sum())))
    for method,v in [('old',x),('new',y)]:matched.append(dict(comparison=comparison,method=method,design=tier,i=i,j=j,**measure(v[common])))
  matched=pd.DataFrame(matched);attrition=pd.DataFrame(attrition)
  assert len(matched)==1264 and len(attrition)==632
  assert (attrition.old_eligible==attrition.common+attrition.lost).all()
  assert (attrition.new_eligible==attrition.common+attrition.gained).all()
  wide=matched.pivot(index=['comparison']+KEYS,columns='method',values=METRICS)
  assert (wide['n_features_per_pair']['old']==wide['n_features_per_pair']['new']).all()
  write(raw,out/'pairwise_native_features.tsv');write(summarize(raw,['profile','design']),out/'native_feature_summary.tsv')
  write(matched,out/'pairwise_matched_features.tsv');write(summarize(matched,['comparison','method','design']),out/'matched_feature_summary.tsv')
  write(attrition,out/'pairwise_feature_attrition.tsv');write(attrition.groupby(['comparison','design'],sort=False).median(numeric_only=True).reset_index(),out/'feature_attrition_summary.tsv')
  cov=pd.read_csv(NEW/'coverage.tsv',sep='\t');cov=cov[cov.level.eq('KO')].copy();assert len(cov)==116
  cov['allocated_intensity']=cov.single_eligible_term+cov.multiple_eligible_terms
  assert np.allclose(cov.allocated_intensity+cov.no_eligible_term+cov.human_withheld,cov.total_LFQ,rtol=1e-9)
  for field in ['allocated_intensity','no_eligible_term','human_withheld']:cov[field+'_fraction']=cov[field]/cov.total_LFQ
  cov=cov.merge(design[['sample_id','group']],left_on='sample',right_on='sample_id',validate='m:1')
  write(cov,out/'coverage.tsv')
  write(cov.groupby(['profile','group'],sort=False)[[c for c in cov if c.endswith('_fraction')]].median().reset_index(),out/'coverage_summary.tsv')
  checks=dict(pair_count_per_profile=316,new_profiles=2,matched_comparisons=2,attrition_partitions_checked=632,coverage_rows=116,undefined_r_pairs=int(raw.r_per_pair.isna().sum()),protein_groups_refit=False)
 write(design,out/'design.tsv')
 dump(out/'COMPLETE.json',dict(status='PASS',utc=datetime.now(timezone.utc).isoformat(),mode=mode,job_id=os.environ.get('SLURM_JOB_ID'),checks=checks,input_manifest_sha256=sha(H/'HSTOOL_PRECISION_INPUTS.json'),code_manifest_sha256=sha(H/'HSTOOL_PRECISION_CODE.json'),outputs=[meta(p) for p in sorted(out.iterdir())],scope='Descriptive KO precision on the frozen58acquisition roster; no species attribution, inference or manuscript adoption.'))
 print(json.dumps(checks,indent=2))
if __name__=='__main__':main(sys.argv[1])
