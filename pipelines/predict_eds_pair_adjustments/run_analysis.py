"""Fixed-patient diagnostic decomposition; never select a preferred model by significance."""
from pathlib import Path
import os, sys, json, hashlib, itertools
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr, t
from statsmodels.api import OLS
H=Path(__file__).resolve().parent; A=H.parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(name,x): (H/name).write_text(json.dumps(x,indent=2)+'\n')
frozen=json.loads((H/'FREEZE.json').read_text())
for v in frozen['files']: assert sha(v['path'])==v['sha256'],v['path']
inputs=json.loads((A/'INPUTS.json').read_text())
for v in inputs.values(): assert sha(v['path'])==v['sha256'],v['path']
def read(k,**kw):return pd.read_csv(inputs[k]['path'],sep='\t',**kw)
order=read('order'); z=np.load(inputs['heatmap']['path']); order=order.set_index('sample_id').loc[z['samples']].reset_index()
assert len(order)==order.patient_id.nunique()==85
meta=read('metabolites',index_col=0).loc[order.microbiome_id]; meta.index=order.sample_id
mat={k:read(k,index_col=0).drop(index='FunctionallyUnassigned',errors='ignore').loc[:,order.sample_id].T for k in ['KO','Pathway']}
assert order[['Human','pg_parsimony','antibiotic_state','clinical_group']].notna().all().all()
assert np.isfinite(order[['Human','pg_parsimony']]).all().all() and (order.pg_parsimony>0).all()
cats={k:pd.get_dummies(order[k],drop_first=True).to_numpy(float) for k in ['antibiotic_state','clinical_group']}
names=['host','depth','antibiotics','clinical']
subsets=[tuple(c) for n in range(5) for c in itertools.combinations(names,n)]
models=['none' if not c else '+'.join(c) for c in subsets]
pairs=json.loads((H/'TARGETS.json').read_text())
prior=pd.read_csv(A/'augmented_pair_screen.tsv',sep='\t')

def blocks(idx):
 return {'host':rankdata(order.Human.to_numpy()[idx])[:,None], 'depth':rankdata(np.log2(order.pg_parsimony.to_numpy()[idx]))[:,None], 'antibiotics':cats['antibiotic_state'][idx], 'clinical':cats['clinical_group'][idx]}
def estimate(x,y,b,subset):
 Z=np.column_stack([np.ones(len(x))]+[b[k] for k in subset]); ranks=np.column_stack([rankdata(x),rankdata(y)])
 res=ranks-Z @ np.linalg.lstsq(Z,ranks,rcond=None)[0]
 r=float(np.corrcoef(res.T)[0,1]); rank=int(np.linalg.matrix_rank(Z)); df=len(x)-rank-1
 p=float(2*t.sf(abs(r)*np.sqrt(df/max(1e-16,1-r*r)),df))
 return r,p,rank,res,Z
masks={}
for pair in pairs:
 x=mat[pair['layer']][pair['protein']].to_numpy(float); y=meta[pair['metabolite']].to_numpy(float)
 masks[pair['id']]=np.isfinite(x)&np.isfinite(y)&(x>0)&(y>0)
common=np.logical_and.reduce(list(masks.values())); assert common.sum()>=43
rosters=[];checks=[]; results=[]; bootstrap_n=1000 if sys.argv[1]=='full' else 0
if sys.argv[1]=='full': assert json.loads((H/'PILOT.json').read_text())['script_sha256']==sha(__file__)
for pi,pair in enumerate(pairs):
 for scope in ['pair_complete','all_pairs_complete']:
  idx=np.where(masks[pair['id']] if scope=='pair_complete' else common)[0]; n=len(idx)
  for j in idx:rosters.append(dict(pair=pair['id'],scope=scope,sample_id=order.sample_id.iloc[j],patient_id=order.patient_id.iloc[j]))
  b=blocks(idx)
  for norm in ['native','annotated_relative']:
   xx=mat[pair['layer']][pair['protein']].to_numpy(float)
   if norm!='native':xx=xx/mat[pair['layer']].sum(axis=1).to_numpy(float)
   yy=meta[pair['metabolite']].to_numpy(float); x=xx[idx];y=yy[idx]
   estimates=[estimate(x,y,b,s) for s in subsets]; base=estimates[0][0]
   assert np.isclose(base,spearmanr(x,y).statistic,atol=1e-12)
   for r,p,rank,res,Z in estimates:
    # Independent projection and regression-t check.
    U,S,_=np.linalg.svd(Z,full_matrices=False);U=U[:,:rank];ranks=np.column_stack([rankdata(x),rankdata(y)]);ind=ranks-U@(U.T@ranks)
    assert np.isclose(r,np.corrcoef(ind.T)[0,1],atol=1e-10)
    fit=OLS(rankdata(y),np.column_stack([Z,rankdata(x)])).fit()
    assert np.isclose(p,fit.pvalues[-1],atol=1e-9)
   if scope=='pair_complete':
    old=prior[(prior.layer==pair['layer'])&(prior.metabolite==pair['metabolite'])&(prior.protein_feature==pair['protein'])&(prior.normalization==norm)]
    assert len(old)==3
    assert np.isclose(base,old[old.model=='unadjusted'].rho.iloc[0],atol=1e-12)
    # Exact original replay uses covariate ranks calculated on all 85 patients.
    oldblocks=blocks(np.arange(85));oldblocks={k:v[idx] for k,v in oldblocks.items()}
    replay=estimate(x,y,oldblocks,tuple(names))[0]
    assert np.isclose(replay,old[old.model=='host_depth_antibiotic_clinical'].rho.iloc[0],atol=1e-12)
    checks.append(dict(pair=pair['id'],normalization=norm,n=n,original_full_adjusted_rho=replay,subset_rank_full_adjusted_rho=estimates[-1][0]))
   draws=[]
   if bootstrap_n and scope=='pair_complete':
    rng=np.random.default_rng(2026093032+pi)
    for _ in range(bootstrap_n):
     jj=rng.integers(0,n,n);bb=blocks(idx[jj]);draws.append([estimate(x[jj],y[jj],bb,s)[0] for s in subsets])
    draws=np.array(draws);assert np.isfinite(draws).all()
   for j,(s,(r,p,rank,_,_)) in enumerate(zip(subsets,estimates)):
    row=dict(**pair,scope=scope,normalization=norm,model=models[j],n=n,rho=r,delta_from_unadjusted=r-base,covariate_df=rank-1,nominal_p=p)
    if len(draws):
     lo,hi=np.quantile(draws[:,j],[.025,.975]);dl,dh=np.quantile(draws[:,j]-draws[:,0],[.025,.975]);row.update(lower95=lo,upper95=hi,delta_lower95=dl,delta_upper95=dh,bootstrap_draws=bootstrap_n)
    results.append(row)
  print(pair['id'],scope,n,'complete',flush=True)
res=pd.DataFrame(results)
assert len(res)==len(pairs)*2*2*16
assert res.groupby(['id','scope']).n.nunique().eq(1).all()
pd.DataFrame(checks).to_csv(H/('pilot_replay.tsv' if not bootstrap_n else 'original_replay.tsv'),sep='\t',index=False)
if not bootstrap_n:
 save('PILOT.json',dict(status='PASS',script_sha256=sha(__file__),patients=85,common_patients=int(common.sum()),models_checked=len(res),exact_original_replay_checks=len(checks),independent_projection_and_p_checks=True,job_id=os.environ.get('SLURM_JOB_ID')));sys.exit()
res.to_csv(H/'all_models.tsv',sep='\t',index=False)
pd.DataFrame(rosters).to_csv(H/'patient_rosters.tsv',sep='\t',index=False)
order.to_csv(H/'patient_covariates.tsv',sep='\t',index=False)
order[['Human','pg_parsimony']].corr(method='spearman').to_csv(H/'host_depth_correlation.tsv',sep='\t')
counts=[]
for pid,mask in masks.items():
 for var in cats:
  for category,num in order.loc[mask,var].value_counts(dropna=False).items():counts.append(dict(pair=pid,covariate=var,category=category,patients=int(num)))
pd.DataFrame(counts).to_csv(H/'category_counts.tsv',sep='\t',index=False)
save('COMPLETE.json',dict(status='COMPUTATION_COMPLETE_PENDING_REVIEW',job_id=os.environ.get('SLURM_JOB_ID'),script_sha256=sha(__file__),pairs=len(pairs),models=len(res),bootstrap_draws=bootstrap_n,all_pairs_common_patients=int(common.sum()),original_replays=len(checks),scope='Post hoc descriptive sensitivity; intervals pointwise, nominal p values not discovery evidence; all models retained, no preferred model selected. Original global covariate ranks replayed separately; new models rank within the fixed analysis subset.'))
