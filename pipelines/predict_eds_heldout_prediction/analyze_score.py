"""Fixed-roster, patient-separated EDS candidates; no raw searches or score calibration."""
from pathlib import Path
import argparse, hashlib, json, time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
from scipy.spatial.distance import cdist
from sklearn.linear_model import Ridge
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, mean_absolute_error
from statsmodels.stats.multitest import multipletests
import statsmodels.api as sm

H=Path(__file__).resolve().parent; ROOT=Path(__import__('os').environ['FASTALAKE_DATA_DIR']); D=H/'source_data'
A19=ROOT/'MM_review_analysis/A19_4_6_biology_robustness/results'
ALPHAS=np.array([.01,.1,1,10,100]); SEEDS=[20260915,20260916,20260917]
PROV={}; FIT_AUDIT=[]
def read(p,**kw):
 p=Path(p); PROV[str(p)]={'path':str(p),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}
 return pd.read_csv(p,sep='\t',**kw)
def write(x,name): x.to_csv(D/(name+'.tsv'),sep='\t',index=False)
def inputs():
 base=A19/'PREDICT__simplified/source_data/Figure4'
 tracks=read(base/'EDS_context_patient_roster.tsv',index_col=0)
 order=read(base/'heatmap_85_patient_order.tsv').set_index('microbiome_id')
 assert len(order)==len(tracks)==tracks.idpatient.nunique()==85
 tracks=tracks.loc[order.index].copy(); tracks['specimen']=tracks.index; tracks['sample_id']=order.sample_id
 historical=read(ROOT/'analyses/enterotypes_2026-09-07/specimen_predictions_and_metadata.tsv',index_col=0)
 assert np.allclose(tracks.ET_dysbiosis,historical.loc[tracks.index,'ET_dysbiosis'])
 frames={}
 for grouping in ['simplified','original']:
  b=A19/('PREDICT__'+grouping)/'source_data'
  raw=read(b/'PREDICT_noclust/microbial_OG.tsv',index_col=0).loc[tracks.sample_id]
  ac=read(b/'Figure4/acquisition_accounting.tsv').set_index('sample_id').loc[tracks.sample_id]
  assert np.array_equal(ac.patient_id.astype(int),tracks.idpatient.astype(int))
  assert (raw.to_numpy()>=0).all() and raw.sum(axis=1).gt(0).all()
  assert np.allclose(ac[['Human','Microbial','Mixed','Unknown']].sum(axis=1),1)
  for field in ['Human','Microbial','Mixed','Unknown','pg_parsimony']:
   tracks[grouping+'_'+field]=ac[field].to_numpy()
  tracks['antibiotic_broad']=ac.antibiotic_state.map({'On':1,'Off':0}).to_numpy()
  frames[grouping]=(raw.to_numpy(float),np.column_stack([np.log10(ac.pg_parsimony),ac.Human]))
 return tracks,frames,historical

def features(raw,z,tr,te,model):
 keep=(raw[tr]>0).mean(axis=0)>=.1
 a=raw[np.ix_(tr,keep)]; b=raw[np.ix_(te,keep)]
 assert keep.any() and (a.sum(1)>0).all() and (b.sum(1)>0).all()
 a=np.sqrt(a/a.sum(1)[:,None]); b=np.sqrt(b/b.sum(1)[:,None])
 zs=z[tr].std(0); zs=np.where(zs>0,zs,1)
 za=(z[tr]-z[tr].mean(0))/zs; zb=(z[te]-z[tr].mean(0))/zs
 if model=='baseline':a,b=za,zb
 elif model=='combined':a,b=np.c_[a,za],np.c_[b,zb]
 return a,b,int(keep.sum())

def hats(a,b):
 ac=a-a.mean(0); bc=b-a.mean(0)
 vals,vec=np.linalg.eigh(ac@ac.T); vals=np.maximum(vals,0)
 cross=bc@ac.T@vec; out=[]; n=len(a)
 for alpha in ALPHAS:
  w=(cross/(vals+alpha))@vec.T
  w=w-w.mean(1)[:,None]+1/n
  out.append(w)
 return np.array(out)

def make_cache(raw,z,seed,model):
 folds=list(KFold(5,shuffle=True,random_state=seed).split(raw)); cache=[]
 for fold,(tr,te) in enumerate(folds):
  assert not set(tr)&set(te)
  inner=[]
  for it,iv in KFold(3).split(tr):
   itr,ival=tr[it],tr[iv]; a,b,n=features(raw,z,itr,ival,model)
   inner.append((itr,ival,hats(a,b)))
   FIT_AUDIT.append(dict(seed=seed,model=model,outer_fold=fold,stage='inner',train=';'.join(map(str,itr)),test=';'.join(map(str,ival)),features=n))
  a,b,n=features(raw,z,tr,te,model); cache.append((tr,te,inner,hats(a,b)))
  FIT_AUDIT.append(dict(seed=seed,model=model,outer_fold=fold,stage='outer',train=';'.join(map(str,tr)),test=';'.join(map(str,te)),features=n))
 return cache

def predict(cache,Y):
 Y=np.asarray(Y); Y=Y[:,None] if Y.ndim==1 else Y
 out=np.full(Y.shape,np.nan); choices=[]
 for tr,te,inner,mat in cache:
  loss=[]
  for itr,ival,ws in inner:
   pp=ws@Y[itr]; loss.append(((pp-Y[ival][None,:,:])**2).mean(1))
  best=np.stack(loss).mean(0).argmin(0)
  pp=mat@Y[tr]; out[te]=pp[best,:,np.arange(Y.shape[1])].T
  choices.append(best)
 assert np.isfinite(out).all()
 return out,np.array(choices)

def distance_predict(raw,z,seed,Y):
 out=np.full(Y.shape,np.nan); audit=[]
 for f,(tr,te) in enumerate(KFold(5,shuffle=True,random_state=seed).split(raw)):
  a,b,n=features(raw,z,tr,te,'protein'); dist=cdist(b,a)/np.sqrt(2)
  for j in range(Y.shape[1]):
   threshold=np.quantile(Y[tr,j],.25); refs=np.flatnonzero(Y[tr,j]<=threshold)
   assert len(refs)>=5
   out[te,j]=np.sort(dist[:,refs],axis=1)[:,:5].mean(1)
   if j==0:audit.append(dict(seed=seed,fold=f,training_EDS_cutoff=threshold,reference_rows=';'.join(map(str,tr[refs])),n_reference=len(refs)))
 return out,audit

def direct_predict(raw,z,seed,model,Y):
 Y=np.asarray(Y); Y=Y[:,None] if Y.ndim==1 else Y
 out=np.full(Y.shape,np.nan); choices=[]
 for tr,te in KFold(5,shuffle=True,random_state=seed).split(raw):
  loss=[]
  for it,iv in KFold(3).split(tr):
   a,b,_=features(raw,z,tr[it],tr[iv],model)
   pp=np.array([Ridge(alpha=alpha,solver='cholesky').fit(a,Y[tr[it]]).predict(b) for alpha in ALPHAS])
   loss.append(((pp-Y[tr[iv]][None,:,:])**2).mean(1))
  best=np.stack(loss).mean(0).argmin(0); choices.append(best)
  a,b,_=features(raw,z,tr,te,model)
  pp=np.array([Ridge(alpha=alpha,solver='cholesky').fit(a,Y[tr]).predict(b) for alpha in ALPHAS])
  out[te]=pp[best,:,np.arange(Y.shape[1])].T
 return out,np.array(choices)

def main(pilot=False):
 start=time.monotonic(); D.mkdir(exist_ok=True)
 meta,frames,hist=inputs(); y=meta.ET_dysbiosis.to_numpy(); rng=np.random.default_rng(SEEDS[0])
 perms=3 if pilot else 199; Y=np.column_stack([y]+[rng.permutation(y) for _ in range(perms)])
 rows=[]; scores=[]; selections=[]; nulls=[]; references=[]; checks=[]; primary={}
 for grouping,(raw,z) in frames.items():
  for seed in SEEDS[:1 if pilot else 3]:
   yy=Y if seed==SEEDS[0] else y[:,None]
   for model in ['protein','baseline','combined','distance']:
    if model=='distance':pred,refs=distance_predict(raw,z,seed,yy);references.extend([dict(grouping=grouping,**r) for r in refs])
    else:
     cache=make_cache(raw,z,seed,model); pred,sel=predict(cache,yy)
     if seed==SEEDS[0]:
      direct,ds=direct_predict(raw,z,seed,model,yy[:,:4])
      err=float(np.max(abs(direct-pred[:,:4])));assert err<1e-7 and np.array_equal(ds,sel[:,:4]),(grouping,model,err)
      checks.append(dict(grouping=grouping,model=model,max_direct_sklearn_error=err,alpha_choices_equal=True))
     for f,best in enumerate(sel[:,0]):selections.append(dict(grouping=grouping,seed=seed,model=model,fold=f,alpha=ALPHAS[best]))
    val=pred[:,0]; rho=spearmanr(y,val).statistic
    row=dict(grouping=grouping,seed=seed,model=model,n=85,rho=rho,MAE=mean_absolute_error(y,val) if model!='distance' else np.nan,r2=r2_score(y,val) if model!='distance' else np.nan,minimum=val.min(),maximum=val.max(),outside_unit_interval=int(((val<0)|(val>1)).sum()))
    if yy.shape[1]>1:
     nr=np.array([spearmanr(yy[:,i],pred[:,i]).statistic for i in range(1,len(yy.T))]);row['permutation_p']=(1+(abs(nr)>=abs(rho)-1e-12).sum())/(1+perms)
     for j,v in enumerate(nr):nulls.append(dict(grouping=grouping,model=model,permutation=j+1,rho=v))
     primary[grouping,model]=val
    rows.append(row)
    scores.append(meta[['specimen','idpatient','ET_dysbiosis','enterotype']].assign(grouping=grouping,seed=seed,model=model,score=val))
   print(grouping,seed,'complete',round(time.monotonic()-start,1),'seconds',flush=True)
 if pilot:
  report=dict(status='PASS',seconds=time.monotonic()-start,patients=85,permutations=perms,direct_checks=checks,scope='pilot only; not final performance')
  (D/'PILOT_VALIDATION.json').write_text(json.dumps(report,indent=2)); print(json.dumps(report,indent=2)); return
 metrics=pd.DataFrame(rows); ok=metrics.permutation_p.notna();metrics.loc[ok,'permutation_q_8']=multipletests(metrics.loc[ok,'permutation_p'],method='fdr_bh')[1]
 boots=np.random.default_rng(20260918).integers(0,85,(2000,85)); improvements=[]; ranges=[]
 for grouping in frames:
  base=primary[grouping,'baseline']
  for model in ['protein','baseline','combined','distance']:
   val=primary[grouping,model]; rs=np.array([spearmanr(y[ix],val[ix]).statistic for ix in boots]);lo,hi=np.quantile(rs,[.025,.975])
   ranges.append(dict(grouping=grouping,model=model,rho_ci_low=lo,rho_ci_high=hi))
   if model in ['protein','combined']:
    delta=abs(base-y)-abs(val-y);ci=np.quantile(delta[boots].mean(1),[.025,.975])
    improvements.append(dict(grouping=grouping,model=model,MAE_baseline_minus_model=delta.mean(),ci_low=ci[0],ci_high=ci[1],bootstrap='2000 patients; fixed held-out predictions'))
 metrics=metrics.merge(pd.DataFrame(ranges),on=['grouping','model']); metrics.loc[metrics.seed.ne(SEEDS[0]),['rho_ci_low','rho_ci_high']]=np.nan
 grouping_check=[]
 for model in ['protein','baseline','combined','distance']:
  a,b=primary['original',model],primary['simplified',model]
  grouping_check.append(dict(model=model,rho_predictions=spearmanr(a,b).statistic,mean_absolute_difference=np.mean(abs(a-b))))
 context=[]
 for (grouping,model),val in primary.items():
  for field,vector in [('host_fraction',meta[grouping+'_Human']),('log_depth',np.log10(meta[grouping+'_pg_parsimony']))]:
   context.append(dict(grouping=grouping,model=model,context=field,n=85,rho=spearmanr(val,vector).statistic))
  for typ,g in meta.assign(score=val).groupby('enterotype'):
   context.append(dict(grouping=grouping,model=model,context=typ,n=len(g),rho=spearmanr(g.ET_dysbiosis,g.score).statistic))
 # Adjusted descriptions of observed phenotype, separate from prediction performance.
 associations=[]
 for grouping in frames:
  cov=meta[['age','sex','antibiotic_broad']].copy();cov['host']=meta[grouping+'_Human'];cov['depth']=np.log10(meta[grouping+'_pg_parsimony'])
  for outcome in ['Human','Microbial','pg_parsimony']:
   for model,cs in [('unadjusted',[]),('age_sex_antibiotics',['age','sex','antibiotic_broad']),('mutual_host_depth',['age','sex','antibiotic_broad']+(['host'] if outcome=='pg_parsimony' else ['depth']))]:
    tab=cov[cs].assign(EDS=y,outcome=meta[grouping+'_'+outcome]).dropna()
    for c in tab:
     if c not in ['sex','antibiotic_broad']:tab[c]=(rankdata(tab[c])-.5)/len(tab)
    fit=sm.OLS(tab.outcome,sm.add_constant(tab.drop(columns='outcome'))).fit(cov_type='HC3')
    associations.append(dict(grouping=grouping,outcome=outcome,model=model,n=len(tab),beta=fit.params.EDS,p=fit.pvalues.EDS,rho=spearmanr(y,meta[grouping+'_'+outcome]).statistic))
 assoc=pd.DataFrame(associations);assoc['q_18']=multipletests(assoc.p,method='fdr_bh')[1]
 # Historical metrics retained separately from the new model/population.
 histrows=[]
 for scope,table in [('all_262',hist),('Figure4_85',hist.loc[meta.index])]:
  for model in ['og','depth_host_baseline']:
   v=table[model+'_predicted_ET_dysbiosis'];histrows.append(dict(scope=scope,model=model,n=len(table),patients=table.idpatient.nunique(),rho=spearmanr(table.ET_dysbiosis,v).statistic,r2=r2_score(table.ET_dysbiosis,v),MAE=mean_absolute_error(table.ET_dysbiosis,v)))
 write(meta.reset_index(drop=True),'patient_roster');write(metrics,'model_metrics');write(pd.concat(scores),'heldout_scores');write(pd.DataFrame(improvements),'increment_over_baseline');write(pd.DataFrame(grouping_check),'grouping_score_agreement');write(pd.DataFrame(context),'score_context');write(assoc,'phenotype_associations');write(pd.DataFrame(histrows),'historical_model_reproduction');write(pd.DataFrame(selections),'selected_alphas');write(pd.DataFrame(nulls),'permutation_results');write(pd.DataFrame(references),'distance_reference_folds');write(pd.DataFrame(FIT_AUDIT),'fit_boundary_audit')
 foldrows=[]
 for seed in SEEDS:
  for f,(tr,te) in enumerate(KFold(5,shuffle=True,random_state=seed).split(y)):
   for i in te:foldrows.append(dict(seed=seed,fold=f,row=i,specimen=meta.iloc[i].specimen,patient=int(meta.iloc[i].idpatient)))
 write(pd.DataFrame(foldrows),'outer_folds')
 report=dict(status='PASS',seconds=time.monotonic()-start,patients=85,groupings=list(frames),split_seeds=SEEDS,permutations=199,direct_sklearn_checks=checks,new_external_validation=False,scope='Fixed 85 patients; EDS-only target; cached transforms are label-independent and trained inside folds')
 (D/'COMPUTATION_VALIDATION.json').write_text(json.dumps(report,indent=2));(D/'INPUT_PROVENANCE.json').write_text(json.dumps(list(PROV.values()),indent=2))
 print(metrics[metrics.seed.eq(SEEDS[0])].to_string(index=False));print(pd.DataFrame(improvements).to_string(index=False));print('seconds',time.monotonic()-start,flush=True)
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('--pilot',action='store_true');main(parser.parse_args().pilot)
