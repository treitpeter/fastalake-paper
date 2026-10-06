from pathlib import Path
import ast,json,hashlib,os
import numpy as np,pandas as pd
from scipy.stats import rankdata,spearmanr,t
from statsmodels.stats.multitest import multipletests
H=Path(__file__).resolve().parent
# Reuse frozen statistical definitions only, without rerunning the original screen.
source=ast.parse((H/'analyze.py').read_text());defs=[n for n in source.body if isinstance(n,ast.FunctionDef) and n.name in ['corr','bh']];exec(compile(ast.Module(body=defs,type_ignores=[]),'<frozen statistics>','exec'))
m=json.loads((H/'INPUTS.json').read_text());p={k:Path(v['path']) for k,v in m.items()}
for k,v in m.items():assert hashlib.sha256(p[k].read_bytes()).hexdigest()==v['sha256']
order=pd.read_csv(p['order'],sep='\t');z=np.load(p['heatmap']);order=order.set_index('sample_id').loc[z['samples']].reset_index();eds=order.dysbiosis_score.to_numpy()
meta=pd.read_csv(p['metabolites'],sep='\t',index_col=0).loc[order.microbiome_id];ko=pd.read_csv(p['KO'],sep='\t',index_col=0).drop(index='FunctionallyUnassigned').loc[:,order.sample_id].T
C=np.column_stack([rankdata(order.Human),rankdata(np.log2(order.pg_parsimony)),pd.get_dummies(order[['antibiotic_state','clinical_group']],drop_first=True).to_numpy(float)])
links=json.loads((H/'ADDITIONAL_LINKS.json').read_text());rows=[]
for link in links:
 for norm in ['native','annotated_relative']:
  if link['KO'] in ko:
   x=ko[link['KO']].to_numpy(float);x=x if norm=='native' else x/ko.sum(axis=1).to_numpy();x=np.where(x>0,x,np.nan)
  else:x=np.full(85,np.nan)
  y=meta[link['metabolite']].to_numpy(float);y=np.where(y>0,y,np.nan)
  for model,c in [('unadjusted',None),('host_depth_antibiotic_clinical',C),('host_depth_antibiotic_clinical_EDS',np.column_stack([C,rankdata(eds)]))]:
   n,r,pairp=corr(x,y,c);rows.append(dict(**link,layer='KO',normalization=norm,model=model,protein_feature=link['KO'],n=n,rho=r,p=pairp))
new=pd.DataFrame(rows);prior=pd.read_csv(H/'all_linked_associations.tsv',sep='\t');allpairs=pd.concat([prior,new],ignore_index=True);assert not allpairs.duplicated(['layer','normalization','model','metabolite','protein_feature']).any()
allpairs['q_augmented']=allpairs.groupby(['layer','normalization','model']).p.transform(lambda v:bh(v))
for key,group in allpairs.groupby(['layer','normalization','model']):
 good=group.p.notna();assert np.allclose(group.loc[good,'q_augmented'],multipletests(group.loc[good,'p'],method='fdr_bh')[1])
allpairs.to_csv(H/'augmented_pair_screen.tsv',sep='\t',index=False);allpairs[allpairs.source.notna()].to_csv(H/'individual_followup_results.tsv',sep='\t',index=False)
# Independent numerical spot validation of all heatmap correlations and partial rank residual formula.
features=pd.read_csv(p['features'],sep='\t');checks=[]
for feature in features[features.block=='Metabolomics'].feature_id:
 y=meta[feature].to_numpy(float);valid=np.isfinite(y)&(y>0)
 rho,pv=spearmanr(y[valid],eds[valid]);checks.append(dict(feature=feature,rho=float(rho)))
 n,r,pv=corr(y,eds,C)
 Z=np.column_stack([np.ones(len(C)),C]);Q=np.linalg.qr(Z)[0];xr=rankdata(y)-Q@(Q.T@rankdata(y));yr=rankdata(eds)-Q@(Q.T@rankdata(eds));assert np.isclose(r,np.corrcoef(xr,yr)[0,1],atol=1e-10)
# Bootstrap uncertainties for one literature-supported consuming enzyme and its two substrates.
rng=np.random.default_rng(2026093031);ci=[]
for feature in ['C18_F1906','C18_F1073']:
 if feature not in meta:continue
 x=ko['K01667'].to_numpy(float);y=meta[feature].to_numpy(float);valid=(x>0)&(y>0)&np.isfinite(x)&np.isfinite(y);xx=x[valid];yy=y[valid];n=len(xx);boot=[]
 for _ in range(2000):
  idx=rng.integers(0,n,n);boot.append(spearmanr(xx[idx],yy[idx]).statistic)
 lo,hi=np.quantile(boot,[.025,.975]);ci.append(dict(feature=feature,KO='K01667',n=n,rho=float(spearmanr(xx,yy).statistic),lower95=float(lo),upper95=float(hi)))
(H/'FOLLOWUP_COMPLETE.json').write_text(json.dumps(dict(status='PASS_NUMERICAL_REQUIRES_BIOCHEMICAL_REVIEW',added_links=links,bootstrap=ci,independent_BH_verified=True,partial_rank_QR_verified=True,job_id=os.environ.get('SLURM_JOB_ID'),scope='Literature and stereochemistry follow-ups are explicitly post hoc; q values use the augmented entire connected-pair family. Host-depth adjustment and relative normalization are sensitivities, not a causal identification model.'),indent=2))
print(new[['metabolite','KO','normalization','model','n','rho','p']].to_string(index=False))
