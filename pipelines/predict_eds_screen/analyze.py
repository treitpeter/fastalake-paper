"""Current-manuscript, patient-matched exploratory reaction/compound screen."""
from pathlib import Path
import json,hashlib,sys,os
import numpy as np
import pandas as pd
from scipy.stats import rankdata,spearmanr,t
H=Path(__file__).resolve().parent
ROOT=None  # unused in the original
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def read(key,**kw):return pd.read_csv(inputs[key],sep='\t',**kw)
def bh(p):
 p=np.asarray(p,float);v=np.full(len(p),np.nan);idx=np.where(np.isfinite(p))[0];s=idx[np.argsort(p[idx])];v[s]=np.minimum.accumulate((p[s]*len(s)/np.arange(1,len(s)+1))[::-1])[::-1].clip(0,1);return v
def corr(x,y,c=None):
 ok=np.isfinite(x)&np.isfinite(y)
 if c is not None:ok &= np.isfinite(c).all(axis=1)
 n=int(ok.sum())
 if n<43 or np.unique(x[ok]).size<3 or np.unique(y[ok]).size<3:return n,np.nan,np.nan
 a=rankdata(x[ok]);b=rankdata(y[ok]);k=0
 if c is not None:
  C=np.column_stack([np.ones(n),c[ok]]);k=np.linalg.matrix_rank(C)-1
  a=a-C @ np.linalg.lstsq(C,a,rcond=None)[0];b=b-C @ np.linalg.lstsq(C,b,rcond=None)[0]
 r=np.corrcoef(a,b)[0,1];df=n-k-2
 return n,float(r),float(2*t.sf(abs(r)*np.sqrt(df/max(1e-16,1-r*r)),df))
manifest=json.loads((H/'INPUTS.json').read_text());inputs={k:Path(v['path']) for k,v in manifest.items()}
for k,v in manifest.items():assert sha(v['path'])==v['sha256'],k
order=read('order');assert len(order)==order.patient_id.nunique()==85
z=np.load(inputs['heatmap']);features=read('features');order=order.set_index('sample_id').loc[z['samples']].reset_index();eds=order.dysbiosis_score.to_numpy()
meta=read('metabolites',index_col=0).loc[order.microbiome_id];meta.index=order.sample_id
ann=read('identity',keep_default_na=False).set_index('feature');native=read('annotation',keep_default_na=False).set_index('variableid')
assert all(ann.loc[f,'name']==native.loc[f,'annotation'] for f in meta.columns)
mat={level:read(level,index_col=0).drop(index='FunctionallyUnassigned',errors='ignore').loc[:,order.sample_id].T for level in ['KO','Pathway']}
for j,f in features.iterrows():
 v=z['raw'][j];positive=(v>0)&np.isfinite(v)
 if positive.sum()>=3:assert np.isclose(spearmanr(v[positive],eds[positive]).statistic,f.rho_EDS_85,atol=1e-12)
 if f.block=='Microbial KEGG pathway':assert np.allclose(v,mat['Pathway'][f.feature_id])
 if f.block=='Metabolomics':assert np.allclose(v,meta[f.feature_id])
assert all(x.index.equals(meta.index) for x in mat.values())
C=np.column_stack([rankdata(order.Human),rankdata(np.log2(order.pg_parsimony)),pd.get_dummies(order[['antibiotic_state','clinical_group']],drop_first=True).to_numpy(float)])
cr=read('compound_reaction',header=None,names=['compound','reaction']);rk=read('reaction_ko',header=None,names=['reaction','KO']);cp=read('compound_pathway',header=None,names=['compound','Pathway'])
for frame in [cr,rk,cp]:
 for col in frame:frame[col]=frame[col].str.split(':').str[-1]
cp.Pathway=cp.Pathway.str.replace('map','map',regex=False)
a=ann.reset_index();a=a[a.feature.isin(meta.columns)&a.compound.str.fullmatch('C[0-9]{5}')]
links=a[['feature','name','compound']].merge(cr,on='compound').merge(rk,on='reaction');links['measured_KO']=links.KO.isin(mat['KO'].columns)
# Keep absent targets and every reaction link; no favorable-direction filtering.
links.to_csv(H/'reaction_links.tsv',sep='\t',index=False)
if sys.argv[1]=='pilot':
 n,r,p=corr(np.arange(85),np.arange(85));assert n==85 and r>.99999
 assert np.allclose(bh([.01,.04,.03]),[.03,.04,.04])
 save(H/'PILOT.json',dict(status='PASS',patients=85,heatmap_features_reproduced=20,metabolites=len(meta.columns),KO=len(mat['KO'].columns),pathways=len(mat['Pathway'].columns),reaction_links=len(links),script_sha256=sha(__file__)))
 print('Pilot PASS',flush=True);sys.exit()
assert json.loads((H/'PILOT.json').read_text())['script_sha256']==sha(__file__)
rows=[];pairrows=[]
for level,frame in [('Metabolite',meta),*mat.items()]:
 for norm in (['native'] if level=='Metabolite' else ['native','annotated_relative']):
  values=frame.to_numpy(float).copy()
  if norm=='annotated_relative':values/=values.sum(axis=1)[:,None]
  values=np.where(values>0,values,np.nan)
  for j,f in enumerate(frame.columns):
   x=values[:,j]
   for model,c in [('unadjusted',None),('host_depth_antibiotic_clinical',C)]:
    n,r,p=corr(x,eds,c);rows.append(dict(layer=level,normalization=norm,model=model,feature=f,n=n,rho=r,p=p))
  if level=='Metabolite':continue
  l=links[links.measured_KO] if level=='KO' else a[['feature','name','compound']].merge(cp,on='compound')
  l=l[l[level].isin(frame.columns)].drop_duplicates(['feature',level])
  for _,link in l.iterrows():
   x=values[:,frame.columns.get_loc(link[level])];y=meta[link.feature].to_numpy(float);y=np.where(y>0,y,np.nan)
   for model,c in [('unadjusted',None),('host_depth_antibiotic_clinical',C),('host_depth_antibiotic_clinical_EDS',np.column_stack([C,rankdata(eds)]))]:
    n,r,p=corr(x,y,c);pairrows.append(dict(layer=level,normalization=norm,model=model,metabolite=link.feature,name=link['name'],compound=link.compound,protein_feature=link[level],n=n,rho=r,p=p))
results=pd.DataFrame(rows);pairs=pd.DataFrame(pairrows)
for df in [results,pairs]:df['q']=df.groupby(['layer','normalization','model']).p.transform(lambda v:bh(v))
results.to_csv(H/'all_EDS_associations.tsv',sep='\t',index=False);pairs.to_csv(H/'all_linked_associations.tsv',sep='\t',index=False)
joined=pairs.merge(results[results.layer=='Metabolite'].drop(columns=['layer','normalization']).rename(columns={'feature':'metabolite','rho':'metabolite_EDS_rho','q':'metabolite_EDS_q','p':'metabolite_EDS_p','n':'metabolite_EDS_n'}),on=['metabolite','model'],how='left').merge(results.rename(columns={'feature':'protein_feature','rho':'protein_EDS_rho','q':'protein_EDS_q','p':'protein_EDS_p','n':'protein_EDS_n'}),on=['layer','normalization','model','protein_feature'],how='left')
joined['opposing_EDS_directions']=(joined.metabolite_EDS_rho>0)&(joined.protein_EDS_rho<0)
joined['both_EDS_q05_and_inverse_pair_q05']=joined.opposing_EDS_directions & (joined.metabolite_EDS_q<.05)&(joined.protein_EDS_q<.05)&(joined.rho<0)&(joined.q<.05)
joined.to_csv(H/'joined_evidence.tsv',sep='\t',index=False)
selected=features[features.block=='Metabolomics'].feature_id
joined[joined.metabolite.isin(selected)].to_csv(H/'heatmap_metabolite_linked_evidence.tsv',sep='\t',index=False)
save(H/'COMPLETE.json',dict(status='COMPUTATION_COMPLETE_REQUIRES_REACTION_DIRECTION_AND_SCIENTIFIC_REVIEW',patients=85,EDS_tests=len(results),pair_tests=len(pairs),reaction_links=len(links),job_id=os.environ.get('SLURM_JOB_ID'),script_sha256=sha(__file__),scope='Continuous patient-matched exploratory screen, whole measured feature families, two-sided tests; selected display is not independent validation. Reaction membership does not establish consumption direction or flux.'))
print('Analysis complete',flush=True)
