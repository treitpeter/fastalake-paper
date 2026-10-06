"""M20: vary prevalence only, preserving production rosters, transformations and families."""
from pathlib import Path
from itertools import combinations
import sys,json,hashlib,os
import numpy as np,pandas as pd
from scipy.spatial.distance import pdist,squareform
from scipy.stats import rankdata
from statsmodels.stats.multitest import multipletests
O=Path(__file__).resolve().parent
ROOT=Path.cwd().resolve()
MODE=sys.argv[1]; assert MODE in ['pilot','full']
if MODE=='full':
 pilot=json.loads((O/'pilot/COMPLETE.json').read_text());assert pilot['status']=='PASS' and pilot['code_sha256']==hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
OUT=O/MODE; OUT.mkdir(exist_ok=True)
SEED=20260907; NPERM=0 if MODE=='pilot' else 1999
OLD=ROOT/'analyses/enterotypes_2026-09-07'
DATA=ROOT/'MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data'
CAP=ROOT/'fable_5.1_review_2026-09-16/scoping_review_2026-09-19/figure_comments_2026-09-21/capscan_dna/results'
AGG=ROOT/'analyses/capscan_metag_analysis_2026-09-10/aggregate'
NATIVE=ROOT/'00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23'
sources={};coverage=[];results=[];matrices={};rosters=[]
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rd(p,**kw):
 p=Path(p);sources[str(p)]=dict(sha256=sha(p),bytes=p.stat().st_size)
 return pd.read_csv(p,sep='\t',**kw)
def relative(d):return d.div(d.sum(axis=1).replace(0,np.nan),axis=0)
def logmet(d):return np.log10(d.where(d>0))
def fresh(co,lev):return rd(DATA/(co+'_noclust')/('host_representative_intensity.tsv' if lev=='host' else 'microbial_'+lev+'.tsv'),index_col=0)
def collapse(d):
 ids=d.index.to_series().str.extract(r'(MicrobPred_\d+)',expand=False)
 assert ids.notna().all()
 return d.groupby(ids).mean()
def transform(d,kind,prev,co,lev,scope):
 x=d.to_numpy(float)
 if kind=='metabolites':
  keep=np.isfinite(x).mean(0)>=.7
  y=x[:,keep];med=np.nanmedian(y,axis=0);y=np.where(np.isfinite(y),y,med[None,:])
  sd=y.std(0);y=y[:,sd>1e-8];sd=sd[sd>1e-8];y=(y-y.mean(0))/sd
  median=minimum=np.nan
 else:
  assert np.isfinite(x).all() and (x>=0).all()
  keep=(x>0).mean(0)>=prev if prev>0 else (x>0).any(0)
  y=x[:,keep];frac=y.sum(1)/x.sum(1);median=float(np.median(frac));minimum=float(frac.min())
  assert (y.sum(1)>0).all()
  y=np.sqrt(y/np.maximum(y.sum(1,keepdims=True),1e-30))
  for idx,v in zip(d.index,frac):rosters.append(dict(cohort=co,scope=scope,level=lev,layer=kind,prevalence=prev,unit=str(idx),signal_fraction=float(v)))
 assert np.isfinite(y).all()
 coverage.append(dict(cohort=co,scope=scope,level=lev,layer=kind,prevalence=prev,n_units=len(d),input_features=x.shape[1],retained_features=y.shape[1],median_signal_fraction=median,minimum_signal_fraction=minimum))
 return y
def evaluate(co,scope,lev,frames,ci,prev,pairset=None):
 mats={k:transform(v.loc[ci],k,prev,co,lev,scope) for k,v in frames.items()}
 n=len(ci);rng=np.random.default_rng(SEED)
 if co=='CAPSCAN':
  regions=np.tile(np.arange(5),n//5)
  design=np.column_stack([np.ones(n),pd.get_dummies(regions,drop_first=True,dtype=float)])
  mats={k:v-design @ np.linalg.lstsq(design,v,rcond=None)[0] for k,v in mats.items()}
  blocks=np.arange(n).reshape(n//5,5);perms=np.array([blocks[rng.permutation(n//5)].ravel() for _ in range(NPERM)])
 else:perms=np.array([rng.permutation(n) for _ in range(NPERM)])
 ix=np.triu_indices(n,1);ranks={}
 for k,v in mats.items():
  dist=pdist(v);r=rankdata(dist);r-=r.mean();r/=np.linalg.norm(r)
  ranks[k]=r
  matrices[co+'__'+scope+'__'+lev+'__'+str(prev)+'__'+k]=squareform(dist)
 for a,b in pairset or list(combinations(mats,2)):
  if lev!='OG' and {a,b}=={'host_metaP','metabolites'}:continue
  rho=float(ranks[a]@ranks[b]);p=np.nan
  if NPERM:
   bm=squareform(ranks[b]);ext=0
   for start in range(0,NPERM,200):
    pp=perms[start:start+200];null=bm[pp[:,ix[0]],pp[:,ix[1]]]@ranks[a];ext+=np.count_nonzero(np.abs(null)>=abs(rho))
   p=(1+ext)/(1+NPERM)
  results.append(dict(cohort=co,scope=scope,level=lev,layer_a=a,layer_b=b,prevalence=prev,n_units=n,n_donors=n if co=='PREDICT' else n//5,rho=rho,p=p))
first=rd(OLD/'PREDICT_fullstack_primary_samples.tsv',index_col=0);assert len(first)==77 and first.idpatient.nunique()==77
host=collapse(fresh('PREDICT','host'));met=logmet(rd(OLD/'PREDICT_metabolite_intensity.tsv',index_col=0))
for lev in ['OG','KO','Pathway']:
 frames=dict(metaG=rd(OLD/('metaG_OG_TPM.tsv' if lev=='OG' else 'PREDICT_metaG_'+lev+'.tsv'),index_col=0),metaT=rd(OLD/('metaT_OG_TPM.tsv' if lev=='OG' else 'PREDICT_metaT_'+lev+'.tsv'),index_col=0),microbial_metaP=collapse(fresh('PREDICT',lev)),host_metaP=host,metabolites=met)
 for prev in [0.,.1]:evaluate('PREDICT','displayed77',lev,frames,first.index,prev)
 print('PREDICT',lev,flush=True)
des=rd(OLD/'CAPSCAN_specimen_design.tsv',index_col=0);des['specimen']=des.specimen.astype(str)
met=rd(OLD/'CAPSCAN_metabolite_intensity.tsv',index_col=0);met.index=met.index.astype(str)
og=fresh('CAPSCAN','OG').reindex(des.index);valid=set(des.loc[og.sum(axis=1)>0,'specimen']);met=met.loc[met.index.isin(valid)]
dspec=des.groupby('specimen').first().loc[met.index]
host=relative(fresh('CAPSCAN','host').loc[des.index]).groupby(des.specimen).median()
cells={}
for lev in ['OG','KO','Pathway']:
 micro=relative(fresh('CAPSCAN',lev).loc[des.index]).groupby(des.specimen).median()
 frames=dict(microbial_metaP=micro,host_metaP=host,metabolites=logmet(met));cells[lev]={}
 for k,frame in frames.items():
  si=frame.index.intersection(dspec.index)
  cells[lev][k]=frame.loc[si].groupby([dspec.loc[si,'subject_id'],dspec.loc[si,'group']]).mean()
common=set.intersection(*(set(x.index) for fs in cells.values() for x in fs.values()))
reg=['duodenum','jejunum','ileum','colon','stool']
donors14=sorted(d for d in {c[0] for c in common} if all((d,r) in common for r in reg));assert len(donors14)==14
runs=rd(ROOT/'analyses/capscan_metag_analysis_2026-09-10/manifests/runs.tsv');runs['specimen']=runs.specimen.astype(str)
dna={}
for lev in ['OG','KO','Pathway']:
 d=rd(AGG/('metaG_'+lev+'.tsv'),index_col=0);d=d.loc[d.sum(axis=1)>0]
 d.index=runs.set_index('run_accession').loc[d.index,'specimen'].values
 d=d.groupby(level=0).mean();si=d.index.intersection(dspec.index)
 dna[lev]=relative(d.loc[si]).groupby([dspec.loc[si,'subject_id'],dspec.loc[si,'group']]).mean()
have=set.intersection(*(set(v.index) for v in dna.values()))
donors12=sorted(d for d in donors14 if all((d,r) in have for r in reg));assert len(donors12)==12
assert donors12==json.loads((CAP/'DNA_EXTENSION.json').read_text())['dna_complete_donors']
ci=lambda ds:pd.MultiIndex.from_tuples([(d,r) for d in ds for r in reg],names=['subject_id','group'])
for lev in ['OG','KO','Pathway']:
 for prev in [0.,.1]:
  evaluate('CAPSCAN','displayed_mixed12_14',lev,cells[lev],ci(donors14),prev)
  frames=dict(cells[lev],metaG=dna[lev])
  evaluate('CAPSCAN','displayed_mixed12_14',lev,frames,ci(donors12),prev,[('metaG',k) for k in ['microbial_metaP','host_metaP','metabolites']])
  frames=dict(metaG=dna[lev],**cells[lev])
  evaluate('CAPSCAN','matched12',lev,frames,ci(donors12),prev)
 print('CAPSCAN',lev,flush=True)
res=pd.DataFrame(results);res['q']=np.nan
for key,g in res.groupby(['cohort','scope','prevalence']):
 assert len(g)==(28 if key[0]=='PREDICT' else 16),(key,len(g))
 if NPERM:res.loc[g.index,'q']=multipletests(g.p,method='fdr_bh')[1]
refpred=rd(DATA/'crossomics/fullstack_profile_correlations.tsv');refpred=refpred[refpred.cohort.eq('PREDICT')]
refcap=rd(CAP/'capscan_profile_correlations_with_DNA.tsv');refcap=refcap[~refcap.status.str.startswith('same-donor')].copy();refcap['cohort']='CAPSCAN';refcap['q']=refcap.q_family
ref=pd.concat([refpred,refcap],ignore_index=True)
key=lambda r:(r['cohort'],r['level'],frozenset((r['layer_a'],r['layer_b'])))
refs={key(r):r for r in ref.to_dict('records')};replayed=res[(res.prevalence==.1)&res.scope.str.startswith('displayed')]
assert len(replayed)==len(refs)==44
diffs=[]
for r in replayed.to_dict('records'):
 b=refs[key(r)]
 diffs.append(dict(cohort=r['cohort'],level=r['level'],layer_a=r['layer_a'],layer_b=r['layer_b'],rho_difference=abs(r['rho']-b['rho']),p_difference=abs(r['p']-b['p']) if NPERM else None,q_difference=abs(r['q']-b['q']) if NPERM else None))
assert max(x['rho_difference'] for x in diffs)<1e-10
if NPERM:
 assert max(x['p_difference'] for x in diffs)<1e-10
 assert max(x['q_difference'] for x in diffs)<1e-10
res.to_csv(OUT/'functional_correlations.tsv',sep='\t',index=False)
pd.DataFrame(coverage).drop_duplicates().to_csv(OUT/'feature_coverage.tsv',sep='\t',index=False)
pd.DataFrame(rosters).drop_duplicates().to_csv(OUT/'signal_retention_by_profile.tsv',sep='\t',index=False)
pd.DataFrame(diffs).to_csv(OUT/'production_replay_differences.tsv',sep='\t',index=False)
np.savez_compressed(OUT/'distance_matrices.npz',**matrices)
for name in ['SOURCE_DATA_correlations.tsv','SOURCE_DATA_feature_summary.tsv']:
 tab=rd(NATIVE/name);tab.to_csv(OUT/('native_'+name),sep='\t',index=False)
record=dict(status='PASS',mode=MODE,permutations=NPERM,seed=SEED,production_endpoints_reproduced=44,functional_result_rows=len(res),sources=sources,code_sha256=sha(__file__),job=os.environ.get('SLURM_JOB_ID'),scope='Functional prevalence-only sensitivity and CapScan population check; native results reused from verified saved analyses. No raw search or gene mapping rerun.')
(OUT/'COMPLETE.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:v for k,v in record.items() if k!='sources'}),flush=True)
