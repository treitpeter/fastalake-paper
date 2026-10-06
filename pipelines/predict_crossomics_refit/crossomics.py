from common import *
from itertools import combinations
from scipy.spatial.distance import pdist,squareform
from scipy.stats import rankdata,spearmanr
from statsmodels.stats.multitest import multipletests
import math,importlib.util
mode=sys.argv[1];assert mode in ['baseline','refit'];verify()
for p,h in load(H/'CROSS_INPUTS.json').items():assert sha(p)==h,p
assert sha(__file__)==load(H/'CROSS_CODE.json')['crossomics.py']
if mode=='refit':gate(H/'cross_baseline');gate(H/'groups');gate(H/'groups/PREDICT');gate(H/'groups/CAPSCAN')
out=H/('cross_'+mode);out.mkdir(exist_ok=False)
DATA=A33/'results/crossomics_repaired__simplified/source_data'
SEED=20260907;NPERM=1999;coverage=[];results=[];matrices={};rosters=[]
functions(X/'provenance/check_prevalence.py',['relative','logmet','transform','evaluate'],globals())
# The legacy mixed12/14 archive reused a key for60and70cells. Statistics were
# correct; persist BOTH matrices with explicit sample size to prevent overwriting.
_original_evaluate=evaluate
def evaluate(co,scope,lev,frames,ci,prev,pairset=None):
 _original_evaluate(co,scope,lev,frames,ci,prev,pairset)
 prefix='__'.join([co,scope,lev,str(prev)])
 for layer in frames:
  matrices[prefix+'__n'+str(len(ci))+'__'+layer]=matrices.pop(prefix+'__'+layer)

ns=dict(np=np,pd=pd,sparse=sparse,math=math,SEED=SEED,combinations=combinations,rankdata=rankdata,squareform=squareform,spearmanr=spearmanr)
functions(X/'native/native_pyramids.py',['profile_distance','bh','correlations'],ns)
def fresh(co,lev):
 if mode=='baseline':return rd(DATA/(co+'_noclust')/('host_representative_intensity.tsv' if lev=='host' else 'microbial_'+lev+'.tsv'),index_col=0)
 return frame(H/'groups'/co,'supported_human_groups') if lev=='host' else frame(C/'profiles'/(co+'_noclust'),lev+'_other_unresolved')
def collapse(d):
 ids=d.index.to_series().str.extract(r'(MicrobPred_\d+)',expand=False);assert ids.notna().all()
 return d.groupby(ids).mean()
first=rd(OLD/'PREDICT_fullstack_primary_samples.tsv',index_col=0);assert len(first)==first.idpatient.nunique()==77
host=collapse(fresh('PREDICT','host'));met=logmet(rd(OLD/'PREDICT_metabolite_intensity.tsv',index_col=0));predictframes={}
for lev in ['OG','KO','Pathway']:
 fs=dict(metaG=rd(OLD/('metaG_OG_TPM.tsv' if lev=='OG' else 'PREDICT_metaG_'+lev+'.tsv'),index_col=0),metaT=rd(OLD/('metaT_OG_TPM.tsv' if lev=='OG' else 'PREDICT_metaT_'+lev+'.tsv'),index_col=0),microbial_metaP=collapse(fresh('PREDICT',lev)),host_metaP=host,metabolites=met)
 predictframes[lev]=fs
 for prev in [0.,.1]:evaluate('PREDICT','displayed77',lev,fs,first.index,prev)
 print(mode,'functional PREDICT',lev,flush=True)
des=rd(OLD/'CAPSCAN_specimen_design.tsv',index_col=0);des.specimen=des.specimen.astype(str)
met=rd(OLD/'CAPSCAN_metabolite_intensity.tsv',index_col=0);met.index=met.index.astype(str)
oldog=rd(DATA/'CAPSCAN_noclust/microbial_OG.tsv',index_col=0).reindex(des.index)
newog=fresh('CAPSCAN','OG').reindex(des.index)
wr(pd.DataFrame(dict(acquisition=des.index,previous_positive=oldog.sum(1).gt(0).to_numpy(),current_positive=newog.sum(1).gt(0).to_numpy())),out/'eligibility_audit.tsv')
# Preserve the exact prior contributing-specimen roster for the method comparison.
valid=set(des.loc[oldog.sum(1)>0,'specimen']);met=met.loc[met.index.isin(valid)]
dspec=des.groupby('specimen').first().loc[met.index];host=relative(fresh('CAPSCAN','host').loc[des.index]).groupby(des.specimen).median();cells={}
for lev in ['OG','KO','Pathway']:
 micro=relative(fresh('CAPSCAN',lev).loc[des.index]).groupby(des.specimen).median();cells[lev]={}
 for k,f in dict(microbial_metaP=micro,host_metaP=host,metabolites=logmet(met)).items():
  si=f.index.intersection(dspec.index);cells[lev][k]=f.loc[si].groupby([dspec.loc[si,'subject_id'],dspec.loc[si,'group']]).mean()
common=set.intersection(*(set(x.index) for fs in cells.values() for x in fs.values()))
donors14=sorted(d for d in {c[0] for c in common} if all((d,r) in common for r in REG));assert len(donors14)==14
runs=rd(R/'analyses/capscan_metag_analysis_2026-09-10/manifests/runs.tsv');runs.specimen=runs.specimen.astype(str);dna={}
for lev in ['OG','KO','Pathway']:
 d=rd(R/'analyses/capscan_metag_analysis_2026-09-10/aggregate'/('metaG_'+lev+'.tsv'),index_col=0);d=d.loc[d.sum(1)>0];d.index=runs.set_index('run_accession').loc[d.index,'specimen'].values;d=d.groupby(level=0).mean();si=d.index.intersection(dspec.index);dna[lev]=relative(d.loc[si]).groupby([dspec.loc[si,'subject_id'],dspec.loc[si,'group']]).mean()
ci=lambda ds:pd.MultiIndex.from_product([ds,REG],names=['subject_id','group'])
assert all(set(ci(DONORS))<=set(v.index) for v in dna.values())
for lev in ['OG','KO','Pathway']:
 for prev in [0.,.1]:
  evaluate('CAPSCAN','displayed_mixed12_14',lev,cells[lev],ci(donors14),prev)
  evaluate('CAPSCAN','displayed_mixed12_14',lev,dict(cells[lev],metaG=dna[lev]),ci(DONORS),prev,[('metaG',k) for k in ['microbial_metaP','host_metaP','metabolites']])
  evaluate('CAPSCAN','matched12',lev,dict(metaG=dna[lev],**cells[lev]),ci(DONORS),prev)
 print(mode,'functional CAPSCAN',lev,flush=True)
res=pd.DataFrame(results);res['q']=np.nan
for key,g in res.groupby(['cohort','scope','prevalence']):
 assert len(g)==(28 if key[0]=='PREDICT' else 16);res.loc[g.index,'q']=multipletests(g.p,method='fdr_bh')[1]
wr(res,out/'functional_correlations.tsv');wr(pd.DataFrame(coverage).drop_duplicates(),out/'feature_coverage.tsv');wr(pd.DataFrame(rosters).drop_duplicates(),out/'signal_retention.tsv');np.savez_compressed(out/'functional_distances.npz',**matrices)
wr(dspec.reset_index(),out/'capscan_specimen_roster.tsv');wr(first.reset_index(),out/'predict_patient_roster.tsv')
# Native DNA/RNA sequence axes are unchanged; distances come from hash-verified adopted sources.
nativerows=[];nativemats={};qualities={};controlerrs={}
for co in ['PREDICT','CAPSCAN']:
 native=np.load(X/('native_'+co+'_distance_matrices.npz'))
 if mode=='baseline':
  m=rd(groupdir(co)/'reaggregated/simplified_study_group_matrix.tsv.gz',index_col=0).fillna(0).T
  hg=m.loc[:,m.columns.str.contains('_HUMAN',regex=False)];ng=m.loc[:,~m.columns.str.contains('_HUMAN',regex=False)]
  assert hg.equals(fresh(co,'host')) or np.allclose(hg,fresh(co,'host'))
 else:hg=frame(H/'groups'/co,'supported_human_groups');ng=frame(H/'groups'/co,'other_unresolved_groups')
 if co=='PREDICT':
  fs={'Non-host protein':collapse(ng).loc[first.index],'Host protein':collapse(hg).loc[first.index],'Metabolites':logmet(rd(OLD/'PREDICT_metabolite_intensity.tsv',index_col=0)).loc[first.index]}
 else:
  fs={}
  for k,f in [('Non-host protein',ng),('Host protein',hg)]:
   byspec=relative(f.loc[des.index]).groupby(des.specimen).median();si=byspec.index.intersection(dspec.index)
   fs[k]=byspec.loc[si].groupby([dspec.loc[si,'subject_id'],dspec.loc[si,'group']]).mean().loc[ci(DONORS)]
  fs['Metabolites']=logmet(met).groupby([dspec.subject_id,dspec.group]).mean().loc[ci(DONORS)]
 for label,prev in [('primary_prevalence10',.1),('sensitivity_all_features',0.)]:
  ds={k:native[label+'__'+k] for k in (['DNA','RNA'] if co=='PREDICT' else ['DNA'])};qualities[co+'__'+label]={}
  for k,f in fs.items():
   ds[k],q=ns['profile_distance'](f.to_numpy() if k=='Metabolites' else sparse.csr_matrix(f.to_numpy()),k,prev,capscan=co=='CAPSCAN');qualities[co+'__'+label][k]=q
   if mode=='baseline' or k=='Metabolites':
    err=float(np.max(abs(ds[k]-native[label+'__'+k])));assert err<1e-7,(co,k,err);controlerrs[co+'__'+label+'__'+k]=err
  tab=ns['correlations'](ds,co);tab.insert(1,'analysis',label);nativerows.append(tab)
  nativemats.update({co+'__'+label+'__'+k:v for k,v in ds.items()})
 print(mode,'native',co,flush=True)
native=pd.concat(nativerows,ignore_index=True);wr(native,out/'native_correlations.tsv');np.savez_compressed(out/'native_distances.npz',**nativemats);dump(out/'native_quality.json',qualities)
checks={}
for name,oldpath,keys in [('functional',X/'full/functional_correlations.tsv',['cohort','scope','level','layer_a','layer_b','prevalence']),('native',X/'full/native_SOURCE_DATA_correlations.tsv',['cohort','analysis','layer_a','layer_b'])]:
 old=rd(oldpath if mode=='baseline' else H/'cross_baseline'/(name+'_correlations.tsv'));got=res if name=='functional' else native
 merged=old.merge(got,on=keys,suffixes=('_previous','_revised'),validate='one_to_one');assert len(merged)==len(old)==len(got)
 for c in ['rho','p','q']:
  err=float(abs(merged[c+'_previous']-merged[c+'_revised']).max());checks[name+'_'+c]=err
  if mode=='baseline':assert err<1e-10,(name,c,err)
 wr(merged,out/(name+'_method_comparison.tsv'))
if mode=='refit':
 for f in ['capscan_specimen_roster.tsv','predict_patient_roster.tsv']:assert rd(out/f).equals(rd(H/'cross_baseline'/f)),f
 for fn in ['functional_distances.npz','native_distances.npz']:
  old=np.load(H/'cross_baseline'/fn);new=np.load(out/fn)
  for k in old.files:
   if k.split('__')[-1] in ['metaG','metaT','metabolites','DNA','RNA','Metabolites']:assert np.allclose(old[k],new[k],rtol=1e-10,atol=1e-10),(fn,k)
receipt(out,mode=mode,functional_rows=len(res),native_rows=len(native),baseline_errors=checks,native_controls=controlerrs,primary_prevalence=0,capscan_labels='PROVISIONAL_FROZEN_ROSTER_KEY_MISSING',driver_sha256=sha(__file__))
