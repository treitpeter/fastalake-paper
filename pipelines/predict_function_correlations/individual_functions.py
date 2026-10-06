from common import *
from itertools import combinations
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests
verify();assert sha(__file__)==load(H/'INDIVIDUAL_CODE.json')['individual_functions.py']
for p,h in load(H/'INDIVIDUAL_INPUTS.json').items():assert sha(p)==h,p
first=rd(OLD/'PREDICT_fullstack_primary_samples.tsv',index_col=0);assert len(first)==first.idpatient.nunique()==77
D=A33/'results/crossomics_repaired__simplified/source_data'
def relative(d):return d.div(d.sum(1).replace(0,np.nan),axis=0)
def collapse(d):
 ids=d.index.to_series().str.extract(r'(MicrobPred_\d+)',expand=False);assert ids.notna().all();return d.groupby(ids).mean()
for mode in ['baseline','refit']:
 out=H/('individual_'+mode);out.mkdir(exist_ok=False);rows=[];elig=[]
 for lev in ['OG','KO','Pathway']:
  raw=rd(D/'PREDICT_noclust'/('microbial_'+lev+'.tsv'),index_col=0) if mode=='baseline' else frame(C/'profiles/PREDICT_noclust',lev+'_other_unresolved')
  fs=dict(metaG=rd(OLD/('metaG_OG_TPM.tsv' if lev=='OG' else 'PREDICT_metaG_'+lev+'.tsv'),index_col=0),metaT=rd(OLD/('metaT_OG_TPM.tsv' if lev=='OG' else 'PREDICT_metaT_'+lev+'.tsv'),index_col=0),microbial_metaP=collapse(raw))
  for a,b in combinations(fs,2):
   terms=sorted(set(fs[a])&set(fs[b]));aa,bb=[relative(fs[k].loc[first.index,terms]) for k in [a,b]]
   for term in terms:
    ok=(aa[term]>0)&(bb[term]>0);n=int(ok.sum());r=p=np.nan
    if n>=35:r,p=spearmanr(aa.loc[ok,term],bb.loc[ok,term])
    elig.append(dict(level=lev,layer_a=a,layer_b=b,feature=term,n_patients=n,eligible=bool(np.isfinite(p))))
    if np.isfinite(p):rows.append(dict(level=lev,layer_a=a,layer_b=b,feature=term,n_patients=n,rho=r,p=p))
  print(mode,lev,flush=True)
 res=pd.DataFrame(rows);res['q']=res.groupby(['level','layer_a','layer_b']).p.transform(lambda v:multipletests(v,method='fdr_bh')[1]);wr(res,out/'crossomics_same_function_correlations.tsv');wr(pd.DataFrame(elig),out/'eligibility.tsv')
 summary=res.groupby(['level','layer_a','layer_b']).agg(features=('rho','size'),median_rho=('rho','median'),significant=('q',lambda v:int((v<.05).sum()))).reset_index();wr(summary,out/'same_function_summary.tsv')
 old=rd(D/'crossomics/crossomics_same_function_correlations.tsv' if mode=='baseline' else H/'individual_baseline/crossomics_same_function_correlations.tsv');keys=['level','layer_a','layer_b','feature'];comp=old.merge(res,on=keys,how='outer',suffixes=('_previous','_revised'),indicator=True,validate='one_to_one');wr(comp,out/'method_comparison.tsv')
 if mode=='baseline':
  assert len(old)==len(res)==17914 and comp._merge.eq('both').all()
  for col in ['rho','p','q','n_patients']:assert np.allclose(comp[col+'_previous'],comp[col+'_revised'],atol=1e-10),(col,abs(comp[col+'_previous']-comp[col+'_revised']).max())
 else:
  control=comp[comp.layer_a.eq('metaG')&comp.layer_b.eq('metaT')];assert control._merge.eq('both').all()
  for col in ['rho','p','q']:assert np.allclose(control[col+'_previous'],control[col+'_revised'],atol=1e-10)
 receipt(out,mode=mode,tests=len(res),patients=77,full_BH_families=9,driver_sha256=sha(__file__),normalization='Relative signal within pairwise shared functional vocabulary; positive in both layers in at least35patients')
