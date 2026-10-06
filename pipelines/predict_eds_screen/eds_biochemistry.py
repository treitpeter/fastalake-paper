from common import *
from scipy.stats import spearmanr
import shutil,subprocess
mode=sys.argv[1];assert mode in ['baseline','refit'];verify()
for p,h in load(H/'EDS_INPUTS.json').items():assert sha(p)==h,p
assert sha(__file__)==load(H/'EDS_CODE.json')['eds_biochemistry.py']
if mode=='refit':gate(H/'eds_baseline');gate(H/'groups/PREDICT')
out=H/('eds_'+mode);out.mkdir(exist_ok=False)
manifest=load(E/'INPUTS.json')
if mode=='refit':
 data=out/'input_data';data.mkdir();old=Path(manifest['heatmap']['path']).parent
 feat=rd(old/'heatmap_features.tsv');cols=rd(old/'heatmap_columns.tsv');order=rd(old/'heatmap_85_patient_order.tsv');raw=np.load(old/'heatmap_values.npz')['raw'].copy()
 assert raw.shape==(20,259)
 cov=rd(C/'clinical_refit/acquisition_accounting.tsv').set_index('sample_id')
 for d in [cols,order]:
  d['Human']=cov.loc[d.sample_id,'Human'].to_numpy();d['pg_parsimony']=cov.loc[d.sample_id,'pg_parsimony'].to_numpy()
 hp=frame(H/'groups/PREDICT','supported_human_groups');path=frame(C/'profiles/PREDICT_noclust','Pathway_other_unresolved')
 for j,r in feat.iterrows():
  if r.block=='Host proteome':raw[j]=hp[r.feature_id].reindex(cols.sample_id).to_numpy() if r.feature_id in hp else 0
  elif r.block=='Microbial KEGG pathway':raw[j]=path[r.feature_id].loc[cols.sample_id].to_numpy()
 ix=pd.Index(cols.sample_id).get_indexer(order.sample_id);assert (ix>=0).all();val=raw[:,ix]
 pos=np.where((val>0)&np.isfinite(val),val,np.nan)
 with np.errstate(invalid='ignore',divide='ignore'):
  zz=(pos-np.nanmean(pos,axis=1)[:,None])/np.nanstd(pos,axis=1)[:,None]
 feat['n_positive_85']=np.isfinite(pos).sum(1);feat['rho_EDS_85']=[spearmanr(v[np.isfinite(v)],order.dysbiosis_score.to_numpy()[np.isfinite(v)]).statistic if np.isfinite(v).sum()>=3 else np.nan for v in pos]
 wr(feat,data/'heatmap_features.tsv');wr(cols,data/'heatmap_columns.tsv');wr(order,data/'heatmap_85_patient_order.tsv')
 np.savez_compressed(data/'heatmap_values.npz',raw=raw,samples=cols.sample_id.to_numpy(dtype=str))
 np.savez_compressed(data/'heatmap_85_values.npz',raw=val,z=zz,features=feat.feature_id.to_numpy(dtype=str),samples=order.sample_id.to_numpy(dtype=str))
 for key,name in [('heatmap','heatmap_85_values.npz'),('features','heatmap_features.tsv'),('order','heatmap_85_patient_order.tsv')]:manifest[key]=dict(path=str(data/name),sha256=sha(data/name))
 for lev in ['KO','Pathway']:
  # Original driver expects feature-by-acquisition matrices.
  p=data/(lev+'_profile.tsv.gz');adapter=frame(C/'profiles/PREDICT_noclust',lev+'_other_unresolved').T;adapter.loc['FunctionallyUnassigned']=0.;adapter.to_csv(p,sep='\t');manifest[lev]=dict(path=str(p),sha256=sha(p))
 dump(data/'HEATMAP_CHECK.json',dict(patients=len(order),observations=len(cols),fixed_features=len(feat),positive_features=int((feat.n_positive_85>0).sum()),no_reselection=True,host_quantity='Human-assigned peptide signal within stable simplified reporting-group coordinates',metabolite_values_unchanged=bool(np.array_equal(raw[feat.block.eq('Metabolomics')],np.load(old/'heatmap_values.npz')['raw'][feat.block.eq('Metabolomics')]))))
for name in ['analyze.py','individual_followup.py','ADDITIONAL_LINKS.json']:shutil.copy2(E/name,out/name)
dump(out/'INPUTS.json',manifest)
for args in [['analyze.py','pilot'],['analyze.py','full'],['individual_followup.py']]:subprocess.run([sys.executable,str(out/args[0]),*args[1:]],check=True)
sep=out/'separate_adjustments';sep.mkdir()
for name in ['run_analysis.py','TARGETS.json']:shutil.copy2(E/'separate_adjustments'/name,sep/name)
dump(sep/'FREEZE.json',{'files':[dict(path=str(sep/n),sha256=sha(sep/n)) for n in ['run_analysis.py','TARGETS.json']]})
for phase in ['pilot','full']:subprocess.run([sys.executable,str(sep/'run_analysis.py'),phase],check=True)
checks={}
if mode=='baseline':
 for f in ['all_EDS_associations.tsv','all_linked_associations.tsv','augmented_pair_screen.tsv','separate_adjustments/all_models.tsv']:
  a=rd(out/f);b=rd(E/f);assert a.shape==b.shape and list(a)==list(b),f
  for col in a:
   if pd.api.types.is_numeric_dtype(b[col]):assert np.allclose(a[col],b[col],atol=1e-9,rtol=1e-8,equal_nan=True),(f,col)
   else:assert a[col].fillna('').equals(b[col].fillna('')),(f,col)
  checks[f]=len(a)
else:
 for f,keys in [('all_EDS_associations.tsv',['layer','normalization','model','feature']),('augmented_pair_screen.tsv',['layer','normalization','model','metabolite','protein_feature']),('separate_adjustments/all_models.tsv',['id','scope','normalization','model'])]:
  a=rd(out/f);b=rd(H/'eds_baseline'/f);m=b.merge(a,on=keys,how='outer',suffixes=('_previous','_revised'),indicator=True,validate='one_to_one');wr(m,out/(Path(f).stem+'_comparison.tsv'))
receipt(out,mode=mode,baseline_replayed=checks,all_families_retained=True,driver_sha256=sha(__file__),fixed_roster=85,scope='Exploratory biochemical associations; no reaction-flux or causal inference; metabolite-identity qualifications retained.')
