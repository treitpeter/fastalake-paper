from common import *
from itertools import combinations
from scipy.stats import rankdata,spearmanr,pearsonr
from scipy.spatial.distance import squareform,pdist
from types import SimpleNamespace
import math
M=R/'00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23'
mode=sys.argv[1];assert mode in ['baseline','refit'];verify()
for p,h in load(H/'DISTANCE_INPUTS.json').items():assert sha(p)==h,p
assert sha(__file__)==load(H/'DISTANCE_CODE.json')['distance_sensitivity.py']
if mode=='refit':gate(H/'distance_baseline');gate(H/'cross_refit')
out=H/('distance_'+mode);out.mkdir(exist_ok=False)
env=dict(np=np,pd=pd,sparse=sparse,math=math,rankdata=rankdata,spearmanr=spearmanr,pearsonr=pearsonr,pdist=pdist,squareform=squareform,combinations=combinations)
functions(X/'native/native_pyramids.py',['bh'],env);env['native']=SimpleNamespace(SEED=20260907,bh=env['bh'])
functions(M/'experiment.py',['projection','relative','rank_gram','gram_distance','met_matrix','gram_for','rank_checks','correlate'],env);env['rank_checks']()
alltables=[];errors={};archives={}
for co in ['PREDICT','CAPSCAN']:
 old=np.load(M/'distances'/co/'distance_matrices.npz');oldtable=rd(M/'distances'/co/'comparisons.tsv');prepared=M/'prepared'/co
 if mode=='baseline':layers={k:sparse.load_npz(prepared/(k.replace(' ','_')+'.npz')) for k in ['Non-host protein','Host protein']}
 else:
  layers={};pred=rd(OLD/'PREDICT_fullstack_primary_samples.tsv',index_col=0)
  if co=='CAPSCAN':
   des=rd(OLD/'CAPSCAN_specimen_design.tsv',index_col=0);des.specimen=des.specimen.astype(str);dspec=rd(H/'cross_refit/capscan_specimen_roster.tsv',index_col=0);dspec.index=dspec.index.astype(str);ci=pd.MultiIndex.from_product([DONORS,REG],names=['subject_id','group'])
  for key,folderkey in [('Non-host protein','other_unresolved_groups'),('Host protein','supported_human_groups')]:
   raw=frame(H/'groups'/co,folderkey)
   if co=='PREDICT':
    ids=raw.index.to_series().str.extract(r'(MicrobPred_\d+)',expand=False);raw=raw.groupby(ids).mean().loc[pred.index]
   else:
    x=raw.loc[des.index];rel=x.div(x.sum(1).replace(0,np.nan),axis=0);spec=rel.groupby(des.specimen).median();si=spec.index.intersection(dspec.index);raw=spec.loc[si].groupby([dspec.loc[si,'subject_id'],dspec.loc[si,'group']]).mean().loc[ci]
   layers[key]=sparse.csr_matrix(raw.to_numpy())
 layers['Metabolites']=np.load(prepared/'Metabolites_log_profiles.npy');rows=[];quality={};archive={}
 for variant in ['all_observed','prevalence10','wuyts_recipe']:
  for metric in ['hellinger','rank']:
   if variant=='wuyts_recipe' and metric=='hellinger':continue
   grams={};quality[variant+'__'+metric]={}
   for key,x in layers.items():grams[key],quality[variant+'__'+metric][key]=env['gram_for'](x,key,variant,metric)
   for adj in ['unadjusted']+(['region_adjusted'] if co=='CAPSCAN' else []):
    prefix='__'.join([variant,metric,adj]);ds={k:old[prefix+'__'+k] for k in (['DNA','RNA'] if co=='PREDICT' else ['DNA'])}
    for key,g in grams.items():
     d=env['gram_distance'](g,adj=='region_adjusted',metric=='rank');ds[key]=d
     if mode=='baseline' or key=='Metabolites':
      err=float(abs(d-old[prefix+'__'+key]).max());assert err<1e-7,(co,key,variant,metric,adj,err);errors[co+'__'+prefix+'__'+key]=err
    archive.update({prefix+'__'+k:d for k,d in ds.items()})
    for outer in ['pearson','spearman']:
     t=env['correlate'](ds,co,outer);t.insert(1,'features',variant);t.insert(2,'distance',metric);t.insert(3,'adjustment',adj);t.insert(4,'outer_correlation',outer);rows.append(t)
 table=pd.concat(rows,ignore_index=True);keys=['cohort','features','distance','adjustment','outer_correlation','layer_a','layer_b'];comp=oldtable.merge(table,on=keys,validate='one_to_one',suffixes=('_previous','_revised'));assert len(comp)==len(table)==len(oldtable)
 if mode=='baseline':
  for c in ['coefficient','p','p_positive_tail','q']:assert np.allclose(comp[c+'_previous'],comp[c+'_revised'],atol=1e-10),(co,c)
 wr(comp,out/(co+'_method_comparison.tsv'));np.savez_compressed(out/(co+'_distances.npz'),**archive);dump(out/(co+'_quality.json'),quality);alltables.append(table);print(mode,co,len(table),'PASS',flush=True)
table=pd.concat(alltables,ignore_index=True);wr(table,out/'comparisons.tsv');controlled=table[table.features!='wuyts_recipe'];assert len(controlled)==176;wr(controlled,out/'controlled_176_comparisons.tsv')
receipt(out,mode=mode,all_comparisons=len(table),controlled_comparisons=len(controlled),cached_gene_distances='Unchanged sequence-layer source and exact axes verified independently',native_control_errors=errors,driver_sha256=sha(__file__),scope='Full two-sided native method grid and exploratory recipe sensitivity; no MOFA fit or causal ranking of methods')
