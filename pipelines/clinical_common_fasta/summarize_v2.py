from pathlib import Path
import json,hashlib,csv,sys,datetime,math
import numpy as np,pandas as pd
B=Path(__file__).resolve().parents[1];E=B.parent;co=sys.argv[1];outname=sys.argv[2] if len(sys.argv)>2 else co
# Derivative of summarize.py (sha256 51a67c7f270dc8e52b291d9cdcc88bce7846aab1b62dff4324a73fad6437380e): acquisitions without a participant mapping
# (canonical metadata 'UNMATCHED-no-patient-mapping', participant '') stay in every acquisition-level table and
# check, are exported separately, and are excluded from specimen/participant averages, effects and bootstrap.
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def meta(p):return dict(path=str(Path(p).resolve()),sha256=sha(p))
def load(p):return json.loads(Path(p).read_text())
def write(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
paired=E/'common_cache/comparison/results/full'/co;native=E/'original_stitch/comparison/full'/co
pv=load(paired/'COMPLETE.json');nv=load(native/'COMPLETE.json');assert pv['status']==nv['status']=='PASS'
for rec in pv['outputs']:assert sha(rec['path'])==rec['sha256']
for rec in nv['sources']:assert sha(rec['path'])==rec['sha256']
f=load(E/'clinical_expansion/inputs/FROZEN.json');manifest=Path(f['run'])/'manifest.json';expected=next(x for x in f['files'] if Path(x['path'])==manifest);assert sha(manifest)==expected['sha256'];tasks=[x for x in load(manifest)['tasks'] if x['cohort']==co]
md=pd.DataFrame([dict(sample=t['sample'],participant=str(t['design']['participant']),specimen=str(t['design']['specimen_id']),region=str(t['design']['region'])) for t in tasks]);assert md['sample'].is_unique and not md.isna().any().any()
d=pd.read_csv(paired/'paired_acquisition_comparison.tsv',sep='\t');assert len(d)==len(md)=={'PREDICT':306,'CAPSCAN':331}[co];assert set(d['sample'])==set(md['sample']);assert not d['sample'].duplicated().any()
for prefix,suffix in [('','peptides'),('positive_','features')]:
 old=d['original_'+prefix+suffix];new=d['common_'+prefix+suffix];shared=d['shared_'+prefix+suffix]
 np.testing.assert_array_equal(old,shared+d['lost_'+prefix+suffix]);np.testing.assert_array_equal(new,shared+d['gained_'+prefix+suffix])
assert d.baseline_count_difference_from_common_receipt.eq(0).all()
n=pd.read_csv(native/'native_per_acquisition.tsv',sep='\t');assert not n.duplicated(['sample','arm','method']).any();assert len(n)==4*len(d)
w=n.pivot(index='sample',columns=['arm','method'],values='positive_groups');w.columns=[a+'_'+m+'_groups' for a,m in w.columns];d=d.merge(w,left_on='sample',right_index=True,validate='one_to_one').merge(md,on='sample',validate='one_to_one')
cv=pd.read_csv(native/'native_coverage.tsv',sep='\t');
for r in cv.itertuples(index=False):
 assert int(d[r.arm+'_'+r.method+'_groups'].sum())==r.positive_cells
 assert math.isclose(1-r.positive_cells/(r.native_rows*len(d)),r.native_missing_fraction,abs_tol=1e-12)
for arm in ['original','common']:
 d[arm+'_directlfq_group_cells_lost']=d[arm+'_sum_groups']-d[arm+'_directlfq_groups'];assert d[arm+'_directlfq_group_cells_lost'].ge(0).all()
cols=['original_peptides','common_peptides','gained_peptides','lost_peptides','original_positive_features','common_positive_features','gained_positive_features','lost_positive_features','original_sum_groups','common_sum_groups','original_directlfq_groups','common_directlfq_groups','original_directlfq_group_cells_lost','common_directlfq_group_cells_lost']
un=d[d.participant.eq('')].copy();dp=d[d.participant.ne('')].copy();assert len(un)+len(dp)==len(d)
spec=dp.groupby(['participant','region','specimen'],sort=True)[cols].mean().reset_index();don=spec.groupby('participant',sort=True)[cols].mean().reset_index();reg=spec.groupby(['participant','region'],sort=True)[cols].mean().reset_index()
metrics=['peptides','positive_features','sum_groups','directlfq_groups'];summ=[]
for level,x in [('participant',don),('participant_region',reg)]:
 for metric in metrics:
  old=x['original_'+metric];new=x['common_'+metric];x[metric+'_difference']=new-old;x[metric+'_percent_difference']=100*(new/old.where(old>0)-1)
 if level=='participant':groups=[('all_observed_specimens',x)]
 else:groups=list(x.groupby('region',sort=True))
 for region,g in groups:
  for metric in metrics:
   for kind in ['difference','percent_difference']:
    values=g[metric+'_'+kind].dropna().to_numpy();nn=len(values);assert np.isfinite(values).all();rng=np.random.default_rng(20260915);weights=rng.multinomial(nn,np.full(nn,1/nn),size=20000) if nn else None
    ci=np.quantile(weights@values/nn,[.025,.975]) if nn>=8 else [None,None]
    summ.append(dict(cohort=co,level=level,region=region,metric=metric,effect=kind,participants=len(g),defined_participants=nn,mean=float(values.mean()) if nn else None,median=float(np.median(values)) if nn else None,minimum=float(values.min()) if nn else None,maximum=float(values.max()) if nn else None,positive=int((values>0).sum()),negative=int((values<0).sum()),zero=int((values==0).sum()),CI95_low=ci[0],CI95_high=ci[1]))
# Independent ordinary-Python reconstruction of hierarchical participant means.
for pid,g in dp.groupby('participant'):
 for c in cols:
  byspec={}
  for r in g[['region','specimen',c]].itertuples(index=False,name=None):byspec.setdefault(r[:2],[]).append(r[2])
  expected=sum(sum(v)/len(v) for v in byspec.values())/len(byspec)
  assert math.isclose(float(don.loc[don.participant.eq(pid),c].iloc[0]),expected,rel_tol=1e-12,abs_tol=1e-9),(pid,c)
# Cross-check observed-identity denominators independently from frequency exports.
source_summary=load(paired/'SUMMARY.json')
for filename,key in [('peptide_observation_counts.tsv','peptide_coverage'),('feature_observation_counts.tsv','feature_coverage')]:
 counts={'original':dict(vocabulary=0,any=0,half=0,pct80=0,all=0,positive_cells=0),'common':dict(vocabulary=0,any=0,half=0,pct80=0,all=0,positive_cells=0)}
 with (paired/filename).open() as stream:
  for r in csv.DictReader(stream,delimiter='\t'):
   for arm in counts:
    v=int(r[arm]);assert 0<=v<=len(d);c=counts[arm];c['vocabulary']+=1;c['any']+=v>0;c['half']+=v>=math.ceil(len(d)/2);c['pct80']+=v>=math.ceil(.8*len(d));c['all']+=v==len(d);c['positive_cells']+=v
 for arm,c in counts.items():
  for k,v in c.items():assert source_summary[key][arm][k]==v,(key,arm,k)
  assert math.isclose(source_summary[key][arm]['missing_fraction'],1-c['positive_cells']/(c['vocabulary']*len(d)),abs_tol=1e-12)
out=B/outname;out.mkdir(exist_ok=False)
for name,x in [('per_acquisition',d),('per_specimen',spec),('per_participant',don),('per_participant_region',reg),('participant_summary',pd.DataFrame(summ)),('native_coverage',cv)]:x.to_csv(out/(name+'.tsv'),sep='\t',index=False)
zero=d.loc[d.original_peptides.eq(0)|d.original_sum_groups.eq(0)|d.original_directlfq_groups.eq(0)];zero.to_csv(out/'zero_original_denominators.tsv',sep='\t',index=False)
un.to_csv(out/'unassigned_participant_acquisitions.tsv',sep='\t',index=False)
acq=[]
for metric in metrics:
 old=d['original_'+metric];new=d['common_'+metric];rat=new/old.where(old>0)
 acq.append(dict(metric=metric,original_median=float(old.median()),common_median=float(new.median()),paired_ratio_median=float(rat.median()),ratio_n=int(rat.notna().sum()),higher=int((new>old).sum()),lower=int((new<old).sum()),equal=int((new==old).sum())))
agr=pd.read_csv(native/'exact_assignment_agreement.tsv',sep='\t');med=agr.groupby('method').median(numeric_only=True).reset_index();med.to_csv(out/'exact_assignment_medians.tsv',sep='\t',index=False)
write(out/'SUMMARY.json',dict(cohort=co,acquisitions=len(d),specimens=len(spec),participants=len(don),unassigned_participant_acquisitions=len(un),unassigned_samples=sorted(un['sample']),participant_level_acquisitions=len(dp),acquisition_summaries=acq,zero_denominator_acquisitions=len(zero),paired_summary=source_summary,native_summary=nv,limits='Descriptive observed-specimen mixture, independently inferred native dictionaries; no joint Sage-search, biological truth, clinical technical precision, or calibrated FDR claim.'))
write(out/'VALIDATION.json',dict(status='PASS',cohort=co,acquisitions=len(d),participants=len(don),unassigned_participant_acquisitions=len(un),participant_level_acquisitions=len(dp),derived_from_summarize_sha256='51a67c7f270dc8e52b291d9cdcc88bce7846aab1b62dff4324a73fad6437380e',independent_hierarchical_mean_checks=len(don)*len(cols),independent_frequency_denominators=True,source_records=[meta(paired/'COMPLETE.json'),meta(native/'COMPLETE.json'),meta(manifest)],code=meta(__file__),outputs=[meta(p) for p in sorted(out.iterdir())],scientific_review='PENDING',checked_utc=datetime.datetime.now(datetime.timezone.utc).isoformat()))
print('CLINICAL QUANT SYNTHESIS PASS',co,len(d),len(don),'unassigned',len(un),flush=True)
