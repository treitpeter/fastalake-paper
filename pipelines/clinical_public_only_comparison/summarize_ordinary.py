"""clinical_public_only_comparison: summarize ordinary. Moved from the original analysis as it ran, apart from path scrubbing; see README.md."""
from pathlib import Path
import json,hashlib,datetime
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
HERE=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
m=json.load(open(HERE/'ordinary_manifest.json'));frames=[];records=[];gates={};inputs=[]
for i,chunk in enumerate(m['chunks']):
 folder=HERE/'ordinary_outputs'/str(i);r=json.load(open(folder/'COMPLETE.json'));assert r['status']=='PASS' and r['n_acquisitions']==len(chunk) and r['manifest_sha256']==sha(HERE/'ordinary_manifest.json')
 assert sha(folder/'counts.tsv')==r['counts_sha256'];records.extend(r['records']);frames.append(pd.read_csv(folder/'counts.tsv',sep='\t',dtype={'subject_id':str,'specimen_id':str}))
 for rec in r['records']:
  for arm,x in rec['inputs'].items():
   for k,v in x['gates'].items():gates[k]=gates.get(k,0)+v
 inputs.append(dict(path=str(folder/'COMPLETE.json'),sha256=sha(folder/'COMPLETE.json')))
# Substitute only the five source-repaired searches already used by v15.
cor=HERE/'ordinary_current_five';cr=json.load(open(cor/'ordinary_outputs/0/COMPLETE.json'));assert cr['status']=='PASS' and cr['n_acquisitions']==5
assert cr['manifest_sha256']==sha(cor/'ordinary_manifest.json');assert cr['counts_sha256']==sha(cor/'ordinary_outputs/0/counts.tsv')
replace_ids={x['index'] for x in cr['records']};assert len(replace_ids)==5
frames=[g[~g['index'].isin(replace_ids)] for g in frames];frames.append(pd.read_csv(cor/'ordinary_outputs/0/counts.tsv',sep='\t',dtype={'subject_id':str,'specimen_id':str}))
records=[r for r in records if r['index'] not in replace_ids]+cr['records'];gates={}
for r in records:
 for arm,x in r['inputs'].items():
  for k,v in x['gates'].items():gates[k]=gates.get(k,0)+v
inputs.append(dict(path=str(cor/'ordinary_outputs/0/COMPLETE.json'),sha256=sha(cor/'ordinary_outputs/0/COMPLETE.json')))
assert sorted(r['index'] for r in records)==list(range(578))
f=pd.concat(frames,ignore_index=True);assert len(f)==5780 and not f.duplicated(['index','q','min_length']).any()
assert (f.shared+f.fl_only==f.n_fl).all() and (f.shared+f.mg_only==f.n_mg).all()
f['ratio']=f.n_fl/f.n_mg.replace(0,np.nan);f['gain_pct']=(f.ratio-1)*100
base=f[f.q.eq(.01)].set_index(['index','min_length']);keys=pd.MultiIndex.from_frame(f[['index','min_length']])
f['fl_retained_pct']=100*f.n_fl.to_numpy()/base.n_fl.reindex(keys).replace(0,np.nan).to_numpy();f['mg_retained_pct']=100*f.n_mg.to_numpy()/base.n_mg.reindex(keys).replace(0,np.nan).to_numpy()
f.to_csv(HERE/'source_data/ordinary_all_pairs.tsv',sep='\t',index=False)
ROOT=next(p for p in HERE.parents if (p/'00_START_HERE_V7.md').exists())
paper=ROOT/'00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/paired_benchmark_rows.tsv'
pap=pd.read_csv(paper,sep='\t');joined=f[f.q.eq(.01)&f.min_length.eq(5)].merge(pap,on=['cohort','sample'],validate='one_to_one')
assert len(joined)==578 and (joined.n_fl==joined.n_test).all() and (joined.n_mg==joined.n_reference).all() and (joined.shared==joined.overlap).all()
inputs.append(dict(path=str(paper),sha256=sha(paper)))
# Constant denominator-defined subset across all thresholds, to expose population changes.
common=f.groupby(['cohort','index']).n_mg.min();good={idx for (co,idx),v in common.items() if v>0};fixed=f[f['index'].isin(good)];fixed.to_csv(HERE/'source_data/ordinary_fixed_defined_pairs.tsv',sep='\t',index=False)
fixed.groupby(['cohort','q','min_length']).agg(pairs=('index','size'),median_gain_pct=('gain_pct','median'),median_fl_retained_pct=('fl_retained_pct','median'),median_mg_retained_pct=('mg_retained_pct','median')).to_csv(HERE/'source_data/ordinary_fixed_defined_summary.tsv',sep='\t')

rows=[]
for (co,q,lo),g in f.groupby(['cohort','q','min_length']):
 rows.append(dict(cohort=co,level='acquisition',stratum='ALL',q=q,min_length=int(lo),n=len(g),defined_ratios=int(g.ratio.notna().sum()),zero_mg=int(g.n_mg.eq(0).sum()),zero_fl=int(g.n_fl.eq(0).sum()),median_gain_pct=g.gain_pct.median(),median_fl_retained_pct=g.fl_retained_pct.median(),median_mg_retained_pct=g.mg_retained_pct.median(),median_fl_only=g.fl_only.median(),median_mg_only=g.mg_only.median(),fl_more=int(g.n_fl.gt(g.n_mg).sum()),equal=int(g.n_fl.eq(g.n_mg).sum()),fl_fewer=int(g.n_fl.lt(g.n_mg).sum())))
metrics=['gain_pct','fl_retained_pct','mg_retained_pct','fl_only','mg_only'];eligible=f[f.subject_id.notna()&f.subject_id.ne('')];excluded=f[f.subject_id.isna()|f.subject_id.eq('')];excluded[['cohort','sample','index']].drop_duplicates().to_csv(HERE/'source_data/ordinary_participant_exclusions.tsv',sep='\t',index=False)
spec=eligible.groupby(['cohort','q','min_length','subject_id','specimen_id','group'])[metrics].mean().reset_index();part=spec.groupby(['cohort','q','min_length','subject_id'])[metrics].mean().reset_index();part.to_csv(HERE/'source_data/ordinary_participant_sensitivity.tsv',sep='\t',index=False)
reg=spec.groupby(['cohort','q','min_length','group','subject_id'])[metrics].mean().reset_index();reg.to_csv(HERE/'source_data/ordinary_participant_strata.tsv',sep='\t',index=False)
for frame,stratify in [(part,False),(reg,True)]:
 keys=['cohort','q','min_length']+(['group'] if stratify else [])
 for key,g in frame.groupby(keys):
  co,q,lo=key[:3];r=dict(cohort=co,level='participant',stratum=key[3] if stratify else 'ALL',q=q,min_length=int(lo),n=len(g),defined_ratios=int(g.gain_pct.notna().sum()))
  r.update({'median_'+col:g[col].median() for col in metrics});rows.append(r)
s=pd.DataFrame(rows);s.to_csv(HERE/'source_data/ordinary_summary.tsv',sep='\t',index=False)
# Plot all prespecified q values, without selecting the most favourable point.
plt.rcParams.update({'font.size':8,'axes.labelsize':8,'axes.titlesize':8,'xtick.labelsize':7,'ytick.labelsize':7,'legend.fontsize':7,'pdf.fonttype':42,'svg.fonttype':'none'})
fig,axs=plt.subplots(2,2,figsize=(7.2,5.0),sharex=True)
qs=[.001,.005,.01,.02,.05];colors=['#297ea0','#b95332']
for j,co in enumerate(['PREDICT','CAPSCAN']):
 for lo,c in zip([5,9],colors):
  g=s[(s.cohort.eq(co))&s.level.eq('acquisition')&s.min_length.eq(lo)].set_index('q').loc[qs]
  axs[0,j].plot(range(5),g.median_gain_pct,'o-',color=c,ms=3,lw=1,label=str(lo)+'–50 residues')
  axs[1,j].plot(range(5),g.median_fl_retained_pct,'o-',color=c,ms=3,lw=1,label='FastaLake, '+str(lo)+'–50')
  axs[1,j].plot(range(5),g.median_mg_retained_pct,'s--',color=c,ms=3,lw=1,label='Metagenome, '+str(lo)+'–50')
 axs[0,j].text(.04,.06,'q ≤ 0.001: '+(str(int(s[s.cohort.eq(co)&s.level.eq('acquisition')&s.q.eq(.001)&s.min_length.eq(5)].defined_ratios.iloc[0]))+'/'+str(296 if co=='PREDICT' else 282))+' defined ratios',transform=axs[0,j].transAxes,fontsize=7);axs[0,j].set_title(co+' (%d matched acquisitions)'%(296 if co=='PREDICT' else 282));axs[0,j].axhline(0,c='.6',lw=.6);axs[1,j].axhline(100,c='.6',lw=.6);axs[1,j].set_xticks(range(5),['0.001','0.005','0.01','0.02','0.05']);axs[1,j].set_xlabel('Reported peptide q threshold')
 for ax in axs[:,j]:ax.spines[['top','right']].set_visible(False)
axs[0,0].set_ylabel('Median paired identification\ngain (%)');axs[1,0].set_ylabel('Median yield retained (%)\nvs same-length q ≤ 0.01')
axs[0,0].legend(frameon=False);axs[1,0].legend(frameon=False,fontsize=7)
for k in [0,1]:
 limits=[ax.get_ylim() for ax in axs[k,:]];low=min(x[0] for x in limits);high=max(x[1] for x in limits)
 for ax in axs[k,:]:ax.set_ylim(low,high)
fig.tight_layout()
for ext in ['pdf','svg','png']:fig.savefig(HERE/('figures/Ordinary_threshold_sensitivity.'+ext),dpi=220)
record=dict(status='PASS',utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),pairs=578,cohorts=f.groupby('cohort')['index'].nunique().to_dict(),rows=len(f),chunk_receipts=inputs,all_q01_identity_sets_equal=True,all578_v15_main_count_and_overlap_rows_equal=True,five_repaired_sources_substituted=sorted(replace_ids),all_gain_loss_accounting_passed=True,gates=gates,zero_denominator_rows=int(f.n_mg.eq(0).sum()),workflow_fdr_established=False,scope='Reporting filter sensitivity of existing searches; no q re-estimation; matched FL_public versus MG',aggregation='Mean acquisition paired ratios within specimen then mean specimens within participant; median participant means; no acquisition-level inferential tests')
(HERE/'ORDINARY_VALIDATION.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps(record['gates']));print(s[s.stratum.eq('ALL')&s.q.isin([.001,.01])][['cohort','level','q','min_length','n','median_gain_pct','median_fl_retained_pct','median_mg_retained_pct']].round(3).to_string(index=False))
