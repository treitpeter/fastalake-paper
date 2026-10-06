"""Recount frozen donor and fusion controls; preserve their different estimands."""
import collections,csv,hashlib,json,math,os,re,sys,time
from pathlib import Path
import numpy as np
import pandas as pd
ROOT=Path(sys.argv[1]).resolve()
def sha(p):
 h=hashlib.sha256()
 with Path(p).open('rb') as f:
  for b in iter(lambda:f.read(8<<20),b''):h.update(b)
 return h.hexdigest()
def dump(p,v):Path(p).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def table(p,rows):
 rows=list(rows)
 with Path(p).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
def verify(p,d):assert sha(p)==d,(str(p),'hash changed')
def donor(i):
 m=json.loads((ROOT/'inputs/manifest.json').read_text());task=m['donor'][i];base=Path(task['output_directory']);out=ROOT/'results/donor'/str(i);out.mkdir(parents=True,exist_ok=False)
 for name,key in [('RESULT.json','result_record_sha256'),('canonical_scores.tsv','canonical_scores_sha256'),('search/results.sage.tsv','result_sha256'),('pairs.tsv','pairs_sha256'),('search.fasta','search_fasta_sha256'),('config.json','search_config_sha256')]:verify(base/name,task[key])
 old=json.loads((base/'RESULT.json').read_text());assert old['status']=='COMPLETE'
 for k in ['target','donor','seed','tier','group']:assert str(old[k])==str(task[k]),k
 import donor_original as original
 seqs,pairs=original.read_pairs(base/'pairs.tsv')
 best={};nrows=0
 with (base/'search/results.sage.tsv').open() as f:
  for row in csv.DictReader(f,delimiter='\t'):
   nrows+=1;p=re.sub(r'\[[^\]]*\]','',row['peptide']).replace('I','L');kind,pair=seqs[p];label=int(row['label']);assert label in [-1,1] and (label==-1)==(kind in ['decoy','p_decoy'])
   score=(float(row['peptide_q']),float(row['posterior_error']));assert all(map(math.isfinite,score)) and 0<=score[0]<=1
   if label==1 and (p not in best or score<best[p]):best[p]=score
 with (base/'canonical_scores.tsv').open() as f:
  saved={r['peptide']:(float(r['q_value']),float(r['posterior_error'])) for r in csv.DictReader(f,delimiter='\t')}
 assert saved==best,'raw-score reconstruction differs'
 rows=[]
 for low in [5,9]:
  for q in [.001,.005,.01,.02,.05]:
   acc={p for p,v in best.items() if low<=len(p)<=50 and v[0]<=q};T=E=E0=B=J=0
   for pair in pairs.values():
    t,e=pair['target'],pair['p_target'];ta,ea=t in acc,e in acc;T+=ta;E+=ea
    if ea:
     E0+=not ta
     if ta:B+=best[e]<best[t];J+=best[e]==best[t]
   N=T+E;assert N==len(acc);A=E+E0+2*B
   r={k:task[k] for k in ['index','target','donor','tier','seed','group']};r.update(min_length=low,q_threshold=q,n_target=T,n_entrapment=E,partner_undiscovered=E0,entrapment_ranks_higher=B,exact_score_ties=J,paired_numerator=A,lower_fraction=E/N if N else None,combined_estimate=2*E/N if N else None,paired_estimate=A/N if N else None,paired_ties_upper=(A+2*J)/N if N else None)
   ref=next(v for v in old['summaries'] if v['min_length']==low and v['q_threshold']==q)
   for k in ref:assert r[k]==ref[k],(k,r[k],ref[k])
   rows.append(r)
 table(out/'recount.tsv',rows)
 dump(out/'COMPLETE.json',dict(status='PASS',task=task,raw_rows=nrows,quartets=len(pairs),canonical_scores=len(best),original_summaries_match=True,outputs={'recount.tsv':sha(out/'recount.tsv')},job_id=os.getenv('SLURM_JOB_ID')))
 print('PASS donor',i,flush=True)
def fusion(i):
 import audit_fusion as audit
 m=json.loads((ROOT/'inputs/manifest.json').read_text());task=m['fusion'][i]
 verify(Path(task['directory'])/'COMPLETE.json',task['completion_sha256'])
 verify(ROOT/'frozen/audit_fusion.py',m['frozen_sha256']['audit_fusion.py'])
 audit.source_commit=lambda:'frozen-sha256:'+m['frozen_sha256']['audit_fusion.py']
 audit.run_task(task,ROOT/'results/fusion'/str(i))
def collect():
 import matplotlib;matplotlib.use('Agg')
 import matplotlib.pyplot as plt
 m=json.loads((ROOT/'inputs/manifest.json').read_text());out=ROOT/'results';don=[];fus=[];proof=[]
 for arm,tasks,name,rows in [('donor',m['donor'],'recount.tsv',don),('fusion',m['fusion'],'decomposition.tsv',fus)]:
  for i,t in enumerate(tasks):
   p=out/arm/str(i);c=json.loads((p/'COMPLETE.json').read_text())
   assert c.get('status')=='PASS' if arm=='donor' else c['independent_accounting_matches'] and c['independent_final_pairing_matches']
   verify(p/name,c['outputs'][name]);x=pd.read_csv(p/name,sep='\t');rows.append(x);proof.append(dict(path=str(p/'COMPLETE.json'),sha256=sha(p/'COMPLETE.json')))
 d=pd.concat(don,ignore_index=True);f=pd.concat(fus,ignore_index=True);assert len(d)==1890 and len(f)==21300
 source=pd.read_csv(ROOT/'inputs/donor_searches.tsv',sep='\t');keys=['index','min_length','q_threshold']
 merged=d.merge(source,on=keys,suffixes=('_new','_old'),validate='one_to_one');assert len(merged)==len(d)
 for field in ['n_target','n_entrapment','partner_undiscovered','entrapment_ranks_higher','exact_score_ties','paired_numerator','paired_estimate']:assert np.allclose(merged[field+'_new'],merged[field+'_old'],equal_nan=True,atol=1e-12,rtol=1e-12)
 d.to_csv(out/'donor_recount.tsv',sep='\t',index=False);f.to_csv(out/'fusion_recount.tsv',sep='\t',index=False)
 gkeys=['target','donor','tier','group','min_length','q_threshold'];g=d.groupby(gkeys,dropna=False);assert (g.size()==3).all();ds=g[['n_target','n_entrapment','paired_estimate','lower_fraction','combined_estimate']].mean().reset_index();ds.to_csv(out/'donor_seed_means.tsv',sep='\t',index=False)
 summary=[]
 for k,x in ds.groupby(['tier','min_length','q_threshold']):
  z=x.paired_estimate;summary.append(dict(tier=k[0],min_length=k[1],q_threshold=k[2],acquisitions=len(x),median=z.median(),q25=z.quantile(.25),q75=z.quantile(.75),minimum=z.min(),maximum=z.max()))
 sd=pd.DataFrame(summary);sd.to_csv(out/'donor_summary.tsv',sep='\t',index=False)
 fk=['cohort','seed','q_threshold','min_residues','max_residues'];ff=f.groupby(fk)[['n_canonical','paired_numerator','n_unpaired_accepted','n_entrapment_q_tied']].sum().reset_index();ff['diagnostic']=(ff.paired_numerator+ff.n_unpaired_accepted)/ff.n_canonical.replace(0,np.nan);ff['paired_component']=ff.paired_numerator/ff.n_canonical.replace(0,np.nan);ff['unpaired_component']=ff.n_unpaired_accepted/ff.n_canonical.replace(0,np.nan);ff.to_csv(out/'fusion_summary.tsv',sep='\t',index=False)
 # Compute a common descriptive aggregation as a sensitivity, without pretending the assays coincide.
 fs=f.query('max_residues==50').groupby(['cohort','sample','q_threshold','min_residues']).paired_diagnostic.mean().reset_index();fs.to_csv(out/'fusion_acquisition_seed_means.tsv',sep='\t',index=False)
 matched=ds.query('tier=="T0_self"').merge(ds.query('tier=="T1_within_group"'),on=['target','group','min_length','q_threshold'],suffixes=('_own','_transfer'),validate='one_to_one');matched['paired_delta']=matched.paired_estimate_transfer-matched.paired_estimate_own;assert len(matched)==50;matched.to_csv(out/'five_exact_own_transfer_comparisons.tsv',sep='\t',index=False)
 fig,axes=plt.subplots(1,2,figsize=(10,4));colors={5:'#b05d37',9:'#237a89'}
 for low,x in sd.query('tier=="T1_within_group"').groupby('min_length'):
  axes[0].plot(x.q_threshold*100,x['median']*100,'o-',color=colors[low],label=f'{low}–50 residues');axes[0].fill_between(x.q_threshold.to_numpy()*100,x.q25.to_numpy()*100,x.q75.to_numpy()*100,color=colors[low],alpha=.15)
 axes[0].plot([0,5],[0,5],':',color='grey');axes[0].set(xlabel='Nominal peptide q threshold (%)',ylabel='Donor paired diagnostic (%)',title='58 transfers: seed means, median and IQR');axes[0].legend(frameon=False)
 nom=ff.query('q_threshold==.01 and min_residues==9 and max_residues==50');labels=[]
 for i,(_,x) in enumerate(nom.iterrows()):
  axes[1].bar(i,100*x.paired_component,color='#237a89');axes[1].bar(i,100*x.unpaired_component,bottom=100*x.paired_component,color='#b05d37');labels.append(x.cohort+'\n'+str(int(x.seed))[-2:])
 axes[1].axhline(1,ls=':',color='grey');axes[1].set_xticks(range(len(labels)),labels,rotation=45,ha='right');axes[1].set(ylabel='Fusion diagnostic (%)',title='9–50 residues: paired + unresolved allowance');fig.tight_layout()
 for ext in ['png','pdf','svg']:fig.savefig(out/('control_reconciliation.'+ext),dpi=180)
 lines=['# Control reconciliation — fresh accounting','', 'All 189 donor searches and 1,420 fusion searches were recounted from their original PSM outputs. Donor canonical scores and all saved estimator components reproduce; fusion final protein-pair digestion and peptide counts reproduce. This validates accounting, not whole-procedure FDR calibration.','', '| Assay | Control entry | Score ordering | Diagnostic | Aggregation |','|---|---|---|---|---|','| Frozen donor | Peptide controls after donor protein selection | Ordered (peptide q, posterior error), strict comparison; ties separate | (E + E₀ + 2B)/(T+E) | Three ratios averaged per acquisition, then median/IQR |','| Linked fusion | Linked target/control proteins before evidence selection | Minimum peptide q; equal-q comparisons included conservatively | (paired numerator + unresolved accepted sequences)/accepted canonical sequences | Sum of counts within cohort and seed |','| Full-route unlinked biological/shuffle controls | Before selection; jointly competing control classes | Original production search | Observed class-label fractions | Separate descriptive evidence; not paired FDP |','', 'The same reporting lengths and an acquisition-weighted fusion sensitivity are supplied, but shared plotting axes cannot make these distinct experiments equivalent. No control hits overlapping the biological reservoir are subtracted. Five own-database examples are compared only with the same five transferred recipients. The 9–50 view is nested in the original 5–50 search; confidence is not recomputed.','', '## Recomputed donor result at q ≤ 0.01','']
 for _,r in sd.query('tier=="T1_within_group" and q_threshold==.01').iterrows():lines.append(f'- {int(r.min_length)}–50: median {100*r["median"]:.3f}%; IQR {100*r.q25:.3f}–{100*r.q75:.3f}%; {int(r.acquisitions)} acquisition means.')
 lines+=['','CAMPI fusion remains above 1% in the archived 9–50 diagnostic; the paired contribution and unresolved allowance remain visible separately. This must be discussed alongside the donor result rather than replaced by it. The completed known-truth simulation and unlinked controls are retained as separately labelled source tables in inputs; their truth/label definitions are not transferred to ordinary production identifications.','', 'Source table arithmetic, estimator definitions, per-search hashes, full length/threshold grids, and exact own/transfer comparisons are retained here. Supplement reconciliation and visual review remain required before manuscript text is finalized.']
 (out/'REPORT.md').write_text('\n'.join(lines)+'\n');dump(out/'AGGREGATE_VALIDATION.json',dict(status='PASS',donor_searches=189,fusion_searches=1420,donor_archived_table_matches=True,sources=proof))
 print('PASS combined accounting report')
if __name__=='__main__':
 mode=sys.argv[2]
 if mode=='donor':donor(int(sys.argv[3]))
 elif mode=='fusion':fusion(int(sys.argv[3]))
 elif mode=='collect':collect()
 else:raise ValueError(mode)
