from pathlib import Path
from datetime import datetime,timezone
import json,hashlib,subprocess,shutil
import pandas as pd,numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
H=Path(__file__).resolve().parent;J=H.parents[1];ROOT=J.parents[2];base=ROOT/'MM_review_analysis/A04_4_2_expansion/alphanovo_bags_2026-09-23'
T=ROOT.parent/'thesis/reviews/paragraph_walkthrough_2026-09-29/resolution_sources';source=T/'ALPHANOVO_FASTALAKE_BAGS_CLUSTER_REVIEW_2026-09-30.json';invdir=T/'bags_slurm_review_2026-09-30/validated_inventory'
load=lambda p:json.loads(p.read_text());sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,d):p.write_text(json.dumps(d,ensure_ascii=False,indent=2)+'\n')
d=load(source);hashes=[]
for ref in d['sources']+d['source_receipts']:
 p=Path(ref['path']);actual=sha(p);same=actual==ref['sha256'];mutable=p==base/'monitor_2026-09-26/STATUS.md';assert same or mutable,p
 hashes.append(dict(path=str(p),sha256=actual,prior_sha256=ref['sha256'],match=same,mutable_monitor=mutable))
receipts=sorted((base/'bounded_sweep_2026-09-25/outputs').glob('*/*/SEARCHES_COMPLETE.json'));assert len(receipts)==48;rows=[]
for p in receipts:
 r=load(p);assert r['status']=='MATCHED_SEARCHES_COMPLETE' and r['comparison_reference_arm']=='exact';rows+=r['rows']
f=pd.DataFrame(rows);assert len(f)==336 and not f.duplicated(['cohort','sample','arm']).any();assert (f.gained-f.lost==f.net).all();summary=f.groupby('arm')[['gained','lost','net']].sum()
for arm,expect in [('bags_all',(6961,5941,1020)),('bags_le8',(5051,3318,1733)),('mass_all',(31111,33306,-2195))]:assert tuple(summary.loc[arm])==expect
f.to_csv(H/'bounded48_per_acquisition.tsv',sep='\t',index=False);summary.to_csv(H/'bounded48_summary.tsv',sep='\t')
fa=base/'gain_loss_2026-09-24/fragment_validation/audit';v=load(fa/'VALIDATION.json');assert len(v['cases'])==6
for ref in v['sources']:assert sha(ref['path'])==ref['sha256'];hashes.append(dict(path=ref['path'],sha256=ref['sha256'],match=True))
fr=pd.read_csv(fa/'LOCAL_DIAGNOSTIC_SUMMARY.tsv',sep='\t');fr=fr[(fr.contrast=='bags_vs_exact')&(fr.change=='gained')&(fr.relation=='composition_preserving_rearrangement')];assert fr.spectrum_events.sum()==215 and fr.zero_diagnostic.sum()==193;fr.to_csv(H/'fragment_order_subset.tsv',sep='\t',index=False)
inv=load(invdir/'RECEIPT_INVENTORY.json')
control=load(Path(inv['entrapment_complete_union']['source']['path']));assert len(control['reports'])==58 and len({r['sample'] for r in control['reports']})==58
for r in control['reports']:assert sha(r['path'])==r['sha256'];hashes.append(dict(path=r['path'],sha256=r['sha256'],match=True))
for key in ['native_COMPLETE','entrapment_COMPLETE']:
 r=inv[key];assert sha(r['source']['path'])==r['source']['sha256'];assert r['receipt']['status']=='PASS' and r['receipt']['scope']=='ONE_ACQUISITION_FULL_RESERVOIR_PILOT'
patterns=[('Full-reservoir order/flip',base.parent/'full_reservoir_2026-09-23/outputs','*/*/searches/COMPLETE.json',710,'MATCHED_SEARCHES_COMPLETE','Five malformed-bag failures; incomplete cohort'),('Native bag/mass',base/'attempt_02/outputs','*/*/SEARCHES_COMPLETE.json',710,'MATCHED_SEARCHES_COMPLETE','Running; baseline mismatches remain'),('Matched shuffled route controls',base/'paired_nulls_2026-09-23/outputs','*/*/NULL_COMPARISON.json',710,'NULL_DIAGNOSTIC_COMPLETE_REVIEW_REQUIRED','Dependent on native inputs; diagnostic review required'),('Substitution-bound searches',base/'sub_bounds_2026-09-25/outputs','*/*/SEARCHES_COMPLETE.json',53,'MATCHED_SEARCHES_COMPLETE','Two pattern-limit failures; eleven CAMPI candidates held')]
status=[];allreceipts={}
for label,folder,pattern,expected,required,note in patterns:
 files=set(p.resolve() for p in folder.glob(pattern))
 key={'Full-reservoir order/flip':'order_searches','Native bag/mass':'native_searches','Matched shuffled route controls':'null_reports','Substitution-bound searches':'subbounds_searches'}[label]
 for ref in inv[key]['receipts']:
  pp=Path(ref['path']).resolve();assert sha(pp)==ref['sha256'];files.add(pp)
 files=sorted(files);states=[load(p).get('status') for p in files];assert all(s==required for s in states),(label,set(states));allreceipts[label]=[dict(path=str(p),sha256=sha(p)) for p in files];status.append(dict(experiment=label,completed=len(files),planned=expected,interpretation=note))
status += [dict(experiment='Bounded-bag comparison',completed=48,planned=48,interpretation='Completed exploratory matched-search comparison; nominal q only'),dict(experiment='Paired fragment audit',completed=6,planned=6,interpretation='Completed pairwise fragment audit; comparators are not known truth'),dict(experiment='Original-parser entrapment controls',completed=58,planned=58,interpretation='Three replacements included; sparse exposure and parser limitations'),dict(experiment='Corrected-parser full-reservoir pilot',completed=1,planned=1,interpretation='Native and augmented controls passed for one acquisition; not cohort validation')]
pd.DataFrame(status).to_csv(H/'experiment_status.tsv',sep='\t',index=False)
# Freeze the inspected scheduler state without mutating or resubmitting jobs.
jobids=list(load(invdir/'JOB_SUMMARY.json'));assert len(jobids)==106
for cmd,name in [(['squeue','-u','treit','-h','-o','%i|%j|%T|%R'],'squeue_snapshot.txt'),(['sacct','-j',','.join(jobids),'-X','-n','-P','--format=JobID,JobName,State,ExitCode'],'sacct_snapshot.tsv')]:
 r=subprocess.run(cmd,capture_output=True,text=True,timeout=55);assert r.returncode==0;rpath=H/name;rpath.write_text(r.stdout)
# Keep source provenance and the original staged-work inventory as separate records.
for p in [source,invdir/'RECEIPT_INVENTORY.json',invdir/'SUBMISSION_INVENTORY.json',invdir/'COMPLETE.json',T/'bags_slurm_review_2026-09-30/REVIEW_FOR_PETER.txt']:shutil.copy2(p,H/('source_'+p.name))
fig,axes=plt.subplots(2,1,figsize=(9,6.8),gridspec_kw={'height_ratios':[4,1.3]});ax=axes[0]
order=['bags_le2','bags_le4','bags_le8','bags_all','mass_all','comb_all'];labels=['Bags: ≤2 orderings','Bags: ≤4 orderings','Bags: ≤8 orderings','Bags: unrestricted','Mass substitutions','Bags + mass'];ss=summary.loc[order];yy=np.arange(len(order));ax.barh(yy,ss.gained,color='#0072B2',label='Gained');ax.barh(yy,-ss.lost,color='#D55E00',label='Lost');ax.axvline(0,color='.2',lw=.8);ax.set_yticks(yy,labels);ax.invert_yaxis();ax.set_xlim(-42000,48000)
for i,(_,r) in enumerate(ss.iterrows()):ax.text(-r.lost-600,i,f'{int(r.lost):,}',ha='right',va='center',fontsize=8);ax.text(r.gained+600,i,f'{int(r.gained):,}',va='center',fontsize=8);ax.text(47500,i,f'{int(r.net):+,}',ha='right',va='center',fontsize=8,fontweight='bold')
ax.set_xticks([-30000,0,30000],['30,000 lost','0','30,000 gained']);ax.set_xlabel('Peptide–acquisition identifications relative to exact lookup');ax.set_title('a  Completed matched comparison: 48 acquisitions',loc='left');ax.legend(loc='upper left',frameon=False,fontsize=8);ax.text(.96,1.035,'Net',transform=ax.transAxes,ha='right',fontsize=9)
bx=axes[1];bx.barh([0],[193],color='#999999',label='No diagnostic ion');bx.barh([0],[22],left=[193],color='#009E73',label='≥1 diagnostic ion');bx.text(96.5,0,'193',ha='center',va='center');bx.text(204,0,'22',ha='center',va='center');bx.set_xlim(0,215);bx.set_yticks([]);bx.set_xticks([0,50,100,150,200,215]);bx.set_xlabel('Order-related same-spectrum comparison pairs');bx.set_title('b  Separate fragment audit: six acquisitions',loc='left');bx.legend(loc='upper center',bbox_to_anchor=(.5,-.55),ncol=2,frameon=False,fontsize=8)
fig.tight_layout(h_pad=2.0,rect=[0,.07,1,1]);fig.savefig(H/'Supplementary_Figure16_bags.pdf');fig.savefig(H/'Supplementary_Figure16_bags.png',dpi=200);plt.close(fig)
review=dict(status='PASS_RECEIPT_RECOUNT_AND_SOURCE_HASH_REVIEW',utc=datetime.now(timezone.utc).isoformat(),source_hashes=hashes,receipt_inventory=allreceipts,job_ids=jobids,experiment_status=status,bounded48_summary=summary.reset_index().to_dict('records'),fragment_pairs=215,fragment_zero_diagnostic=193,scope='Independent receipt/table recount and source hash verification, not a fresh 336-search raw-result reanalysis. Mutable monitor changes recorded; scientific frozen sources unchanged. Pending cohorts not used for completed aggregate claims.',visual_review_pending=True)
save(H/'REVIEW.json',review);print(pd.DataFrame(status).to_string(index=False));print('Static hashes and receipts verified; mutable monitor preserved as changing status only.')
