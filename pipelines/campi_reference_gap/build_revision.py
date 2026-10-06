#!/usr/bin/env python3
"""Portable, deterministic redraw of current-paper results; no search/threshold refit."""
from pathlib import Path
import json, shutil, hashlib, re, os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
P=Path(__file__).resolve().parent
D=Path(os.environ.get('FASTALAKE_SUPPLEMENT_OUTPUT',str(P/'delivery'))); DATA=D/'source_data'; FIG=D/'figures'
for x in [D,DATA,FIG,D/'code',D/'validation']: x.mkdir(exist_ok=True,parents=True)
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':8,'axes.titlesize':9,'axes.labelsize':8,'xtick.labelsize':7,'ytick.labelsize':7,'pdf.fonttype':42,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False,'savefig.dpi':180})
COL=['#2878a3','#df8f44','#619b68','#b35872','#8277ac']
PANELS=[]; FILES=[]; CHECKS=[]; CAPTIONS={}; INPUTS=set()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def src(pattern):
 a=list((P/'source_data').glob(pattern));assert len(a)==1,(pattern,a);INPUTS.add(a[0]);return a[0]
def read(pattern):return pd.read_csv(src(pattern),sep='\t')
def export(df,name,origin,meaning):
 f=DATA/(name+'.tsv');df.to_csv(f,sep='\t',index=False)
 FILES.append(dict(file=str(f.relative_to(D)),origin=str(origin),role=meaning,sha256=sha(f),rows=len(df)))
 return str(f.relative_to(D))
def copy(f,folder='source_data/underlying'):
 f=Path(f); to=D/folder/f.name;to.parent.mkdir(parents=True,exist_ok=True)
 if to.exists():assert sha(to)==sha(f),(to,f)
 else:shutil.copy2(f,to)
 FILES.append(dict(file=str(to.relative_to(D)),origin=str(f),role='unchanged producing input/output; see panel or table registry',sha256=sha(to),rows=''))
 return str(to.relative_to(D))
def panel(fig,letter,endpoint,file,columns,selection,unit,method,code,claim,validation='plotted quantities inspected; upstream inference not independently rerun'):
 PANELS.append(dict(figure=fig,panel=letter,endpoint=endpoint,source_file=file,columns=columns,row_selection=selection,unit=unit,methods=method,producing_code=code,main_text_anchor=claim,validation=validation,redraw_code='code/build_revision.py'))
def check(name,ok,details):
 CHECKS.append(dict(check=name,passed=bool(ok),details=str(details)));assert ok,(name,details)
def save(fig,name,title,caption):
 fig.suptitle(title,fontsize=11,fontweight='bold')
 fig.tight_layout(rect=(0,0,1,.96),pad=1.6)
 if name=='HISTIDINE_MEASURED_SPECTRA':
  fig.canvas.draw();renderer=fig.canvas.get_renderer()
  for axis in fig.axes:
   placed=[]
   labels=sorted([t for t in axis.texts if re.fullmatch(r'[by]\d+(?:²⁺)?',t.get_text())],key=lambda t:t.get_position())
   for label in labels:
    x,y=label.get_position();initial=y;box=label.get_window_extent(renderer).expanded(1.2,1.1)
    for attempt in range(12):
     if not any(box.overlaps(b) for b in placed):break
     y+=10;label.set_position((x,y));box=label.get_window_extent(renderer).expanded(1.2,1.1)
    assert not any(box.overlaps(b) for b in placed),label.get_text()
    assert y<120,(label.get_text(),y)
    if y!=initial:axis.plot([x,x],[initial,y-2],color='.5',lw=.35)
    placed.append(box)
 for ext in ['pdf','svg','png']:fig.savefig(FIG/(name+'.'+ext))
 plt.close(fig);CAPTIONS[name]=caption
def dots(ax,groups,ylabel,percent=False):
 for i,(label,values) in enumerate(groups):
  v=np.asarray(values,dtype=float);v=v[np.isfinite(v)]*(100 if percent else 1)
  x=i+np.linspace(-.16,.16,len(v));ax.scatter(x,v,s=8,alpha=.48,color=COL[i%5],rasterized=len(v)>100)
  if len(v):ax.hlines(np.median(v),i-.25,i+.25,color='black',lw=1.5)
 ax.set_xticks(range(len(groups)),[g[0] for g in groups]);ax.set_ylabel(ylabel)
def title(ax,l,t):ax.set_title(l+'  '+t,loc='left',pad=9)

# Selection: distinct query denominator; sweep does not collapse I/L.
q=read('original_S01/*length_by_arm.tsv');f=export(q,'selection_length_lookup',src('original_S01/*length_by_arm.tsv'),'Exact counts and all five seeds')
check('Lookup fractions',np.allclose(q.fraction_matched,q.matching_queries/q.distinct_queries,equal_nan=True),'matched distinct queries / distinct queries')
w=src('original_S23/*.xlsx');a=pd.read_excel(w,sheet_name='Arm_summary');obs=pd.read_excel(w,sheet_name='Per_acquisition')
obsf=export(obs,'selection_grid_acquisitions',w,'All 19 arms × 58 acquisitions; peptide strings have modifications stripped but I/L NOT collapsed')
af=export(a,'selection_grid_summary',w,'Medians and quartiles; no inference of optimality')
matched=a[a.kind.eq('length')|a.arm.eq('tp0.30')].sort_values('min_length').copy()
mf=export(matched,'selection_length_matched_pool58',w,'Minimum-nine point uses independently audited tp0.30 run, actual upstream 58 inputs; old ml9_pool97 stays in complete table only')
for c in ['peptides_q01','razor_selected']:
 med=obs.groupby('arm')[c].median();check('Selection median '+c,np.allclose(a.set_index('arm')[c+'_median'].sort_index(),med.sort_index()),'Recomputed from 1,102 acquisition-arm observations')
fig,ax=plt.subplots(2,2,figsize=(7.1,6.0));ax=ax.ravel();real=q[q.arm.eq('REAL')]
ax[0].plot(real.length,real.fraction_matched*100,'o-',c=COL[0],label='Authentic')
for _,g in q[~q.arm.eq('REAL')].groupby('arm'):ax[0].plot(g.length,g.fraction_matched*100,c=COL[1],alpha=.6)
ax[0].plot([],[],c=COL[1],label='Five scrambled seeds');ax[0].set(xlabel='Peptide length (residues)',ylabel='Distinct queries matching (%)',xlim=(8.5,50.5));fig.legend(*ax[0].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.947),ncol=2,frameon=False,fontsize=7);title(ax[0],'a','Exact reservoir lookup')
inset=ax[0].inset_axes([.45,.15,.5,.20])
for _,g in q[~q.arm.eq('REAL')].groupby('arm'):inset.plot(g.length,g.fraction_matched*100,'o-',ms=2,c=COL[1],alpha=.6)
inset.set(xlim=(8.8,13.2),ylim=(-.3,10),xticks=[9,11,13],yticks=[0,5,10]);inset.tick_params(labelsize=7);inset.set_title('Scrambled: 9–13 residues',fontsize=7)
ax[1].plot(real.length,real.distinct_queries,'o-',c=COL[0]);ax[1].set(xlabel='Peptide length (residues)',ylabel='Distinct authentic queries',xlim=(8.5,50.5));title(ax[1],'b','Length-specific denominator')
for axis,kind,x in [(ax[2],'length','min_length'),(ax[3],'score','top_percent')]:
 g=(matched if kind=='length' else a[a.kind.eq(kind)]).sort_values(x);xx=g[x].to_numpy()*(100 if kind=='score' else 1)
 axis.errorbar(xx,g.peptides_q01_median/1000,yerr=np.array([g.peptides_q01_median-g.peptides_q01_q1,g.peptides_q01_q3-g.peptides_q01_median])/1000,fmt='o-',c=COL[0],label='Accepted peptides')
 axis.plot(xx,g.razor_selected_median/1000,'s--',c=COL[1],label='Selected proteins (median)');axis.set(xlabel='Minimum lookup length (residues)' if kind=='length' else 'Top-scoring predictions (%)',ylabel='Count (thousands)',ylim=(0,48));axis.legend(frameon=False,fontsize=7,loc='upper right')
 axis.axvline(9 if kind=='length' else 30,c='grey',lw=.8,ls=':')
 title(ax[2],'c','Length sensitivity');title(ax[3],'d','Score-fraction sensitivity')
for l in 'ab':panel('SELECTION_SPECIFICITY',l,'Exact length-specific lookup',f,'length,distinct_queries,matching_queries,fraction_matched','all six arms; displayed lengths 9–50','query','SM2,SM6','code/PRODUCTIONFILTER.py; code/frozen_lookup_plot.py','P024')
for l,k in [('c','length'),('d','score')]:panel('SELECTION_SPECIFICITY',l,'Selection grid',mf if k=='length' else af,'min_length,top_percent,*_median,*_q1,*_q3',('all matched rows; minimum-nine point uses tp0.30' if k=='length' else 'kind=score')+'; underlying '+obsf,'acquisition; n=58 evaluated and upstream inputs per displayed arm','SM2,SM6','code/selection_sensitivity_producer.py','P022–P025; proposed selection-sensitivity sentence','all medians and quartiles recomputed; tp0.30 actual input roster and all 58 raw accepted-peptide totals independently audited; no new search')
save(fig,'SELECTION_SPECIFICITY','Evidence selection and exact-lookup specificity','a,b, Length-specific lookup fractions and authentic-query counts in the worked healthy-stool acquisition; orange traces are five composition-preserving scrambled sets, enlarged for lengths 9–13 in the inset. Denominators are distinct queries. c,d, Selection sensitivity with 58-acquisition upstream pools and evaluated rosters. The nine-residue point in c reuses the verified top-30% run in d. Points show median accepted peptide counts with exclusive-quartile ranges and median selected protein counts. Counts strip modifications without collapsing I/L. Vertical lines mark production settings, not independently optimized choices. All 1,102 original arm–acquisition records, including the older 97-pool reference omitted from the curves, are supplied. Its final FASTAs match the verified 58-pool run in all 58 acquisitions.')

# Error-control architectures remain separate estimands.
fusion=read('original_S24/*fusion_summary.tsv');fusion=fusion[~fusion.cohort.eq('CAMPI_MOCK')].copy()
c=read('current_CAMPI_fusion/pooled_by_seed.tsv').rename(columns={'pooled_diagnostic':'diagnostic'});c['cohort']='CAMPI';fusion=pd.concat([c,fusion],ignore_index=True)
ff=export(fusion,'fusion_pooled_current','current CAMPI pooled_by_seed + human cohorts from S24','Ratios of summed acquisition-level components; not distinct pooled unions')
ind=read('current_CAMPI_fusion/independent_per_acquisition.tsv');human=read('original_S24/*fusion_recount.tsv');human=human[~human.cohort.eq('CAMPI_MOCK')];indf=export(pd.concat([ind,human],ignore_index=True),'fusion_acquisition_components','current CAMPI independent recount + current-paper human controls','All thresholds, length populations, unresolved categories and score-order sensitivity')
check('Fusion arithmetic',np.allclose(fusion.diagnostic,(fusion.paired_numerator+fusion.n_unpaired_accepted)/fusion.n_canonical),'Every cohort/seed/q/length stratum')
for keys,g in ind.groupby(['seed','q_threshold','min_residues','max_residues']):
 r=c.set_index(['seed','q_threshold','min_residues','max_residues']).loc[keys]
 check('CAMPI fusion sum '+str(keys),int(g.n_canonical.sum())==int(r.n_canonical) and int(g.paired_numerator.sum())==int(r.paired_numerator),'15 acquisitions, summed denominator/numerator')
don=read('original_S02/*donor_seed_means.tsv');df=export(don,'donor_entrapment_seed_means',src('original_S02/*donor_seed_means.tsv'),'Acquisition-level mean across three seeds; self-control and donor-selected strata')
copy(src('original_S02/*donor_recount.tsv'))
sim=read('original_S06/*ordinary_outcomes.tsv');sf=export(sim,'simulation_mixture_outcomes',src('original_S06/*ordinary_outcomes.tsv'),'All 72 arm-mixture rows; false discoveries are known from synthetic construction')
u=read('simulation/uncertainty.tsv');uf=export(u,'simulation_uncertainty',src('simulation/uncertainty.tsv'),'Whole-mixture bootstrap intervals; 12 units/regime; 20,000 resamples in frozen producer')
check('Synthetic FDP identity',np.allclose(sim.realised_FDP,sim.false_discoveries/sim.discoveries),'False/discovery counts, all 72 rows')
fig,axes=plt.subplots(3,2,figsize=(7.1,8.4));axes=axes.ravel()
for i,lo in enumerate([5,9]):
 g=fusion[fusion.q_threshold.eq(.01)&fusion.min_residues.eq(lo)&fusion.max_residues.eq(50)].copy();g['label']=g.cohort+'\n'+g.seed.astype(str).str[-2:]
 axes[i].bar(np.arange(len(g)),g.paired_component*100,label='Paired',color=COL[0]);axes[i].bar(np.arange(len(g)),g.unpaired_component*100,bottom=g.paired_component*100,label='Unresolved / unpaired',color=COL[1]);axes[i].set_xticks(np.arange(len(g)),g.label,rotation=45,ha='right');axes[i].set_ylabel('Fusion-linked diagnostic (%)');title(axes[i],'ab'[i],str(lo)+'–50 residues');axes[i].legend(fontsize=7,frameon=False)
 panel('ERROR_CONTROL_ARCHITECTURE','ab'[i],'Fusion-linked component decomposition',ff,'cohort,seed,paired_component,unpaired_component','q_threshold=.01,min_residues='+str(lo)+',max_residues=50; components in '+indf,'sum of acquisition counts within cohort/seed','SM6','code/audit_fusion.py','P042; proposed control-scope sentence','formula and CAMPI pooled components independently recounted from component tables')
g=fusion[fusion.q_threshold.eq(.01)&fusion.min_residues.eq(5)&fusion.max_residues.eq(8)]
dots(axes[2],[(co,g[g.cohort.eq(co)].diagnostic) for co in g.cohort.unique()],'Fusion-linked diagnostic (%)',True);title(axes[2],'c','Short-peptide sensitivity (5–8)')
panel('ERROR_CONTROL_ARCHITECTURE','c','Short accepted-peptide population',ff,'cohort,seed,diagnostic','q_threshold=.01,min_residues=5,max_residues=8','seed; two/cohort','SM6','code/audit_fusion.py','Proposed short-peptide limitation')
g=don[don.q_threshold.eq(.01)];dots(axes[3],[(('Self' if t=='T0_self' else 'Donor')+'\n'+str(lo)+'–50',g[g.tier.eq(t)&g.min_length.eq(lo)].paired_estimate) for t in ['T0_self','T1_within_group'] for lo in [5,9]],'Conditional paired estimate (%)',True);title(axes[3],'d','Peptide-record donor assay');axes[3].tick_params(axis='x',labelsize=7)
panel('ERROR_CONTROL_ARCHITECTURE','d','Explicit-peptide conditional assay',df,'tier,min_length,target,paired_estimate','q_threshold=.01; all self and within-group records','target acquisition; three seeds averaged','SM6','code/CONTROL.py','P042')
for j,reg in enumerate(['baseline','stress']):
 ax=axes[j+4];g=sim[sim.regime.eq(reg)];dots(ax,[(arm,g[g.arm.eq(arm)].realised_FDP) for arm in ['fixed','independent','self']],'Realised FDP (%)',True)
 for k,arm in enumerate(['fixed','independent','self']):
  r=u[u.regime.eq(reg)&u.endpoint.eq(arm)].iloc[0];ax.errorbar(k+.27,r['mean']*100,yerr=[[100*(r['mean']-r.bootstrap_low)],[100*(r.bootstrap_high-r['mean'])]],fmt='s',ms=3,c='black',capsize=3)
 ax.axhline(1,c='grey',ls=':',lw=1);title(ax,'ef'[j],reg.capitalize()+' synthetic mixtures')
 panel('ERROR_CONTROL_ARCHITECTURE','ef'[j],'Known-truth simulation',sf+';'+uf,'regime,arm,unit_id,realised_FDP; mean,bootstrap_low,bootstrap_high','regime='+reg,'independent synthetic mixture; n=12','SM6','code/simulation_evaluate.py; code/simulation_collect.py','P044','FDP counts independently recomputed; bootstrap producer inspected')
save(fig,'ERROR_CONTROL_ARCHITECTURE','What each error control measures','a,b, Fusion-linked diagnostics at nominal peptide q≤0.01, including both paired and unresolved contributions; each bar is one seed. The CAMPI data use the current broad UniProt reservoir. c, The nested 5–8-residue population retains the same search q values. These diagnostic ratios are not calibrated production FDR estimates. d, Explicit-peptide self-selected (five) and donor-selected (58) acquisition controls, each averaged over three seeds. e,f, Known-truth baseline and stress simulations (12 independent mixtures per regime). Dots are mixtures; horizontal bars are medians; squares and intervals are means and 95% whole-mixture bootstrap intervals. The dotted line marks nominal 1%. Controls differ in selection dependence and are not pooled into one error-rate claim.')

# Recount member support from accepted sequence membership observations, including zeros.
members=read('CAMPI/member_summary.tsv');arms=['FL','FL_R','FL_R_M','REF','MG'];names=members.member.unique();records=[];sets={}
obsdir=P/'source_data/CAMPI/member_observations'
for f in sorted(obsdir.glob('*.tsv')):
 t=pd.read_csv(f,sep='\t');copy(f,'source_data/CAMPI_membership')
 for arm in arms:
  a=t[t.arms.fillna('').str.split(';').map(lambda z:arm in z)]
  for lo in [5,9]:
   a1=a[a.peptide.str.len().between(lo,50)];sets[(f.stem,arm,lo)]=set(a1.peptide)
   for member in names:
    v=set(a1.loc[a1.n_members.eq(1)&a1.members.eq(member),'peptide']);records.append(dict(sample=f.stem,arm=arm,min_length=lo,member=member,specific_peptides=len(v)));sets[(f.stem,arm,lo,member)]=v
r=pd.DataFrame(records);rf=export(r,'CAMPI_member_acquisition_counts','15 delivered CAMPI_membership/*.tsv','Member-specific within the supplied eight-proteome reference only; zeros retained')
pooled=[]
for arm in arms:
 for lo in [5,9]:
  for member in names:
   g=r[r.arm.eq(arm)&r.min_length.eq(lo)&r.member.eq(member)];pooled.append(dict(arm=arm,min_length=lo,member=member,pooled_specific_peptides=len(set().union(*(sets[(s,arm,lo,member)] for s in g['sample']))),acquisitions_with_support=int((g.specific_peptides>0).sum()),total_acquisitions=len(g)))
p=pd.DataFrame(pooled);pf=export(p,'CAMPI_member_support_recount','15 delivered CAMPI_membership/*.tsv','Actual counts and observed acquisition support; no recovery classification')
j=p.merge(members,on=['arm','min_length','member'],suffixes=('_recount','_original'));check('Member support full recount',np.array_equal(j.pooled_specific_peptides_recount,j.pooled_specific_peptides_original)&np.array_equal(j.acquisitions_with_support,j.acquisitions_ge1),'80 arm/length/member strata from peptide observations; prior >=10 rule not used')
rec=read('CAMPI/paired_reference_recovery.tsv');recf=export(rec,'CAMPI_reference_recovery',src('CAMPI/paired_reference_recovery.tsv'),'Acquisition-paired identity overlap, reciprocal gain/loss and count ratio')
for x in rec.itertuples():
 aa=sets[(x.sample,x.arm,x.min_length)];bb=sets[(x.sample,'REF',x.min_length)]
 check('CAMPI identity '+x.sample+x.arm+str(x.min_length),len(aa)==x.peptides and len(bb)==x.reference and len(aa&bb)==x.overlap,'Exact set comparison of accepted observations')
ms=read('original_S03/*within_MSFragger_comparisons.tsv');mf=export(ms,'CAMPI_within_MSFragger',src('original_S03/*within_MSFragger_comparisons.tsv'),'Within-engine FastaLake recovery of REF/MG peptide identities; NOT between-engine overlap')
fig,axes=plt.subplots(2,2,figsize=(7.1,7.0));axes=axes.ravel()
for i,lo in enumerate([5,9]):
 sub=p[p.arm.eq('FL')&p.min_length.eq(lo)].set_index('member').loc[names]
 ax=axes[i];ax.barh(np.arange(8),sub.pooled_specific_peptides,color=COL[i]);ax.set_yticks(np.arange(8),[n.split()[0][0]+'. '+' '.join(n.split()[1:]) for n in names],fontsize=7);ax.set_xscale('log');ax.set_xlabel('Pooled member-specific\npeptide count');ax.invert_yaxis()
 for k,x in enumerate(sub.itertuples()):ax.text(max(x.pooled_specific_peptides,1)*1.08,k,f'{x.pooled_specific_peptides:,}; {x.acquisitions_with_support}/15',va='center',fontsize=7)
 ax.set_xlim(.8,1e7);title(ax,'ab'[i],f'FastaLake, {lo}–50 residues')
 panel('CAMPI_RECOVERY','ab'[i],'Eight-member actual support',pf,'member,pooled_specific_peptides,acquisitions_with_support,total_acquisitions',f'arm=FL,min_length={lo}; all-arm observations: '+rf,'distinct peptide union; acquisition support out of 15','SM6','code/CAMPI_aggregate.py; code/build_revision.py','P029–P035','independent set recount from all 15 accepted-membership tables')
g=rec[rec.min_length.eq(5)];dots(axes[2],[(arm,g[g.arm.eq(arm)].reference_recovery) for arm in arms],'REF identities recovered (%)',True);axes[2].set_ylim(0,105);title(axes[2],'c','Sage recovery')
g=ms[ms.min_length.eq(5)];dots(axes[3],[(co,g[g.comparator.eq(co)].reference_recovery) for co in ['REF','MG']],'Reference identities recovered (%)',True);axes[3].set_ylim(0,105);title(axes[3],'d','MSFragger recovery')
panel('CAMPI_RECOVERY','c','Matched REF identity recovery',recf,'sample,arm,reference,overlap,reference_recovery','min_length=5','acquisition; n=15/arm','SM3,SM6','code/CAMPI_aggregate.py','P033','exact accepted-identity intersections independently recomputed')
panel('CAMPI_RECOVERY','d','Independent-engine within-engine recovery',mf,'sample,comparator,reference,intersection,reference_recovery','min_length=5','acquisition; n=15/comparator','SM6','code/SECONDSEARCH.json; frozen engine collector locator in source registry','P036')
save(fig,'CAMPI_RECOVERY','Reference recovery and member-specific evidence','a,b, Actual pooled member-specific peptide counts in the FastaLake arm, with the number of acquisitions contributing at least one such peptide (labels: count; acquisitions/15). Specificity is relative to the eight supplied reference proteomes, not universal species specificity. The 9–50-residue list is nested within the accepted 5–50 list and retains its q values. No ten-peptide recovery criterion is applied. c, Within-acquisition recovery of REF identities by each Sage search arm; FL_R and FL_R_M append REF or REF+MG after selection. d, Within-MSFragger recovery of REF and MG identities by FastaLake. Dots are acquisitions; bars are medians. All five arms, lengths, zeros and acquisition-level member counts are supplied.')

# Common search: identity, group counts, and precision remain different endpoints.
ident=read('original_S08/*identity_comparison.tsv');idf=export(ident,'common_search_identity',src('original_S08/*identity_comparison.tsv'),'All arms/tiers/lengths; matched identity gains and losses')
cv=read('common_search/protein_cv_pairs.tsv');cvf=export(cv,'common_search_pair_CV',src('common_search/protein_cv_pairs.tsv'),'Population-SD pair CV, dependent acquisition pairs, both row populations')
patients=pd.concat([read('original_S11/*'+co+'_participants.tsv').assign(cohort=co) for co in ['PREDICT','CAPSCAN']]);patf=export(patients,'common_search_participants','S11 participant tables','Participant summaries; independent searches against common FASTA, not joint LFQ')
fig,axes=plt.subplots(2,2,figsize=(7.1,6.4));axes=axes.ravel()
g=ident[ident.min_length.eq(5)];dots(axes[0],[(k.replace('_','\n'),v.original_identity_recovery) for k,v in g.groupby('test')],'Original identities recovered (%)',True);title(axes[0],'a','Healthy stool, 58 acquisitions')
for i,(co,g) in enumerate(patients.groupby('cohort')):
 axes[1].scatter(g.peptides_percent_difference,g.directlfq_groups_percent_difference,s=12,alpha=.6,color=COL[i],label=co)
axes[1].axhline(0,c='grey',lw=.6);axes[1].axvline(0,c='grey',lw=.6);axes[1].set(xlabel='Accepted peptides: common vs own (%)',ylabel='DirectLFQ groups: common vs own (%)');axes[1].legend(frameon=False);title(axes[1],'b','Clinical common-FASTA tradeoffs')
for i,t in enumerate(['TechnicalReplicates','WorkflowReplicates']):
 g=cv[cv.protein_class.eq('Microbial')&cv.tier.eq(t)];dots(axes[i+2],[(de.replace('A_arm_alone','Own rows').replace('B_matched_both_arms','Matched')+'\n'+a.replace('joint_corrected','joint'),g[g.definition.eq(de)&g.arm.eq(a)].median_pair_cv_pct) for de in ['A_arm_alone','B_matched_both_arms'] for a in ['original','joint_corrected']],'Median pair CV (%)');title(axes[i+2],'cd'[i],['Technical replicates (66 pairs)','Preparation replicates (105 pairs)'][i])
panel('COMMON_SEARCH_TRADEOFFS','a','Canonical identity retention',idf,'sample,test,original_identity_recovery','min_length=5; all four comparator arms','acquisition; n=58','SM7','code/QUANT.py','P047–P049')
panel('COMMON_SEARCH_TRADEOFFS','b','Identity versus quantified groups',patf,'participant,peptides_percent_difference,directlfq_groups_percent_difference','all 93 PREDICT and 15 CapScan participant summaries','participant','SM7','common-search clinical producer [exact locator pending]','P047; proposed clinical common-search sentence')
for l,t in [('c','TechnicalReplicates'),('d','WorkflowReplicates')]:panel('COMMON_SEARCH_TRADEOFFS',l,'Protein-group pair precision',cvf,'definition,arm,sample_a,sample_b,median_pair_cv_pct','protein_class=Microbial,tier='+t,'dependent acquisition pair','SM7','code/PROTEINCV.py','P049','ddof and row-eligibility verified in producing code; raw intensity replay pending')
save(fig,'COMMON_SEARCH_TRADEOFFS','Common search spaces: identity, yield and precision','a, Healthy-stool original peptide identities retained by repeat, common-FASTA individual, original joint and corrected joint searches. b, Participant-level changes in peptide counts and directLFQ groups for common-FASTA individual searches in PREDICT (93) and CapScan (15). c,d, Microbial protein-group precision under original and corrected joint search: rows present in each arm separately versus rows matched across both arms. Each dot is an acquisition pair (66 technical or 105 preparation pairs); pairs sharing an acquisition are dependent. Pair CV uses population SD/mean. This differs from the separately supplied all-replicate, sample-SD CV. Bars are medians. Increased group counts do not by themselves establish improved identity recovery or precision.')

# Functional attribution: preserve unresolved allocations and conditional labels.
coverage=[]
for co in ['HSTOOL','PREDICT','CAPSCAN']:
 f=src('allocation/'+co+'_KO_all.tsv.gz');t=pd.read_csv(f,sep='\t',index_col=0);copy(f,'source_data/allocation')
 total=t.sum(axis=0);missing=t.loc['FunctionallyUnassigned'];check(co+' allocation matrix',np.isfinite(t.values).all() and (t.values>=0).all() and (total>0).all(),'Nonnegative allocated and remainder rows; per-acquisition denominator includes remainder')
 coverage.extend(dict(cohort=co,sample=s,total_intensity=float(total[s]),unassigned_intensity=float(missing[s]),allocated_fraction=float(1-missing[s]/total[s])) for s in t.columns)
coverage=pd.DataFrame(coverage);covf=export(coverage,'KO_allocation_coverage','source_data/allocation/*_KO_all.tsv.gz','Each acquisition weighted equally in display; denominator allocated + remainder')
leaders=read('original_S09/*function_tier_summary.tsv');lf=export(leaders,'healthy_stool_leader_diversity',src('original_S09/*function_tier_summary.tsv'),'All eligible functions, rarefied to eight resolved acquisitions')
cells=read('original_S19/*contributor_cells.tsv');cellf=export(cells,'CapScan_function_contributors',src('original_S19/*contributor_cells.tsv'),'Every eligible participant-function cell; resolved and unresolved shares distinguished')
agree=read('original_S10/*SUMMARY.tsv');agf=export(agree,'healthy_stool_grouping_leader_sensitivity',src('original_S10/*SUMMARY.tsv'),'All grouping, feature-subset and resolution definitions')
capagree=read('original_S10/*s10c_rule_comparison.tsv');capagree=capagree[capagree.source.str.startswith('A33_')];caf=export(capagree,'CapScan_grouping_leader_sensitivity','S10 current A33 rows only','All four eligibility rules under both current grouping arms')
fig,axes=plt.subplots(3,2,figsize=(7.1,8.4));axes=axes.ravel()
dots(axes[0],[(co,coverage[coverage.cohort.eq(co)].allocated_fraction) for co in ['HSTOOL','PREDICT','CAPSCAN']],'Accepted feature intensity allocated to KO (%)',True);title(axes[0],'a','Functional coverage across cohorts')
dots(axes[1],[(k.replace('Replicates',''),v.expected_distinct_leaders_in_8_resolved_acquisitions) for k,v in leaders.groupby('tier')],'Expected distinct leaders in 8 acquisitions');axes[1].tick_params(axis='x',rotation=35);title(axes[1],'b','Healthy-stool leader diversity')
for i,(col,label) in enumerate([('leader_share','Leader / species-resolved signal'),('species_resolved_fraction','Resolved / non-host functional signal')]):
 for participant,g in cells.groupby('participant'):
  v=np.sort(g[col].dropna()*100);axes[i+2].plot(v,np.arange(1,len(v)+1)/len(v),alpha=.4,lw=.65)
 axes[i+2].set(xlabel=label+' (%)',ylabel='Cumulative fraction of functions',xlim=(0,100));title(axes[i+2],'cd'[i],['Conditional leader share','Resolution coverage'][i])
g=agree[agree.subset.eq('all_eligible')]
for i,method in enumerate(['original','simplified']):
 sub=g[g.method.eq(method)];axes[4].scatter(np.arange(len(sub))+i*.15,sub.equal_KO_agreement*100,label=method,c=COL[i]);axes[4].set_xticks(np.arange(len(sub))+.075,[t.replace('Replicates','')+'\n'+f for t,f in zip(sub.tier,sub['filter'])],rotation=25,ha='right')
axes[4].set_ylabel('Equal-function leader agreement (%)');axes[4].legend(frameon=False);title(axes[4],'e','Grouping and resolution sensitivity')
for i,(source,g) in enumerate(capagree.groupby('source')):axes[5].plot(np.arange(len(g)),g.donor_balanced_agreement*100,'o-',c=COL[i],label=source.split('__')[-1]);axes[5].set_xticks(np.arange(len(g)),g.rule.str.replace('_',' '),rotation=25,ha='right')
axes[5].set_ylabel('Donor-balanced leader agreement (%)');axes[5].legend(frameon=False);title(axes[5],'f','CapScan eligibility sensitivity')
for l,ep,f,cols,sel,unit,code,anchor in [('a','KO allocation coverage',covf,'allocated_fraction,total_intensity,unassigned_intensity','all acquisitions','acquisition','code/PEPBUILD.py; code/selection.py','P050; proposed coverage sentence'),('b','Rarefied leader diversity',lf,'tier,ko,expected_distinct_leaders_in_8_resolved_acquisitions','all eligible functions','function × tier','code/derive_hstool_pc.py','P051'),('c','Conditional leading-species share',cellf,'KO,participant,leader_share','all 20,352 cells; separate participant ECDFs','function × participant','code/derive_capscan_pc.py','P069'),('d','Species resolution of functional signal',cellf,'KO,participant,species_resolved_fraction','all 20,352 cells; separate participant ECDFs','function × participant','code/derive_capscan_pc.py','P069'),('e','Grouping-rule leader sensitivity',agf,'method,tier,filter,equal_KO_agreement','subset=all_eligible','function; equal weighting','code/derive_hstool_pc.py','P051; not P049 own/joint contrast'),('f','Conditional leader sensitivity',caf,'source,rule,donor_balanced_agreement','current A33 source rows','participant-balanced pair summary','code/check_s10c_eligibility.py','P069')]:panel('FUNCTIONAL_ATTRIBUTION',l,ep,f,cols,sel,unit,'SM5,SM8',code,anchor)
save(fig,'FUNCTIONAL_ATTRIBUTION','Functional coverage and conditional species attribution','a, Fraction of accepted quantified feature intensity allocated to a KEGG orthologue; the denominator includes the functionally unassigned remainder. Acquisitions are shown separately by cohort. b, Expected leading-species diversity after rarefaction to eight resolved healthy-stool acquisitions. c,d, Per-participant cumulative distributions across eligible CapScan functions distinguish leader share within species-resolved signal from species-resolved share of non-host functional signal (20,352 cells; 1,416 functions; 15 participants). e,f, Sensitivity to grouping and resolution/eligibility rules, using the complete current analyses. Species labels are representative labels retained only with a compatible same-KO carrier; they do not identify a unique peptide source. The two denominators in c,d must not be interchanged.')

# Computing: per-acquisition stages restricted by explicit array roster; shared stages as run.
jobs=read('compute/original_compute_jobs.tsv');jobs=jobs[jobs.cohort.isin(['HSTOOL','PREDICT','CAPSCAN'])].copy()
per={'05_persample','06_parsimony','07_sage','09_entrap_shuffle'}
hs=read('compute/hstool_sample_list.tsv');design=read('compute/hstool_original_design.tsv');dc=set(design.loc[design.group.str.startswith('DC_'),'sample_id'])
cap=read('compute/capscan_array_roster.tsv');cp=pd.read_csv(P/'code/evidence/CAPSCAN.tsv',sep='\t');keepcap=set(cp.sample_id)
pr=read('compute/predict_array_roster.tsv');rosters={'HSTOOL':list(hs['sample']),'CAPSCAN':list(cap.sample_id),'PREDICT':list(pr.iloc[:,0])}
rows=[]
for x in jobs.to_dict('records'):
 co=x['cohort'];x['sample']='';x['population_scope']='shared stage as run: '+str(len(rosters[co]))+' acquisitions'
 if x['stage'] in per:
  idx=int(x['jobid'].split('_')[1]);sample=rosters[co][idx]
  if co=='HSTOOL' and sample in dc:continue
  if co=='CAPSCAN' and sample not in keepcap:continue
  x['sample']=sample;x['population_scope']='current manuscript acquisition subset'
 x['memory_GiB']=x.pop('maxrss_gb');x['wall_minutes']=x['elapsed_s']/60;x['allocated_CPU_hours']=x['cpu_alloc_s']/3600;rows.append(x)
cj=read('compute/current_CAMPI_task_records.tsv')
for x in cj.to_dict('records'):
 x['cohort']='CAMPI';x['memory_GiB']=x['sampled_peak_gib'];x['wall_minutes']=x['elapsed_s']/60;x['allocated_CPU_hours']=x['elapsed_s']*x['allocated_logical_cpus']/3600;x['population_scope']='current broad-reservoir production: 15 acquisitions';rows.append(x)
jobs=pd.DataFrame(rows);jf=export(jobs,'compute_task_records','human current-paper producing chain plus September broad-reservoir CAMPI; roster files supplied','Wall time is execution only; shared job scope is explicit; missing memory is NA')
for co,n in [('HSTOOL',58),('PREDICT',306),('CAPSCAN',331),('CAMPI',15)]:check('Computing roster '+co,len(jobs[jobs.cohort.eq(co)&jobs.stage.eq('07_sage')])==n,'One task per current manuscript acquisition')
for f in (P/'source_data/compute').glob('*'):
 if f.name in ['hstool_sample_list.tsv','hstool_original_design.tsv','capscan_array_roster.tsv','predict_array_roster.tsv','original_sacct_raw.tsv','current_CAMPI_task_records.tsv','CAMPI_prediction_execution_records.tsv','INFERENCE_SACCT.tsv'] or 'grouping_' in f.name or 'jobids' in f.name:copy(f,'source_data/compute_records')
stages=['02_evidence','05_persample','06_parsimony','07_sage','08_aggregate','09_entrap_shuffle','10_fdp_shuffle'];summary=[]
for co in ['CAMPI','HSTOOL','PREDICT','CAPSCAN']:
 for st in stages:
  g=jobs[jobs.cohort.eq(co)&jobs.stage.eq(st)]
  summary.append(dict(cohort=co,stage=st,n_tasks=len(g),n_memory_sampled=int(g.memory_GiB.notna().sum()),wall_median_min=g.wall_minutes.median(),memory_median_GiB=g.memory_GiB.median(),allocated_CPU_hours=g.allocated_CPU_hours.sum() if len(g) else np.nan,population_scope=';'.join(g.population_scope.unique())))
s=pd.DataFrame(summary);ssf=export(s,'compute_stage_summary',jf,'No addition of acquisition medians; CPU-hours are allocated, not consumed')
fig,axes=plt.subplots(2,2,figsize=(7.1,6.9));axes=axes.ravel();coh=['CAMPI','HSTOOL','PREDICT','CAPSCAN'];namesst=['Shared lookup','Acquisition lookup','Parsimony','Sage','Aggregation','Shuffle search','Shuffle diagnostic']
for i,(col,lab) in enumerate([('wall_median_min','Median execution time (min)'),('memory_median_GiB','Median sampled memory (GiB)'),('allocated_CPU_hours','Allocated CPU-hours (sum)'),('n_memory_sampled','Tasks with usable memory samples')]):
 mat=s.pivot(index='stage',columns='cohort',values=col).reindex(index=stages,columns=coh);vals=mat.values.astype(float);im=axes[i].imshow(np.ma.masked_invalid(np.log10(1+vals)),cmap='Blues',aspect='auto')
 for y,st in enumerate(stages):
  for x,co in enumerate(coh):
   row=s[s.cohort.eq(co)&s.stage.eq(st)].iloc[0];v=vals[y,x];label=('NA' if not np.isfinite(v) or row.n_tasks==0 else f'{v:.2g}') if i<3 else (f'{row.n_memory_sampled}/{row.n_tasks}' if row.n_tasks else 'NA')
   shade=np.log10(1+v)/np.nanmax(np.log10(1+vals)) if np.isfinite(v) else 0
   axes[i].text(x,y,label,ha='center',va='center',fontsize=8,color='white' if shade>.62 and label!='NA' else 'black')
 axes[i].set_xticks(range(4),['CAMPI\n15','Stool\n58','PREDICT\n306','CapScan\n331']);axes[i].set_yticks(range(7),namesst);title(axes[i],'abcd'[i],lab)
 panel('COMPUTE_RESOURCE_MAP','abcd'[i],lab,ssf,'cohort,stage,'+col+',n_tasks,n_memory_sampled,population_scope','all cells; task-level values: '+jf,'job/task; allocation sum for c','SM11','code/collect_sacct.py; code/build_revision.py','P043; proposed cross-cohort resource sentence','array indices mapped to supplied current rosters; temporal/output lineage limitations recorded separately')
save(fig,'COMPUTE_RESOURCE_MAP','Computational resources across manuscript cohorts','Execution and sampled-memory summaries for the producing CPU workflow. Acquisition stages use 15 CAMPI, 58 healthy-stool, 306 PREDICT and 331 CapScan acquisitions. Shared healthy-stool and CapScan stages retain their actual 97- and 368-acquisition execution scopes. CAMPI uses its broad-reservoir September production chain; its matched shuffle timing records are not supplied here and remain NA. a, Median task execution minutes. b, Median usable sampled memory in GiB (2³⁰ bytes). c, Sum of allocated CPU-hours, distinct from consumed CPU time. d, Number of usable memory samples/number of tasks. Cells below the producer’s 0.01-GiB recording threshold are missing, not zero. Shading is log₁₀(1+value); printed values are untransformed. Queueing, conversion, de novo inference and downstream biological reanalysis are outside these stage timings; available CAMPI inference records are supplied separately. This is not a complete workflow speed comparison.')

# Clinical comparator scope and fixed-dictionary loss, without stale spectrum endpoint.
counts=read('original_S13/*curation_counts.tsv');ccf=export(counts,'PREDICT_curation_counts',src('original_S13/*curation_counts.tsv'),'296 shared acquisitions × four comparator arms')
abl=read('original_S18/*plot_values.tsv');abf=export(abl,'clinical_fixed_dictionary_ablation',src('original_S18/*plot_values.tsv'),'All 578 observations, total and removed intensity, KO support and profile distance')
check('Ablation intensity fractions',np.allclose(abl.fraction_intensity_removed,abl.FL_only_LFQ_intensity/abl.all_LFQ_intensity),'Both cohorts, all observations')
check('Ablation KO fraction',np.allclose(abl.KO_support_lost_percent,100*(abl.all_KOs-abl.shared_KOs)/abl.all_KOs),'Both cohorts, all observations')
fig,axes=plt.subplots(2,2,figsize=(7.1,6.5));axes=axes.ravel()
dots(axes[0],[(k.replace('_','\n'),g.peptides/1000) for k,g in counts.groupby('arm')],'Accepted peptides (thousands)');title(axes[0],'a','PREDICT comparator curation')
for i,(col,label) in enumerate([('intensity_removed_percent','LFQ intensity removed (%)'),('KO_support_lost_percent','KO support lost (%)'),('BrayCurtis','Normalized KO Bray–Curtis distance')]):
 dots(axes[i+1],[(co,g[col]) for co,g in abl.groupby('cohort')],label);title(axes[i+1],'bcd'[i],['Removed peptide evidence','Functional support','Profile change'][i])
panel('CLINICAL_COMPARATOR_SCOPE','a','Comparator construction sensitivity',ccf,'sample,arm,peptides','all four arms, n=296 each','paired acquisition','SM3','code/CLINICAL.py','P054–P056')
for l,col in zip('bcd',['intensity_removed_percent','KO_support_lost_percent','BrayCurtis']):panel('CLINICAL_COMPARATOR_SCOPE',l,'Fixed-dictionary removal of FL-only evidence',abf,'sample,cohort,'+col,'all 296 PREDICT and 282 CapScan matched acquisitions','acquisition; unequal participant contributions','SM5,SM8','fixed-dictionary ablation producer [exact locator pending]','P055/P066; new explicit ablation explanation')
save(fig,'CLINICAL_COMPARATOR_SCOPE','Comparator curation and the contribution of additional evidence','a, Accepted peptide counts in four PREDICT database arms on the same 296 acquisitions: full FastaLake including cohort metagenomes, specimen MG, public-only FastaLake and curated specimen MG. b–d, Removing FastaLake-only peptide evidence from a fixed reporting dictionary: intensity removed, KOs losing all support, and Bray–Curtis distance between normalized all-evidence and shared-only KO profiles. Dots are acquisitions (296 PREDICT; 282 CapScan); bars are medians. This ablation is conditional on the FastaLake dictionary and does not reconstruct a metagenome-only workflow. The earlier full-FastaLake versus pooled-MG spectrum comparison is not used as evidence for the current public-only versus specimen-MG comparison.')

# Histidine: six CapScan motif maps plus both PREDICT boundaries and routes.
mot=read('original_S21/*source_data.tsv');motf=export(mot,'CapScan_histidine_motif_display',src('original_S21/*source_data.tsv'),'840 displayed cells; q values from the complete 32,400-test family')
allhist=read('CAPSCAN_histidine/complete_statistical_families.tsv');allhf=export(allhist,'CapScan_histidine_all_tests',src('CAPSCAN_histidine/complete_statistical_families.tsv'),'Complete planned statistical families, including unavailable p=1 hypotheses')
copy(src('CAPSCAN_histidine/participant_motif_differences.tsv.gz'),'source_data/histidine')
ef=read('PREDICT_histidine/acquisition_enrichment.tsv.gz');ef=ef[ef.participant.notna() & ef.participant.astype(str).ne('')].copy();key=['route','side','origin','residue'];ef['O']=ef.observed/ef.denominator.replace(0,np.nan);ef['E']=ef.expected/ef.denominator.replace(0,np.nan)
sp=ef.groupby(key+['cohort','participant','region','specimen_id'],dropna=False)[['O','E']].mean().reset_index();pt=sp.groupby(key+['cohort','participant'],dropna=False)[['O','E']].mean().reset_index();pt['difference']=pt.O-pt.E
ptf=export(pt,'PREDICT_histidine_participant_enrichment','PREDICT acquisition_enrichment.tsv.gz; equal specimens within participant','Observed-minus-source-expected residue frequency; no acquisition pooling')
copy(src('PREDICT_histidine/acquisition_enrichment.tsv.gz'),'source_data/histidine');st=read('PREDICT_histidine/complete_statistical_families.tsv');stf=export(st,'PREDICT_histidine_all_tests',src('PREDICT_histidine/complete_statistical_families.tsv'),'Complete 1,560 residue and 12 reciprocal tests')
fig,axes=plt.subplots(4,2,figsize=(7.1,10.8),gridspec_kw={'height_ratios':[1.15,1.15,1.15,.8]});axes=axes.ravel();scale=float(np.abs(mot['mean']).max()*100)
for i,(region,side) in enumerate((r,s) for r in ['proximal','distal','stool'] for s in ['N','C']):
 g=mot[mot.contrast.eq(region)&mot.side.eq(side)];a=g.pivot(index='residue',columns='position',values='mean').sort_index();qmat=g.pivot(index='residue',columns='position',values='q').reindex_like(a)
 im=axes[i].imshow(a.values*100,vmin=-scale,vmax=scale,cmap='RdBu_r',aspect='auto');axes[i].set_xticks(range(len(a.columns)),[f'{int(x):+}' for x in a.columns]);axes[i].set_yticks(range(len(a)),a.index,fontsize=7);axes[i].set_xlabel('Position relative to boundary');title(axes[i],'abcdef'[i],region.capitalize()+'; '+side+' boundary')
 yy,xx=np.where(qmat.values<.05);axes[i].scatter(xx,yy,s=3,c='black');fig.colorbar(im,ax=axes[i],fraction=.04,pad=.02,label='Observed − expected (pp)')
 panel('HISTIDINE_CONTEXT_REPLICATION','abcdef'[i],'Histidine-conditioned sequence context',motf,'residue,position,mean,q',f'contrast={region},side={side}; full family: '+allhf,'participant; 14 complete-region CapScan donors','SM9','code/collect_revision4.py; code/histidine_helpers.py','P081–P082; correct previous P084 anchor')
for j,route in enumerate(['denovo_top30','sage_unspecific_bothq']):
 g=pt[pt.route.eq(route)&pt.origin.eq('All')&pt.residue.eq('H')];ax=axes[j+6];dots(ax,[(s,g[g.side.eq(s)].difference) for s in ['N','C']],'H enrichment (percentage points)',True)
 for i,side in enumerate(['N','C']):
  r=st[st.route.eq(route)&st.origin.eq('All')&st.residue.eq('H')&st.side.eq(side)&st.endpoint.eq('source_enrichment')].iloc[0];vals=g[g.side.eq(side)].difference.dropna();check('PREDICT H participant mean '+route+side,len(vals)==r.n and np.isclose(vals.mean(),r['mean']),'Specimen means then equal participant means; n='+str(len(vals)))
  ax.errorbar(i+.28,r['mean']*100,yerr=[[100*(r['mean']-r.CI95_low)],[100*(r.CI95_high-r['mean'])]],fmt='s',c='black',ms=3,capsize=3)
 ax.axhline(0,c='grey',lw=.6);title(ax,'gh'[j],'PREDICT: '+('de novo top 30%' if j==0 else 'unspecific Sage; both q gates'))
 panel('HISTIDINE_CONTEXT_REPLICATION','gh'[j],'Independent cohort source-conditioned H enrichment',ptf+';'+stf,'participant,side,difference; mean,CI95_low,CI95_high,q','route='+route+',origin=All,residue=H','participant; n=93 per side/route','SM9','code/collect_revision4.py','P083','93 participant values independently reconstructed from acquisition observations')
save(fig,'HISTIDINE_CONTEXT_REPLICATION','Histidine context and replication across cohorts','a–f, CapScan histidine-conditioned motif deviations from source-protein opportunities for proximal fluid, distal fluid and stool, separately at N and C boundaries. Colours show participant-mean observed-minus-expected frequencies in percentage points; dots mark q<0.05 in the full 32,400-test motif family, not only these 840 displayed cells. The conditioned position is omitted. g,h, PREDICT histidine enrichment for de novo and unspecific-search evidence, retaining both N and C boundaries. Dots are 93 participant means after averaging repeated acquisitions within specimens and specimens within participants; bars are medians, and squares/intervals are means/95% participant-bootstrap intervals. These associations do not establish the cleavage enzyme or when cleavage occurred.')

# Measured spectra with complete observed peaks and bounded alternatives supplied.
case=read('original_S22/*independent_case_review.tsv');casef=export(case,'PREDICT_spectrum_case_review',src('original_S22/*independent_case_review.tsv'),'All four measured cases and independent bounded alternative review')
fig,axes=plt.subplots(4,1,figsize=(7.1,8.8))
for i,(_,r) in enumerate(case.iterrows()):
 cid=r.case_id.replace('PREDICT_','');folder=P/'source_data/spectra'/cid
 peaks=pd.read_csv(folder/'processed_peaks.tsv',sep='\t');frag=pd.read_csv(folder/'matched_fragments.tsv',sep='\t');files=[copy(f,'source_data/spectra/'+cid) for f in folder.iterdir() if f.is_file()]
 ax=axes[i];norm=peaks.intensity.max();ax.vlines(peaks.mz,0,peaks.intensity/norm*100,color='#b5b5b5',lw=.6)
 for series,g in frag.groupby('series'):
  ax.vlines(g.observed_mz,0,g.intensity/norm*100,color=COL[0] if series=='b' else COL[1],lw=.9,label=series+' ions')
  for j,x in enumerate(g.itertuples()):
   if x.intensity/norm>.08:ax.text(x.observed_mz,min(105,x.intensity/norm*100+4+j%2*6),f'{series}{x.ordinal}'+('²⁺' if x.charge==2 else ''),rotation=90,fontsize=7,ha='center',va='bottom')
 ax.set(xlabel='m/z',ylabel='Relative intensity (%)',ylim=(0,137));title(ax,'abcd'[i],cid+'  '+r.peptide+'  ('+r.origin+', '+r.side+' boundary)');ax.legend(frameon=False,loc='upper right',fontsize=7)
 check('Fragment ppm '+cid,np.allclose((frag.observed_mz-frag.calculated_mz)/(frag.calculated_mz-1.0072764)*1e6,frag.ppm,atol=.000001),'Producer ppm uses proton-subtracted m/z (neutral mass/charge); inspected frozen search_revision2.py:91; not a new sequence identification')
 panel('HISTIDINE_MEASURED_SPECTRA','abcd'[i],'Measured spectrum '+cid,';'.join(files)+';'+casef,'mz,intensity; observed_mz,calculated_mz,series,ordinal,charge,ppm','case '+cid,'one selected measured spectrum; not an independent cohort replicate','SM6,SM9','spectrum reconstruction producer [exact locator in PROVENANCE; audit open]','P083; correct previous P084 anchor','ppm arithmetic checked; manual sequence/alternative audit remains open')
save(fig,'HISTIDINE_MEASURED_SPECTRA','Measured PREDICT histidine-boundary spectra','Four selected PREDICT spectra from the unspecific-search analysis. Gray peaks show the processed spectrum; blue and orange peaks mark matched b and y ions. Intensities are normalized to each spectrum’s maximum processed peak. Thin gray leader lines connect displaced ion labels to their peak positions. The source-data package includes observed and processed peaks, all matched fragments, PSM identifiers, source FASTA/RAW provenance and bounded sequence alternatives for every case. N-boundary histidine is inferred from the compatible source context, whereas C-boundary histidine belongs to the observed peptide. Selection is illustrative; four spectra are not four independent estimates of a cleavage mechanism. The bounded alternative screen does not exclude every possible sequence or modification.')

# Cross-omics: all layer pairs, feature sensitivities and method choices together.
native=read('native/SOURCE_DATA_correlations.tsv');nf=export(native,'native_crossomics_correlations',src('native/SOURCE_DATA_correlations.tsv'),'All native pair estimates, both feature populations and complete per-cohort BH families')
ret=read('native/SOURCE_DATA_retained_signal_by_profile.tsv');retf=export(ret,'native_retained_signal',src('native/SOURCE_DATA_retained_signal_by_profile.tsv'),'Per-profile signal retained by prevalence filter, separate from sequence resolution upstream')
gene=read('native/SOURCE_DATA_gene_sequence_coverage.tsv');genf=export(gene,'native_sequence_resolution',src('native/SOURCE_DATA_gene_sequence_coverage.tsv'),'Upstream gene-to-exact-sequence mapping coverage, not feature-filter retention')
method=read('method_comparison/SOURCE_DATA_method_comparisons.tsv');methf=export(method,'native_distance_method_comparisons',src('method_comparison/SOURCE_DATA_method_comparisons.tsv'),'All 220 estimates, both adjustments and tails; displayed subsets specified per panel')
fun=read('multiomics/functional/profile_correlations_checked.tsv');capfun=read('multiomics/functional/capscan_profile_correlations_with_DNA.tsv');fuf=export(fun,'PREDICT_functional_crossomics',src('multiomics/functional/profile_correlations_checked.tsv'),'All 28 tests and original/simplified comparison');cfuf=export(capfun,'CapScan_functional_crossomics',src('multiomics/functional/capscan_profile_correlations_with_DNA.tsv'),'16 primary tests plus seven same-donor reference estimates')
fig,axes=plt.subplots(3,2,figsize=(7.1,9.5));axes=axes.ravel()
layer={'metaG':'DNA','metaT':'RNA','microbial_metaP':'Non-host protein','host_metaP':'Host protein','metabolites':'Metabolites'}
for i,co in enumerate(['PREDICT','CAPSCAN']):
 g=fun.copy() if co=='PREDICT' else capfun[~capfun.status.str.contains('not displayed')].copy();g['pair']=g.layer_a.map(layer)+' / '+g.layer_b.map(layer);val='rho_simplified' if i==0 else 'rho';mat=g.pivot(index='pair',columns='level',values=val).reindex(columns=['OG','KO','Pathway']);im=axes[i].imshow(np.ma.masked_invalid(mat.values),cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto')
 for y in range(len(mat)):
  for x in range(3):axes[i].text(x,y,'NA' if pd.isna(mat.iloc[y,x]) else f'{mat.iloc[y,x]:.2f}',ha='center',va='center',fontsize=7)
 axes[i].set_yticks(range(len(mat)),mat.index,fontsize=7);axes[i].set_xticks(range(3),['OG','KO','Pathway']);title(axes[i],'ab'[i],co+' functional profiles')
 panel('CROSSOMICS_ROBUSTNESS','ab'[i],'Functional profile concordance',fuf if i==0 else cfuf,'layer_a,layer_b,level,'+val+',p,q' if i==0 else 'layer_a,layer_b,level,rho,p,q_family','all primary functional layer pairs; no omission of nulls','77 PREDICT participants; 12 CapScan DNA or 14 non-DNA participants','SM8','code/crossomics_pc.py; current DNA extension producer [locator pending]','P056/P062/P070/P075')
 g=native[native.cohort.eq(co)].copy();g['pair']=g.layer_a+' / '+g.layer_b;paired=g.pivot(index='pair',columns='analysis',values='rho');a=[x for x in paired if x.startswith('primary')][0];b=[x for x in paired if x.startswith('sensitivity')][0]
 for k,row in enumerate(paired.itertuples(index=False)):pass
 for k,(_,row) in enumerate(paired.iterrows()):axes[i+2].plot([row[a],row[b]],[k,k],c='grey',lw=.8)
 axes[i+2].scatter(paired[a],range(len(paired)),c=COL[0],s=14,label='Prevalence ≥10%');axes[i+2].scatter(paired[b],range(len(paired)),c=COL[1],s=14,label='All observed');axes[i+2].set_yticks(range(len(paired)),paired.index.str.replace(' protein',''),fontsize=7);axes[i+2].set(xlim=(-1,1),xlabel='Native distance Spearman correlation');axes[i+2].axvline(0,c='grey',lw=.5);title(axes[i+2],'cd'[i],co+' feature sensitivity')
 panel('CROSSOMICS_ROBUSTNESS','cd'[i],'Native all-observed vs prevalence sensitivity',nf,'layer_a,layer_b,analysis,rho,p,q','cohort='+co+'; all pairs','77 PREDICT profiles or 60 CapScan profiles from 12 participants','SM8','code/native_pyramids.py','P076; proposed native sensitivity explanation')
 g=method[method.cohort.eq(co)&method.adjustment.eq('unadjusted' if i==0 else 'region_adjusted')].copy();g['pair']=g.layer_a+' / '+g.layer_b;g['method']=g.features.str.replace('all_observed','All').str.replace('prevalence10','10%').str.replace('wuyts_recipe','Wuyts')+'\n'+g.distance.str.replace('hellinger','Hell')+'\n'+g.outer_correlation.str.replace('spearman','S').str.replace('pearson','P');mat=g.pivot(index='pair',columns='method',values='coefficient');im=axes[i+4].imshow(mat.values,cmap='RdBu_r',vmin=-1,vmax=1,aspect='auto')
 codes=['M'+str(j+1) for j in range(len(mat.columns))];export(pd.DataFrame({'column':codes,'method':mat.columns.str.replace('\n',' / ')}),co+'_method_display_key',methf,'Exact display column key: filter / distance / outer correlation')
 axes[i+4].set_xticks(range(len(mat.columns)),codes,fontsize=7,rotation=90);axes[i+4].set_yticks(range(len(mat)),mat.index.str.replace(' protein',''),fontsize=7);title(axes[i+4],'ef'[i],co+' method sensitivity');fig.colorbar(im,ax=axes[i+4],fraction=.04,pad=.02,label='Correlation')
 panel('CROSSOMICS_ROBUSTNESS','ef'[i],'Distance-method sensitivity',methf,'features,distance,outer_correlation,layer_a,layer_b,coefficient,p,p_positive_tail,q','cohort='+co+',adjustment='+('unadjusted' if i==0 else 'region_adjusted'),'same frozen complete-stack roster','SM8','code/experiment.py','P076','all 220 delivered; display uses appropriate primary adjustment per cohort')
fig.legend(*axes[2].get_legend_handles_labels(),loc='upper center',bbox_to_anchor=(.5,.957),ncol=2,frameon=False,fontsize=7)
save(fig,'CROSSOMICS_ROBUSTNESS','Cross-cohort functional and native-profile robustness','a,b, Functional-profile concordance across all tested layer pairs and annotation levels; unavailable combinations remain NA. c,d, Native-profile concordance under ≥10% prevalence versus all-observed features. e,f, Full method grid on the same frozen roster, showing PREDICT unadjusted and CapScan region-adjusted results (Hell, Hellinger geometry; rank, within-profile rank geometry; P/S, Pearson/Spearman outer correlation). Colours retain negative and near-zero estimates. PREDICT uses 77 participants; native CapScan uses 60 participant-region profiles from 12 participants. Functional CapScan pairs without DNA use 14 participants and pairs with DNA use 12. Distances sharing profiles are dependent; permutation units are participants, with whole five-region blocks preserved for CapScan. M1–M4 use ≥10% prevalence and M5–M8 all observed features, each ordered Hellinger/Pearson, Hellinger/Spearman, rank/Pearson and rank/Spearman; M9/M10 use the adapted Wuyts filter with rank/Pearson or rank/Spearman. All 220 method estimates, both adjustments, complete p/q families, sequence-resolution coverage and retained-signal fractions are supplied. Exploratory MOFA, including unfinished fits and unstable factors, is retained in a named table.')

# Freeze exact sources used in this redraw and separately preserve all current table-only analyses.
code_map={'PRODUCTIONFILTER.py':'code/evidence/PRODUCTIONFILTER.py','CONTROL.py':'code/evidence/CONTROL.py','SECONDSEARCH.json':'code/evidence/SECONDSEARCH.json','QUANT.py':'code/evidence/QUANT.py','PROTEINCV.py':'code/evidence/PROTEINCV.py','CLINICAL.py':'code/evidence/CLINICAL.py','PEPBUILD.py':'code/evidence/PEPBUILD.py','CAMPI_aggregate.py':'code/CAMPI/aggregate.py','audit_fusion.py':'code/current_CAMPI_fusion/audit_fusion.py','simulation_evaluate.py':'code/simulation/evaluate.py','simulation_collect.py':'code/simulation/collect.py','collect_sacct.py':'source_data/compute/original_collect_sacct.py','native_pyramids.py':'source_data/native/native_pyramids.py','experiment.py':'source_data/method_comparison/experiment.py'}
for f in (P/'code/peptide_allocation').glob('*.py'):code_map[f.name]=str(f.relative_to(P))
for f in (P/'code/histidine').glob('*.py'):code_map[f.name]=str(f.relative_to(P))
code_map['frozen_lookup_plot.py']=str(src('original_S01/*build_figures.py').relative_to(P));code_map['check_s10c_eligibility.py']=str(src('original_S10/*check_s10c_eligibility.py').relative_to(P))
for name,orig in code_map.items():
 f=P/orig;to=D/'code'/name;shutil.copy2(f,to);FILES.append(dict(file=str(to.relative_to(D)),origin=str(f),role='frozen producing code; absolute upstream paths preserved, not a portable full pipeline',sha256=sha(to),rows=''))
TABLES=[]
for label,groups,why in [('CAMPI_FULL_SUPPORT',['original_S03','original_S04','original_S05'],'All second-engine, accession/group support and identity endpoints'),('COMMON_SEARCH_FULL',['original_S08','original_S11'],'Completeness, all-replicate CV, transfer and clinical participant/acquisition tradeoffs'),('COMPARATOR_GROUPS',['original_S14'],'Independent grouping and fixed-dictionary endpoints remain distinct'),('EDS_DESCRIPTIVE',['original_S15'],'Current v11 fixed feature set, all samples and missingness'),('EDS_PREDICTION_SENSITIVITY',['original_S16','original_S17'],'Prediction, internal validation, depth and context controls together'),('TARGETED_FUNCTION',['original_S20'],'All targeted and housekeeping tests, including adjusted nulls')]:
 for group in groups:
  for f in sorted((P/'source_data'/group).glob('*')):
   if f.suffix in ['.tsv','.csv','.xlsx','.gz']:
    dest=copy(f,'source_data/tables/'+label);TABLES.append(dict(table=label,file=dest,reason=why,status='retained current-paper source; table-specific scientific audit in paragraph ledger'))
for f in (P/'source_data/multiomics').rglob('*'):
 if f.is_file() and any(x in f.name.lower() for x in ['eds','dysbio','fixed_feature','25_','s25','predict_participant']):
  dest=copy(f,'source_data/tables/EDS_DESCRIPTIVE');TABLES.append(dict(table='EDS_DESCRIPTIVE',file=dest,reason='Later v11 descriptive EDS data',status='retained'))
for f in (P/'source_data/method_comparison').glob('SOURCE_DATA_MOFA*'):
 dest=copy(f,'source_data/tables/MOFA_COMPLETE');TABLES.append(dict(table='MOFA_COMPLETE',file=dest,reason='All 39 fits, unfinished fits, stability, variance and covariates',status='exploratory; not confirmatory'))
for f in (P/'source_data/common_search').glob('*'):
 if f.suffix in ['.tsv','.gz']:copy(f,'source_data/tables/COMMON_SEARCH_FULL')
for f in (P/'source_data/native').glob('*distance_matrices.npz'):copy(f,'source_data/native_distances')
for f in (P/'source_data/CAPSCAN_histidine').glob('*.tsv'):copy(f,'source_data/histidine')
pd.DataFrame(TABLES).to_csv(D/'TABLE_REGISTRY.tsv',sep='\t',index=False)
pd.DataFrame(PANELS).to_csv(D/'PANEL_SOURCE_CODE_MAP.tsv',sep='\t',index=False)
pd.DataFrame(FILES).drop_duplicates('file').to_csv(D/'FILE_MANIFEST.tsv',sep='\t',index=False)
pd.DataFrame(CHECKS).to_csv(D/'validation/NUMERICAL_CHECKS.tsv',sep='\t',index=False)
(D/'CAPTIONS.md').write_text('\n\n'.join('## '+k+'\n\n'+v for k,v in CAPTIONS.items())+'\n')
shutil.copy2(__file__,D/'code/build_revision.py')
# Package only this builder's selected inputs; never traverse an analysis archive.
for row in FILES:
 f=Path(str(row['origin']))
 if f.is_file() and P in f.parents:INPUTS.add(f)
INPUTS.add(P/'code/evidence/CAPSCAN.tsv')
for f in sorted(INPUTS):
 to=D/'rebuild_inputs'/f.relative_to(P);to.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(f,to)
shutil.copy2(__file__,D/'rebuild_inputs/build_revision.py')
(D/'rebuild_figures.sh').write_text('#!/usr/bin/env bash\nset -euo pipefail\ncd -- "$(dirname -- "$0")"\nexport FASTALAKE_SUPPLEMENT_OUTPUT="$PWD/rebuilt"\npython rebuild_inputs/build_revision.py\n')
print('Built',len(CAPTIONS),'figures,',len(PANELS),'panel mappings;',len(CHECKS),'checks passed; table records',len(TABLES))
