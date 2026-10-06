from pathlib import Path
from collections import Counter
import ast, hashlib, json, os
import ahocorasick
import pandas as pd
import numpy as np
from analyze import HERE, FL, TASKS, load, dump, meta, verify, canonical, table

m=load(HERE/'MECHANISM_INPUTS_v2.json')
for r in m['inputs']+m['code']:verify(r)
old=FL/'MM_review_analysis/A14_4_3_campi_gap/frozen/search.py'
# Import only the audited pure functions, without invoking its CLI/global RUN.
tree=ast.parse(old.read_text());names={'fasta','auto','compatible','scan','fixture'}
pure=ast.Module(body=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in names],type_ignores=[])
assert len(pure.body)==len(names)
exec(compile(pure,str(old),'exec'),globals());fixture()
base=load(HERE/'INPUTS.json');figure=pd.read_csv(HERE/'inputs/figure_metrics.tsv',sep='\t').set_index(['sample','strategy'])
features={}
for arm in ['FL','FLR','FLRM']:
 p=FL/'simplified_grouping_results/manuscript_cohorts_2026-09-10_2200/cohorts'/TASKS[arm]/'direct/feature_assignments.tsv.gz'
 a=pd.read_csv(p,sep='\t',usecols=['sample','canonical_peptide','study_group','intensity'])
 assert a.intensity.gt(0).all()
 features[arm]={s:g for s,g in a.groupby('sample')}
out=HERE/'mechanism_v2';out.mkdir(exist_ok=False)
ordered=sorted(base['tasks'],key=lambda t:(t['sample'] not in ['S01','S11_Fraction4'],t['sample']))
all_gains=[];all_fates=[];all_quant=[];all_groups=[];all_upstream=[]
for index,t in enumerate(ordered):
 sample=t['sample'];frames={};ids={}
 for arm in ['FL','FLR','FLRM']:
  r=t['arms'][arm]['search'];verify(r)
  a=pd.read_csv(r['path'],sep='\t',usecols=['scannr','peptide','rank','charge','label','peptide_q','hyperscore'])
  assert a['rank'].eq(1).all() and a.scannr.is_unique
  a['canonical']=a.peptide.map(canonical);a['accepted']=a.label.eq(1)&a.peptide_q.le(.01)&a.canonical.str.len().between(5,50)
  a=a.set_index('scannr');frames[arm]=a;ids[arm]=set(a.loc[a.accepted,'canonical'])
  assert len(ids[arm])==int(figure.loc[(sample,arm),'peptides'])
  assert features[arm][sample].study_group.nunique()==int(figure.loc[(sample,arm),'quantified_groups'])
 ref=pd.read_csv(t['gap']['path'],sep='\t');verify(t['gap'])
 qsets={arm:set(features[arm][sample].canonical_peptide) for arm in frames}
 wanted=set.union(*ids.values(),*qsets.values(),set(ref.peptide))
 fasta_path=load(t['arms']['FL']['config']['path'])['database']['fasta'];verify(m['fastas'][sample])
 present,digest,audit=scan(fasta_path,wanted)
 assert audit['sha256']==m['fastas'][sample]['sha256']
 assert ref.final_parents.eq(ref.peptide.map(lambda p:present[p])).all()
 assert ref.final_digest_parents.eq(ref.peptide.map(lambda p:digest[p])).all()
 assert all(digest[p]>0 for p in ids['FL'])
 def availability(p):return 'digest_compatible_in_FL_fasta' if digest[p] else 'sequence_present_digest_incompatible' if present[p] else 'sequence_absent_from_FL_fasta'
 for dest in ['FLR','FLRM']:
  gained=ids[dest]-ids['FL'];lost=ids['FL']-ids[dest]
  for cat,n in Counter(availability(p) for p in gained).items():all_gains.append(dict(sample=sample,destination=dest,availability=cat,gained_peptides=n,total_gained=len(gained),total_lost=len(lost),net_change=len(ids[dest])-len(ids['FL'])))
  assert len(gained)-len(lost)==int(figure.loc[(sample,dest),'peptides']-figure.loc[(sample,'FL'),'peptides'])
  y=frames[dest];y=y[y.accepted&y.canonical.isin(gained)];x=frames['FL'].reindex(y.index)
  same=x.label.eq(1)&x.peptide.eq(y.peptide)&x.charge.eq(y.charge)
  categories=pd.Series('',index=y.index,dtype=object)
  categories.loc[same]='same_form_below_threshold_FL'
  categories.loc[x.label.eq(1)&~same&x.accepted.eq(True)]='different_accepted_FL_target'
  categories.loc[x.label.eq(1)&~same&~x.accepted.eq(True)]='different_unaccepted_FL_target'
  categories.loc[x.label.eq(-1)]='FL_decoy_winner'
  categories.loc[x.label.isna()]='no_reported_FL_candidate'
  assert categories.ne('').all();assert not (same&x.accepted.eq(True)).any()
  z=pd.DataFrame(dict(availability=y.canonical.map(availability),previous_FL_result=categories,raw_unchanged=same&x.hyperscore.eq(y.hyperscore)))
  for (avail,cat),g in z.groupby(['availability','previous_FL_result']):all_fates.append(dict(sample=sample,destination=dest,availability=avail,previous_FL_result=cat,spectra=len(g),same_form_raw_unchanged=int(g.raw_unchanged.sum()),total_new_peptide_spectra=len(z)))
  newq=qsets[dest]-qsets['FL']
  def qclass(p):
   if p in qsets['FL']:return 'shared_quantified'
   if p not in ids[dest]:return 'quantified_not_in_destination_accepted_IDs'
   if p in ids['FL']:return 'already_identified_FL_newly_quantified'
   return 'newly_identified_and_quantified'
  for (cat,avail),n in Counter((qclass(p),availability(p)) for p in newq).items():all_quant.append(dict(sample=sample,destination=dest,category=cat,availability=avail,peptides=n,total_destination_only_quantified=len(newq)))
  d=features[dest][sample].copy();d['support_class']=d.canonical_peptide.map(qclass)
  grouping=[]
  for key,g in d.groupby('study_group'):
   grouping.append(dict(group=key,support=';'.join(sorted(set(g.support_class))),intensity=float(g.intensity.sum())))
  gg=pd.DataFrame(grouping);assert len(gg)==int(figure.loc[(sample,dest),'quantified_groups']);assert np.isclose(gg.intensity.sum(),d.intensity.sum(),rtol=1e-10)
  for support,g in gg.groupby('support'):all_groups.append(dict(sample=sample,destination=dest,support=support,groups=len(g),intensity=float(g.intensity.sum()),total_groups=len(gg)))
 # Existing REF-supported losses: prediction filters and candidate-selection stages.
 missing=ref[~ref.accepted_FL.astype(bool)]
 for stage,g in missing.groupby('stage'):
  all_upstream.append(dict(sample=sample,stage=stage,peptides=len(g),direct_prediction_any=int(g.direct_prediction_any.sum()),direct_prediction_qualifying=int(g.direct_prediction_qualifying.sum()),direct_prediction_lookup=int(g.direct_prediction_lookup.sum()),add_ref_recovered=int(g.accepted_FL_R.sum())))
 print('MECHANISM PASS',sample,flush=True)
 if index==1:dump(out/'PILOT_GATE.json',dict(status='PASS',samples=[r['sample'] for r in ordered[:2]],A14_availability_reproduced=True,primary_and_group_counts_reproduced=True))
for name,rows in [('gain_availability.tsv',all_gains),('gain_spectrum_fates.tsv',all_fates),('new_quantified_peptides.tsv',all_quant),('destination_group_support.tsv',all_groups),('upstream_prediction_stages.tsv',all_upstream)]:table(out/name,rows)
dump(out/'COMPLETE.json',dict(status='PASS_COMPUTATION_REQUIRES_REVIEW',acquisitions=15,addback_arms=2,job_id=os.getenv('SLURM_JOB_ID'),current_primary_and_group_counts_reproduced=True,A14_final_availability_reproduced=True,scope='15 current SIHUMIx acquisitions; trace evidence, not biochemical truth',outputs=[meta(p) for p in out.iterdir() if p.is_file()]))
