from pathlib import Path
import json,hashlib,datetime,copy
import numpy as np
import pandas as pd
W=Path(__file__).resolve().parent;R=W.parent  # C347: R = spectrum_split347 (CapScan 347 run)
PRODREV=Path('LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/split_half_repair/mzml_recovery_2026-09-27/scientific_review_2026-09-27')
OLD=PRODREV.parents[1];BASE=OLD.parents[2]  # C347: production locations unchanged
ROSTER347=Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/roster/participant15_added.tsv')
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def js(p):return json.load(open(p))
m=js(R/'MANIFEST.json');audit=js(R/'audit/FULL_AUDIT.json');assert audit['passed'] and audit['verified']==347  # C347
sources=[R/'MANIFEST.json',R/'audit/FULL_AUDIT.json',R/'audit/FULL_SOURCE_DATA.tsv',R/'PROTOCOL.md',PRODREV/'ACQUISITION_REVIEW.tsv',ROSTER347]  # C347
meta={}
for co in ['PREDICT','CAPSCAN']:
 p=BASE/f'supplements/revision_2026-09-24/delivery/source_data/clinical_common/{co}_per_acquisition.tsv';sources.append(p)
 df=pd.read_csv(p,sep='\t');assert not df['sample'].duplicated().any();meta[co]=df.set_index('sample')
# Cross-check CapScan participant labels against separately archived metadata.
p=OLD.parent/'RAW_CAPSCAN_METADATA.tsv';sources.append(p);raw=pd.read_csv(p,sep='\t').set_index('sample_id')
for s,row in meta['CAPSCAN'].iterrows():assert str(int(row['participant']))==str(int(raw.loc[s,'subject_id']))
# C347: participant-15 capsule acquisitions from the fixed roster (participant=subject_id, specimen=main:<workbook Main>)
add=pd.read_csv(ROSTER347,sep='\t');assert len(add)==16 and (add.subject_and_type_agree==True).all()
meta['CAPSCAN']=pd.concat([meta['CAPSCAN'],pd.DataFrame(dict(participant=add.subject_id.values,specimen=['main:'+str(x) for x in add.workbook_Main],region=add.region.values),index=add.sample_id.values)])
assert not meta['CAPSCAN'].index.duplicated().any()
counts={}
for co in ['CAPSCAN']:  # C347: CapScan halves of this run
 p=R/f'inputs/{co}/COMPLETE.json';sources.append(p)
 for r in js(p)['samples']:counts[(co,r['sample'])]=r['counts']
rows=[];configuration_checks=0
for t in m['tasks']:
 co,s=t['cohort'],t['sample'];out=R/'outputs'/co/s;receipt=js(out/'COMPLETE.json');parity=js(out/'PARITY.json')
 sets={a:set((out/a/'ACCEPTED.txt').read_text().splitlines()) for a in ['FULL','HALF_SELF','HALF_DISJOINT']}
 configs=[]
 for a in sets:
  assert sha(out/a/'ACCEPTED.txt')==receipt['metrics'][a]['accepted_sha256']
  cfg=js(out/a/'config.json');configs.append(cfg)
 baseline=js(t['config']);baseline.pop('version',None)
 def normalize(c):
  c=copy.deepcopy(c);c['database'].pop('fasta');c.pop('mzml_paths');c.pop('output_directory',None);return c
 assert all(normalize(c)==normalize(baseline) for c in configs)
 assert len({tuple(c['mzml_paths']) for c in configs})==1
 configuration_checks+=1
 proteins={}
 for p,a in [('even','D'),('odd','S')]:
  pp=list((out/'selection'/p/'parsimony').glob('*_stats.json'));assert len(pp)==1
  proteins[a]=js(pp[0])['selected_proteins']
 c=counts[(co,s)];extra=dict(participant=None,specimen=None)
 if co in meta:
  mr=meta[co].loc[s];extra=dict(participant=mr['participant'],specimen=mr['specimen'])
 for minlen in [5,9]:
  F,S,D=[{p for p in sets[a] if minlen<=len(p)<=50} for a in ['FULL','HALF_SELF','HALF_DISJOINT']]
  for a,x in [('FULL',F),('HALF_SELF',S),('HALF_DISJOINT',D)]:assert len(x)==receipt['metrics'][a][f'peptides_{minlen}_50']
  rows.append(dict(index=t['index'],cohort=co,sample=s,min_length=minlen,full=len(F),self=len(S),disjoint=len(D),D_over_S=len(D)/len(S) if S else np.nan,D_over_F=len(D)/len(F) if F else np.nan,S_over_F=len(S)/len(F) if F else np.nan,lost_from_self=len(S-D),gained_vs_self=len(D-S),shared_with_self=len(D&S),fraction_self_identities_retained=len(D&S)/len(S) if S else np.nan,disjoint_fraction_in_full=len(D&F)/len(D) if D else np.nan,even_predictions=c['even'],odd_predictions=c['odd'],prediction_ratio=c['even']/c['odd'],even_selected_proteins=proteins['D'],odd_selected_proteins=proteins['S'],selected_protein_ratio=proteins['D']/proteins['S'],conversion_recovered=False,**extra))  # C347: the five recovered inputs are PREDICT; no CapScan input recovered
df=pd.DataFrame(rows)
# C347: index = production index for the 331 (379-709), 710+k for the 16 added acquisitions
prod=pd.read_csv(PRODREV/'ACQUISITION_REVIEW.tsv',sep='\t',float_precision='round_trip');pidx=prod[prod.cohort=='CAPSCAN'].drop_duplicates('sample').set_index('sample')['index'].to_dict()
new=sorted(set(df['sample'])-set(pidx));assert len(new)==16 and len(pidx)==331
df['index']=[pidx.get(s,710+new.index(s) if s in new else -1) for s in df['sample']]
# C347: identity of the 331 CapScan rows with production (all columns except index)
old=prod[prod.cohort=='CAPSCAN'].set_index(['sample','min_length']).sort_index();cur=df[df['sample'].isin(pidx)].set_index(['sample','min_length']).sort_index()
cols=[c for c in old.columns if c!='index'];diff=[c for c in cols if not old[c].equals(cur[c].astype(old[c].dtype) if c in cur else None)]
(W/'IDENTITY_331.json').write_text(json.dumps(dict(rows_compared=len(old),identical_columns=[c for c in cols if c not in diff],differing_columns=diff,status='PASS' if not diff else 'DIFFERS'),indent=2)+'\n')
# C347: HSTOOL, CAMPI_MOCK and PREDICT rows reused unchanged from production
df=pd.concat([prod[prod.cohort!='CAPSCAN'],df],ignore_index=True).sort_values(['index','min_length']).reset_index(drop=True)
df.to_csv(W/'ACQUISITION_REVIEW.tsv',sep='\t',index=False)
summary=[]
for (co,length),g in df.groupby(['cohort','min_length']):
 r=dict(cohort=co,min_length=int(length),n=len(g),zero_full=int((g.full==0).sum()),zero_self=int((g.self==0).sum()),zero_disjoint=int((g.disjoint==0).sum()),D_lower_S=int((g.disjoint<g.self).sum()),D_equal_S=int((g.disjoint==g.self).sum()),D_higher_S=int((g.disjoint>g.self).sum()))
 for col in ['D_over_S','D_over_F','S_over_F','full','self','disjoint','lost_from_self','gained_vs_self','fraction_self_identities_retained','disjoint_fraction_in_full','prediction_ratio','selected_protein_ratio']:
  v=g[col];r['median_'+col]=float(v.median());r['q25_'+col]=float(v.quantile(.25));r['q75_'+col]=float(v.quantile(.75));r['min_'+col]=float(v.min());r['max_'+col]=float(v.max())
 summary.append(r)
pd.DataFrame(summary).to_csv(W/'COHORT_SUMMARY.tsv',sep='\t',index=False)
# Descriptive equal-participant sensitivity: median of within-participant median ratios.
# No acquisition-level p-values or independent-observation confidence intervals.
participant=[];ps=[];coverage=[]
for co in ['PREDICT','CAPSCAN']:
 for length in [5,9]:
  g=df[(df.cohort==co)&(df.min_length==length)];valid=g.dropna(subset=['participant'])
  coverage.append(dict(cohort=co,min_length=length,total=len(g),mapped=len(valid),participants=int(valid.participant.nunique()),missing_samples=g[g.participant.isna()]['sample'].tolist()))
  for label,sub in [('all_mapped',valid),('exclude_recovered',valid[~valid.conversion_recovered])]:
   by=sub.groupby('participant')[['D_over_S','D_over_F','S_over_F','fraction_self_identities_retained']].median();by['n_acquisitions']=sub.groupby('participant').size()
   for pid,row in by.iterrows():participant.append(dict(cohort=co,min_length=length,subset=label,participant=pid,**row.to_dict()))
   ps.append(dict(cohort=co,min_length=length,subset=label,participants=len(by),acquisitions=len(sub),median_D_over_S=float(by.D_over_S.median()),median_D_over_F=float(by.D_over_F.median()),median_S_over_F=float(by.S_over_F.median()),participants_D_lower_S=int((by.D_over_S<1).sum())))
pd.DataFrame(participant).to_csv(W/'PARTICIPANT_MEDIANS.tsv',sep='\t',index=False);pd.DataFrame(ps).to_csv(W/'PARTICIPANT_SUMMARY.tsv',sep='\t',index=False)
(W/'METADATA_COVERAGE.json').write_text(json.dumps(coverage,indent=2)+'\n')
result=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),passed=True,acquisitions_checked=configuration_checks,accepted_list_hashes_checked=3*configuration_checks,length_views=[5,9],cohort_summary=summary,participant_summary=ps,metadata_coverage=coverage,clinical_metadata_scope='Existing clinical-common acquisition roster; full exact-name joins, no inferred IDs; CapScan participant labels independently cross-checked against archived RAW_CAPSCAN_METADATA. Four PREDICT unknown IDs excluded only from participant summary.',source_hashes={str(p):sha(p) for p in sources},script_sha256=sha(__file__),limitations=['one split direction only: even selects / odd scores','prediction identities can recur across halves','half-specific lookup cutoffs and selected databases differ','separate native Sage peptide q-values per search arm','no full-workflow FDR or peptide correctness estimate','public-only clinical arms; cannot extrapolate to assembly-inclusive search','five PREDICT alternative-conversion inputs; exclusion sensitivity provided'])
(W/'REVIEW_CHECK.json').write_text(json.dumps(result,indent=2)+'\n')
print(pd.DataFrame(summary)[['cohort','min_length','n','median_D_over_S','median_D_over_F','D_lower_S','D_higher_S','median_prediction_ratio','median_selected_protein_ratio']].to_string(index=False));print(pd.DataFrame(ps).to_string(index=False))
