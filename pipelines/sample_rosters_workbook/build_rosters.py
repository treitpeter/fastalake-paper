"""sample_rosters_workbook: build rosters. Moved from the original analysis as it ran, apart from path scrubbing; see README.md."""
from pathlib import Path
import csv,hashlib,json,shutil
import numpy as np,pandas as pd
B=Path('LAB_ROOT')
W=B/'00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23'
M=W/'mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24'
O=M/'review_decisions/M01_sample_rosters/capscan';O.mkdir(parents=True,exist_ok=True)
A=B/'MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19';Q=A/'results/CAPSCAN__simplified/source_data';X=A/'results/crossomics_repaired__simplified/source_data';OLD=B/'analyses/enterotypes_2026-09-07';MG=B/'analyses/capscan_metag_analysis_2026-09-10';DNA=B/'fable_5.1_review_2026-09-16/scoping_review_2026-09-19/figure_comments_2026-09-21/capscan_dna';N=W/'figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23';A16=B/'MM_review_analysis/A16_4_5_public_comparison';REG=['duodenum','jejunum','ileum','colon','stool']
SOURCES={};CHECKS={};members=[];excluded=[];index=[]
def sha(p):
 h=hashlib.sha256()
 with p.open('rb') as f:
  for c in iter(lambda:f.read(8<<20),b''):h.update(c)
 return h.hexdigest()
def source(p):
 p=Path(p)
 if str(p) not in SOURCES:SOURCES[str(p)]={'source_file':str(p),'sha256':sha(p),'bytes':p.stat().st_size}
 return str(p)
def read(p,**k):source(p);return pd.read_csv(p,sep='\t',**k)
def js(p):source(p);return json.loads(p.read_text())
def write(d,name):d.to_csv(O/name,sep='\t',index=False)
def check(k,v):CHECKS[k]=bool(v);assert v,k
D=read(Q/'Figure5/design.tsv',dtype={'sample_id':str,'specimen_id':str,'participant':str}).rename(columns={'sample_id':'acquisition_id','specimen_id':'specimen_id_biology','participant':'participant_id'})
C=read(OLD/'CAPSCAN_specimen_design.tsv',dtype={'sample_id':str,'specimen':str,'subject_id':str}).rename(columns={'sample_id':'acquisition_id','specimen':'specimen_id_crossomics','subject_id':'participant_crossomics','group':'region_crossomics'})
base=D[['acquisition_id','specimen_id_biology','participant_id','region','aliquot','device_id','aliquot_main_id']].merge(C[['acquisition_id','specimen_id_crossomics','participant_crossomics','region_crossomics']],on='acquisition_id',validate='1:1')
check('base331_292_15',len(base)==331 and base.specimen_id_biology.nunique()==292 and base.participant_id.nunique()==15)
check('exact_acquisition_link_agrees_participant_region',(base.participant_id==base.participant_crossomics).all() and (base.region==base.region_crossomics).all())
check('specimen_identifier_namespaces_one_to_one',base.groupby('specimen_id_biology').specimen_id_crossomics.nunique().max()==1 and base.groupby('specimen_id_crossomics').specimen_id_biology.nunique().max()==1)
base=base.drop(columns=['participant_crossomics','region_crossomics']);base['link_basis']='Exact acquisition_id present in both producing design tables; no identifier stripping or inferred matching.'
write(base,'BASE_ACQUISITIONS.tsv')
S=base.groupby('specimen_id_crossomics',sort=False).first().reset_index();check('specimen_has_single_participant_region',base.groupby('specimen_id_crossomics').participant_id.nunique().max()==1 and base.groupby('specimen_id_crossomics').region.nunique().max()==1)
write(base[['acquisition_id','specimen_id_biology','specimen_id_crossomics','participant_id','region','link_basis']],'IDENTIFIER_CROSSWALK.tsv')
q=read(Q/'Figure5/acquisition_composition_local_candidate.tsv',dtype={'sample_id':str});qids=set(q.sample_id);check('quant330_positive',len(q)==330 and (q.total_LFQ>0).all())
count=read(Q/'Figure5/peptide_counts.tsv');bench=read(A16/'results/all_acquisitions.tsv',dtype={'sample':str});bench=bench[(bench.cohort=='CAPSCAN')&(bench.test_arm=='FL_public')&(bench.reference_arm=='MG')&(bench.min_length==5)].copy();bids=set(bench['sample']);check('public_MG282',len(bench)==282)
manifest=js(A16/'inputs/manifest.json');tasks={r['sample']:r for r in manifest['tasks'] if r['cohort']=='CAPSCAN'};check('A16all331',set(tasks)==set(base.acquisition_id));check('49_exclusions_missing_MG_manifest',all('MG' not in tasks[x]['sets'] for x in set(base.acquisition_id)-bids))
oldbench=read(Q/'Figure5/matched_benchmark.tsv');check('public_and_inclusive_benchmark_roster_equal',set(oldbench.sample_id)==bids)
bb=bench.rename(columns={'sample':'acquisition_id','specimen_id':'specimen_id_A16','subject_id':'participant_A16','group':'region_A16'}).merge(base,on='acquisition_id',validate='1:1');check('benchmark_metadata_agrees',(bb.participant_A16.astype(int).astype(str)==bb.participant_id).all() and (bb.region_A16==bb.region).all());write(bb[['acquisition_id','specimen_id_A16','specimen_id_biology','specimen_id_crossomics','participant_id','region','n_test','n_reference']],'BENCHMARK282_ACQUISITIONS.tsv')
met=read(OLD/'CAPSCAN_metabolite_intensity.tsv',index_col=0);met.index=met.index.astype(str)
og=read(X/'CAPSCAN_noclust/microbial_OG.tsv',index_col=0).reindex(C.acquisition_id)
positive_og=set(og.index[og.sum(1)>0]);validspec=set(base.loc[base.acquisition_id.isin(positive_og),'specimen_id_crossomics']);metspec=set(met.index)&validspec;check('crossomics_valid_met_specimens287',len(metspec)==287)
saved_spec=read(X/'crossomics/CAPSCAN_fullstack_matched_specimens.tsv',dtype=str);check('saved_functional_candidate_specimens_match',set(saved_spec.specimen_main)==metspec)
fcells=read(X/'crossomics/CAPSCAN_fullstack_cells.tsv',dtype=str);fparticipants=set(fcells.subject_id);check('functional70cells14participants',len(fcells)==70 and len(fparticipants)==14)
nspec=read(N/'results/CAPSCAN/protein_metabolite_specimen_roster.tsv',dtype=str).rename(columns={'Unnamed: 0':'specimen_id_crossomics','subject_id':'participant_id','group':'region'});nparticipants=set(nspec.participant_id);ncells=read(DNA/'results/DNA_donor_region_roster.tsv',dtype=str);check('native60cells12participants',len(ncells)==60 and len(nparticipants)==12 and set(ncells.subject_id)==nparticipants)
check('native257specimens_match_rederived',len(nspec)==257 and set(nspec.specimen_id_crossomics)==set(S.loc[S.specimen_id_crossomics.isin(metspec)&S.participant_id.isin(nparticipants),'specimen_id_crossomics']))
runs=read(MG/'manifests/runs.tsv',dtype={'specimen':str,'subject_id':str});active=read(MG/'manifests/active_runs.tsv',dtype={'specimen':str,'subject_id':str});runexcluded=read(MG/'manifests/excluded_runs.tsv',dtype=str);selection=js(MG/'manifests/selection.json');agg=js(MG/'aggregate/aggregation_status.json');check('DNA302_to301',len(runs)==302 and len(active)==301 and set(runs.run_accession)-set(active.run_accession)=={'SRR18585253'} and agg['complete'])
for name in ['runs.tsv','active_runs.tsv','excluded_runs.tsv']:check('selection_sha_'+name,sha(MG/'manifests'/name)==selection['file_sha256'][name])
check('DNA301_completion_receipts',all((MG/'results'/r/'07_eggnog.done.json').exists() for r in active.run_accession))
nruns=read(N/'CAPSCAN_source_roster.tsv',dtype={'id':str,'specimen':str,'participant':str}).rename(columns={'id':'run_accession','participant':'participant_id'});check('native228DNA',len(nruns)==228 and nruns.run_accession.is_unique)
# Reproduce functional DNA eligibility from every level, preserving primary run/specimen IDs.
DNAeligible={};have_by_level={};eligible_run_sets={}
for level in ['OG','KO','Pathway']:
 frame=read(MG/('aggregate/metaG_'+level+'.tsv'),index_col=0);positive=set(frame.index[frame.sum(1)>0]);check('DNA_'+level+'_301_positive',positive==set(active.run_accession))
 e=active[active.run_accession.isin(positive)&active.specimen.isin(metspec)].merge(S[['specimen_id_crossomics','participant_id','region']],left_on='specimen',right_on='specimen_id_crossomics',validate='1:1',suffixes=('_DNA',''))
 check('DNA_'+level+'_join_metadata_agrees',(e.subject_id==e.participant_id).all() and (e.group==e.region).all())
 DNAeligible[level]=e;have_by_level[level]=set(zip(e.participant_id,e.region));eligible_run_sets[level]=set(e.loc[e.participant_id.isin(nparticipants),'run_accession'])
check('native_functional_DNA228_same_runs',all(x==set(nruns.run_accession) for x in eligible_run_sets.values()))
have=set.intersection(*have_by_level.values());missing_cells=[(p,r) for p in sorted(fparticipants,key=int) for r in REG if (p,r) not in have];check('missing_DNA_cells_are5colon7duodenum',set(missing_cells)=={('5','colon'),('7','duodenum')})
# Common row schema: biological and source identifiers remain separate, linked only through explicit designs.
def acqrows(ids):return base[base.acquisition_id.isin(ids)].copy()
def enroll(aid,purpose,unit,rule,acquisitions=None,specimens=None,participants=None,cellpairs=None,dnaruns=None,primary=None,limits='',family=''):
 before=len(members)
 af=acqrows(acquisitions) if acquisitions is not None else base.iloc[0:0].copy()
 sf=S[S.specimen_id_crossomics.isin(specimens)].copy() if specimens is not None else S[S.specimen_id_crossomics.isin(af.specimen_id_crossomics)].copy()
 ps=set(participants if participants is not None else sf.participant_id)
 cp=list(cellpairs if cellpairs is not None else sorted(set(zip(sf.participant_id,sf.region))))
 source_file=source(primary);source_hash=SOURCES[source_file]['sha256']
 def row(kind,entity,**kw):members.append(dict(analysis_id=aid,entity_type=kind,entity_id=str(entity),acquisition_id='',specimen_id_biology='',specimen_id_crossomics='',participant_id='',region='',run_accession='',source_file=source_file,source_sha256=source_hash,inclusion_rule=rule,**kw))
 # Avoid duplicate default keys by update rather than unpack.
 def add(kind,entity,**kw):
  rr=dict(analysis_id=aid,entity_type=kind,entity_id=str(entity),acquisition_id='',specimen_id_biology='',specimen_id_crossomics='',participant_id='',region='',run_accession='',source_file=source_file,source_sha256=source_hash,inclusion_rule=rule);rr.update(kw);members.append(rr)
 for r in af.itertuples():add('acquisition',r.acquisition_id,acquisition_id=r.acquisition_id,specimen_id_biology=r.specimen_id_biology,specimen_id_crossomics=r.specimen_id_crossomics,participant_id=r.participant_id,region=r.region)
 for r in sf.itertuples():add('specimen',r.specimen_id_crossomics,specimen_id_biology=r.specimen_id_biology,specimen_id_crossomics=r.specimen_id_crossomics,participant_id=r.participant_id,region=r.region)
 for p in sorted(ps,key=int):add('participant',p,participant_id=p)
 for p,r in cp:add('participant_region',p+'|'+r,participant_id=p,region=r)
 if dnaruns is not None:
  for r in dnaruns.itertuples():
   known=S[S.specimen_id_crossomics.eq(r.specimen)]
   bio=known.specimen_id_biology.iloc[0] if len(known)==1 else ''
   add('DNA_run',r.run_accession,run_accession=r.run_accession,specimen_id_crossomics=r.specimen,specimen_id_biology=bio,participant_id=str(r.subject_id),region=r.group)
 sub=pd.DataFrame(members[before:]);rosterfile=aid+'.tsv';write(sub,rosterfile)
 index.append(dict(analysis_id=aid,purpose=purpose,biological_unit=unit,inclusion_rule=rule,acquisitions=len(af),specimens=len(sf),participants=len(ps),participant_region_cells=len(cp),DNA_runs=0 if dnaruns is None else len(dnaruns),testing_family=family,source_file=source_file,source_sha256=source_hash,roster_file=str(O/rosterfile),limitations=limits))
ALL=set(base.acquisition_id);P14=set(str(x) for x in range(1,15));Q14=set(base.loc[base.acquisition_id.isin(qids)&base.participant_id.isin(P14),'acquisition_id']);Fspec=set(S.loc[S.specimen_id_crossomics.isin(metspec)&S.participant_id.isin(P14),'specimen_id_crossomics']);Nspec=set(nspec.specimen_id_crossomics);FACQ=set(base.loc[base.specimen_id_crossomics.isin(Fspec),'acquisition_id']);NACQ=set(base.loc[base.specimen_id_crossomics.isin(Nspec),'acquisition_id']);cell14=[(p,r) for p in sorted(P14,key=int) for r in REG];cell12=[(p,r) for p in sorted(nparticipants,key=int) for r in REG]
enroll('CAPSCAN_SEARCH_ALL','All proteomic search acquisitions','Acquisition, nested in specimen and participant','All 331 rows in the producing CapScan design',acquisitions=ALL,primary=Q/'Figure5/design.tsv')
enroll('CAPSCAN_PUBLIC_MG_282','Same-public FastaLake versus specimen-metagenome peptide comparison','Matched acquisition; donor summaries preserve source specimen grouping','A16 CAPSCAN FL_public versus MG, min_length=5; MG comparator exists in task manifest',acquisitions=bids,primary=A16/'results/all_acquisitions.tsv',limits='A16 specimen_id is a participant/device key; preserve it in BENCHMARK282_ACQUISITIONS.tsv. Biology specimen namespace is linked through exact acquisition ID, not substituted into original A16 estimator.')
enroll('CAPSCAN_QUANT_NONEMPTY','Nonempty regional quantitative input','Acquisition','Positive total quantified intensity in the repaired simplified source',acquisitions=qids,primary=Q/'Figure5/acquisition_composition_local_candidate.tsv')
enroll('CAPSCAN_REGIONAL_DESCRIPTIVE15','Host composition summaries including stool-only participant','Participant-region','Nonempty acquisitions averaged within specimen, then participant-region',acquisitions=qids,primary=Q/'Figure5/participant_region_composition_local_candidate.tsv',limits='71 observed cells, including participant 15 stool; not the 70-cell paired design.')
enroll('CAPSCAN_REGIONAL_PAIRED14','Capsule-stool contrasts and own-stool functional geometry','Participant for four paired regions; 70 participant-region profiles','Participants 1–14 with all five regions; nonempty acquisitions, specimen means then region means',acquisitions=Q14,participants=P14,cellpairs=cell14,primary=Q/'Figure5/own_stool_donors.tsv',family='Four host-share comparisons; own-stool label permutation is separate')
enroll('CAPSCAN_FUNCTIONAL_NONDNA14','Functional cross-omics pairs without DNA','Participant, preserving five regions','Valid microbial-OG specimen plus metabolite measurement, and all five participant regions',acquisitions=FACQ,specimens=Fspec,participants=P14,cellpairs=cell14,primary=X/'crossomics/CAPSCAN_fullstack_cells.tsv',family='Seven non-DNA comparisons within the joint 16-test CapScan functional family',limits='Protein acquisition rows are inputs to within-specimen median; feature-specific finite contribution may differ.')
enroll('CAPSCAN_NATIVE12','Native DNA/protein/metabolite profile comparisons','Participant, preserving five regions','Same twelve participants have all five eligible DNA regions; protein/metabolite specimens need not be identical to DNA specimens',acquisitions=NACQ,specimens=Nspec,participants=nparticipants,cellpairs=cell12,dnaruns=active[active.run_accession.isin(nruns.run_accession)],primary=N/'INPUT_MANIFEST.json',family='Six native layer pairs, separately for each prevalence variant',limits='257 protein/metabolite specimens and 228 DNA specimens are separate layer inputs, combined at participant-region; do not assert specimen-level complete cases.')
enroll('CAPSCAN_DNA_COMPLETED301','Completed metagenomic processing runs','DNA sequencing run','Authorized active_runs.tsv after one explicitly documented missing-mate exclusion',participants=set(active.subject_id),cellpairs=sorted(set(zip(active.subject_id,active.group))),dnaruns=active,primary=MG/'manifests/active_runs.tsv',limits='301 distinct DNA specimens; not all have metaproteomics or metabolomics. No proteomic acquisition membership is inferred from DNA sequencing alone.')
enroll('CAPSCAN_FUNCTIONAL_DNA12','Functional DNA versus protein/metabolite profiles','Participant, preserving five regions','Positive OG/KO/pathway DNA profiles; eligible companion specimens; all five regions in twelve participants',acquisitions=NACQ,specimens=Nspec,participants=nparticipants,cellpairs=cell12,dnaruns=active[active.run_accession.isin(nruns.run_accession)],primary=DNA/'results/DNA_donor_region_roster.tsv',family='Nine DNA comparisons within the same 16-test functional CapScan family',limits='Protein/metabolite and DNA specimens can differ within participant-region; native six-test family is separate.')
enroll('CAPSCAN_FUNCTIONAL_NONDNA12_SENS','Matched-participant sensitivity for non-DNA functional pairs','Participant, preserving five regions','Same twelve DNA-complete participants; no DNA layer in the compared pair',acquisitions=NACQ,specimens=Nspec,participants=nparticipants,cellpairs=cell12,primary=DNA/'results/capscan_profile_correlations_with_DNA.tsv',family='Seven sensitivity rows labelled not displayed; no q_family assigned in producing output',limits='Do not count these seven rows as additional members of the displayed 16-test family.')
# Keep sequenced specimen roster independent from proteomic design, with blank biology keys where absent.
dnaall=active[['run_accession','specimen','subject_id','group']].rename(columns={'specimen':'specimen_id_DNA','subject_id':'participant_id_DNA','group':'region_DNA'})
write(dnaall,'DNA301_RUN_SPECIMEN_ROSTER.tsv')
# Exclusions have explicit reference populations and stage-specific reasons.
def excl(aid,kind,entity,reason,primary,**kw):
 row=dict(analysis_id=aid,entity_type=kind,entity_id=str(entity),reason=reason,source_file=source(primary),source_sha256=SOURCES[str(primary)]['sha256']);row.update(kw);excluded.append(row)
for r in base[~base.acquisition_id.isin(bids)].itertuples():excl('CAPSCAN_PUBLIC_MG_282','acquisition',r.acquisition_id,'MG peptide-set/comparator is absent from this acquisition task in A16; public FastaLake result is retained elsewhere.',A16/'inputs/manifest.json',participant_id=r.participant_id,region=r.region)
for r in base[~base.acquisition_id.isin(qids)].itertuples():excl('CAPSCAN_QUANT_NONEMPTY','acquisition',r.acquisition_id,'Zero total quantified intensity; excluded from quantitative profiles, retained among all searched acquisitions.',Q/'Figure5/acquisition_composition_local_candidate.tsv',participant_id=r.participant_id,region=r.region,specimen_id=r.specimen_id_biology)
for aid in ['CAPSCAN_REGIONAL_PAIRED14','CAPSCAN_FUNCTIONAL_NONDNA14']:
 excl(aid,'participant','15','Stool-only in the proteomic design; no capsule regions for paired five-region analysis.',Q/'Figure5/design.tsv')
for p,r in missing_cells:
 for aid in ['CAPSCAN_NATIVE12','CAPSCAN_FUNCTIONAL_DNA12','CAPSCAN_FUNCTIONAL_NONDNA12_SENS']:
  excl(aid,'participant',p,'Excluded as a whole participant: no eligible DNA '+r+' profile after linking positive DNA runs to valid companion specimens; complete five-region blocks required.',DNA/'results/DNA_EXTENSION.json',missing_region=r)
for aid in ['CAPSCAN_NATIVE12','CAPSCAN_FUNCTIONAL_DNA12','CAPSCAN_FUNCTIONAL_NONDNA12_SENS']:
 excl(aid,'participant','15','Stool-only proteomic design; cannot form five-region block.',Q/'Figure5/design.tsv')
for r in runexcluded.itertuples():excl('CAPSCAN_DNA_COMPLETED301','DNA_run',r.run_accession,r.reason,MG/'manifests/excluded_runs.tsv',participant_id=r.subject_id,region=r.group,specimen_id=r.specimen)
for r in S[~S.specimen_id_crossomics.isin(metspec)].itertuples():
 reason='No metabolite row in the producing companion table.' if r.specimen_id_crossomics not in set(met.index) else 'No acquisition with positive microbial OG signal for this specimen.'
 excl('CAPSCAN_FUNCTIONAL_NONDNA14','specimen',r.specimen_id_crossomics,reason,OLD/'CAPSCAN_metabolite_intensity.tsv',participant_id=r.participant_id,region=r.region)
for r in active[~active.run_accession.isin(nruns.run_accession)].itertuples():
 if r.specimen not in metspec:reason='DNA specimen is not in the valid microbial-OG plus metabolite specimen set used by cross-omics.'
 elif r.subject_id not in nparticipants:reason='Participant excluded from complete five-region native/DNA design (5 lacks colon, 7 lacks duodenum, 15 lacks capsule proteomics).'
 else:reason='Not in frozen native source roster; inspect source-specific eligibility.'
 excl('CAPSCAN_NATIVE12','DNA_run',r.run_accession,reason,N/'CAPSCAN_source_roster.tsv',participant_id=r.subject_id,region=r.group,specimen_id=r.specimen)
# Source-specific cell eligibility table is the transparent explanation of 14 ->12.
cellrows=[]
for p,r in cell14:
 for level in ['OG','KO','Pathway']:
  e=DNAeligible[level];sel=e[(e.participant_id==p)&(e.region==r)]
  cellrows.append(dict(participant_id=p,region=r,representation=level,eligible_DNA_runs=len(sel),eligible_DNA_specimens=sel.specimen.nunique(),DNA_region_available=len(sel)>0,in_complete12=p in nparticipants))
write(pd.DataFrame(cellrows),'DNA_CELL_ELIGIBILITY.tsv');write(pd.DataFrame(members),'ANALYSIS_MEMBERSHIP.tsv');write(pd.DataFrame(excluded),'EXCLUSIONS.tsv');write(pd.DataFrame(index),'ANALYSIS_INDEX.tsv')
# Validate actual population cardinalities without using long-table row counts as n.
summary=pd.DataFrame(index).set_index('analysis_id')
check('regional70_native60',summary.loc['CAPSCAN_REGIONAL_PAIRED14','participant_region_cells']==70 and summary.loc['CAPSCAN_NATIVE12','participant_region_cells']==60)
check('native12_functionalDNA12_same_proteomic_inputs',Nspec==set(S.loc[S.specimen_id_crossomics.isin(Fspec)&S.participant_id.isin(nparticipants),'specimen_id_crossomics']))
check('membership_unique',not pd.DataFrame(members).duplicated(['analysis_id','entity_type','entity_id']).any())
# Explicit pairing links, never Cartesian DNA-to-acquisition matching.
write(pd.DataFrame([dict(participant_id=p,region=r,paired_stool_region='stool',source_file=str(Q/'Figure5/regional_host_tests.tsv')) for p in sorted(P14,key=int) for r in REG[:-1]]),'REGIONAL_PAIRED_LINKS.tsv')
for p in [A/'revised/derive_capscan_pc.py',A16/'frozen/search.py',B/'MM_review_analysis/A19_4_6_biology_robustness/frozen/derive_capscan_extras.py',B/'MM_review_analysis/A22_4_6_crossomics_robustness/frozen/derive_crossomics.py',DNA/'extend_current.py',N/'native_pyramids.py',MG/'run_selection.py']:source(p)
write(pd.DataFrame(SOURCES.values()),'SOURCE_MANIFEST.tsv')
(O/'VALIDATION.json').write_text(json.dumps(dict(status='PASS',checks=CHECKS,analyses=index,observations=dict(missing_DNA_cells=missing_cells,matched_public_vs_inclusive_same_acquisition_membership=True,benchmark_excluded=49,quantitative_empty=1,native_DNA_input_runs=228,active_DNA_runs=301,native_protein_metabolite_specimens=257,all_metabolite_specimens=288,valid_crossomics_metabolite_specimens=287),limits=['Local restricted review only.','No raw search, LFQ, metagenome assembly or sequencing inference rerun.','Acquisition-to-specimen links are actual design-table links; no DNA-to-MS acquisition pairing is invented.','The native/DNA layers meet at participant-region; they are not identical specimen sets.']),indent=2)+'\n')
if Path(__file__).resolve() != (O/'build_rosters.py').resolve(): shutil.copy2(Path(__file__),O/'build_rosters.py')
print(pd.DataFrame(index)[['analysis_id','acquisitions','specimens','participants','participant_region_cells','DNA_runs']].to_string(index=False))
print('PASS',len(CHECKS),'checks',len(members),'membership rows',len(excluded),'exclusions')
# Reader-facing refinement: DNA specimens/source categories are independent of MS inputs.
mem=pd.read_csv(O/'ANALYSIS_MEMBERSHIP.tsv',sep='\t',dtype=str,keep_default_na=False)
extra=[]
for (aid,spec),g in mem[mem.entity_type.eq('DNA_run')].groupby(['analysis_id','specimen_id_crossomics'],sort=False):
 rr=g.iloc[0].to_dict();rr.update(entity_type='DNA_specimen',entity_id=spec,run_accession='');extra.append(rr)
mem=pd.concat([mem,pd.DataFrame(extra)],ignore_index=True);write(mem,'ANALYSIS_MEMBERSHIP.tsv')
ix=pd.read_csv(O/'ANALYSIS_INDEX.tsv',sep='\t',keep_default_na=False);ix['proteomic_specimens']=ix.specimens;ix['DNA_specimens']=0;ix['gut_stool_participant_region_cells']=0;ix['saliva_participant_region_cells']=0
for i,r in ix.iterrows():
 z=mem[mem.analysis_id.eq(r.analysis_id)];ds=set(z.loc[z.entity_type.eq('DNA_specimen'),'specimen_id_crossomics']);ps=set(z.loc[z.entity_type.eq('specimen'),'specimen_id_crossomics']);ix.loc[i,'DNA_specimens']=len(ds);ix.loc[i,'specimens']=len(ds|ps);ix.loc[i,'gut_stool_participant_region_cells']=len(z[z.entity_type.eq('participant_region')&z.region.isin(REG)]);ix.loc[i,'saliva_participant_region_cells']=len(z[z.entity_type.eq('participant_region')&z.region.eq('saliva')]);write(z,r.analysis_id+'.tsv')
mask=ix.analysis_id.eq('CAPSCAN_DNA_COMPLETED301');ix.loc[mask,'purpose']='Upstream completed metagenomic processing runs, including saliva';ix.loc[mask,'limitations']='301 distinct DNA specimens comprise 272 gut/stool runs and 29 saliva runs. Eighty-eight source-category cells include 15 saliva cells; only selected matched gut/stool inputs enter later analyses. No proteomic acquisition membership is inferred from DNA alone.'
write(ix,'ANALYSIS_INDEX.tsv')
ex=pd.read_csv(O/'EXCLUSIONS.tsv',sep='\t',dtype=str,keep_default_na=False)
ex['disposition_type']='analysis_eligibility'
ex.loc[ex.analysis_id.eq('CAPSCAN_PUBLIC_MG_282'),'disposition_type']='comparator_unavailable'
ex.loc[ex.analysis_id.eq('CAPSCAN_DNA_COMPLETED301'),'disposition_type']='documented_primary_input_exclusion'
ex.loc[ex.analysis_id.eq('CAPSCAN_QUANT_NONEMPTY'),'disposition_type']='quantitative_zero_input'
ex.loc[ex.entity_type.eq('DNA_run')&ex.analysis_id.eq('CAPSCAN_NATIVE12'),'disposition_type']='not_in_downstream_analysis'
for i,r in ex[ex.entity_type.eq('DNA_run')&ex.analysis_id.eq('CAPSCAN_NATIVE12')].iterrows():
 if r.region=='saliva':ex.loc[i,'reason']='Saliva sequencing run: outside the gut/stool proteomic cross-omics design; retained in the upstream completed DNA roster.'
write(ex,'EXCLUSIONS.tsv')
non=[]
for r in ix.itertuples():
 ids=set(mem.loc[mem.analysis_id.eq(r.analysis_id)&mem.entity_type.eq('acquisition'),'acquisition_id'])
 if not ids:continue
 for q in base[~base.acquisition_id.isin(ids)].itertuples():
  if r.analysis_id=='CAPSCAN_PUBLIC_MG_282':reason='MG comparator absent from A16 task manifest; this is comparator availability, not failed public search.'
  elif r.analysis_id in ['CAPSCAN_QUANT_NONEMPTY','CAPSCAN_REGIONAL_DESCRIPTIVE15']:reason='Zero total quantitative intensity.'
  elif q.participant_id=='15':reason='Participant 15 has stool-only proteomics and cannot form the required five-region block.'
  elif r.analysis_id in ['CAPSCAN_NATIVE12','CAPSCAN_FUNCTIONAL_DNA12','CAPSCAN_FUNCTIONAL_NONDNA12_SENS'] and q.participant_id not in nparticipants:reason='Whole participant outside DNA-complete12 design: participant 5 lacks eligible colon DNA; participant 7 lacks eligible duodenum DNA.'
  elif q.specimen_id_crossomics not in metspec and r.analysis_id.startswith(('CAPSCAN_FUNCTIONAL','CAPSCAN_NATIVE')):reason='Specimen lacks a metabolite row or any positive microbial-OG acquisition; see specimen eligibility exclusions.'
  elif q.acquisition_id not in qids:reason='Zero total quantitative intensity.'
  else:raise AssertionError((r.analysis_id,q.acquisition_id,'Unexplained nonmembership'))
  non.append(dict(analysis_id=r.analysis_id,acquisition_id=q.acquisition_id,specimen_id_biology=q.specimen_id_biology,specimen_id_crossomics=q.specimen_id_crossomics,participant_id=q.participant_id,region=q.region,disposition_type='not_in_analysis',reason=reason))
write(pd.DataFrame(non),'ACQUISITION_NONMEMBERSHIP.tsv')
write(active.groupby('group').agg(DNA_runs=('run_accession','size'),DNA_specimens=('specimen','nunique'),participants=('subject_id','nunique')).reset_index().rename(columns={'group':'source_region_or_material'}),'DNA301_SOURCE_CATEGORIES.tsv')
check('DNA301_includes29saliva_and272gutstool',active.group.eq('saliva').sum()==29 and active.group.isin(REG).sum()==272)
check('complete_membership_unique',not mem.duplicated(['analysis_id','entity_type','entity_id']).any())
check('DNA301_specimen_count_explicit',int(ix.loc[mask,'DNA_specimens'].iloc[0])==301)
val=json.loads((O/'VALIDATION.json').read_text());val['checks']=CHECKS;val['analyses']=ix.to_dict('records');val['observations']['upstream_DNA_saliva_runs']=29;val['observations']['upstream_DNA_gut_stool_runs']=272;val['membership_rows']=len(mem);val['acquisition_nonmembership_rows']=len(non);(O/'VALIDATION.json').write_text(json.dumps(val,indent=2)+'\n')
(O/'HANDOFF.md').write_text('''# CapScan M01 roster delivery — restricted local review

Files are stable. No manuscript, central ledger, repository or public release was changed; no jobs were submitted. `build_rosters.py` reproduces the tables from the actual producing sources; `SOURCE_MANIFEST.tsv` pins those inputs. `VALIDATION.json` reports all checks.

## Exact denominators

- All searches: 331 acquisitions, 292 proteomic specimens, 15 participants.
- Same-public versus specimen-MG benchmark: 282 acquisitions, 250 biology-design specimens, 15 participants. Its original A16 specimen keys use participant/device IDs; retain `BENCHMARK282_ACQUISITIONS.tsv` and the exact-acquisition crosswalk. The old assembly-inclusive benchmark has the same acquisition membership but different peptide results.
- Nonempty quantitative inputs: 330 acquisitions, 291 specimens, 15 participants, 71 observed participant-region cells. One jejunal acquisition from participant 5/specimen main:1181 has zero total quantitative intensity.
- Paired regional/own-stool analysis: 326 acquisitions, 287 specimens, 14 participants, 70 cells. Participant 15 contributes stool only and remains in descriptive summaries.
- Non-DNA functional cross-omics: 322 input acquisitions, 283 protein/metabolite specimens, 14 participants, 70 cells. Its seven displayed comparisons are part of the 16-test functional family.
- Native cross-omics: 290 protein input acquisitions, 257 protein/metabolite specimens, 12 participants, 60 cells, plus 228 DNA runs/specimens. DNA and protein/metabolite specimens can differ within participant-region. Native has a separate six-test family.
- Functional DNA comparisons use the same 12 participants and 228 selected DNA runs; nine displayed comparisons share the 16-test functional family with the seven 14-participant non-DNA comparisons. Seven extra non-DNA comparisons on the 12-person subset are a sensitivity, marked not displayed and without q_family in the source.
- Upstream DNA processing: 301 completed runs from 301 distinct DNA specimens, including 29 saliva and 272 gut/stool runs. Do not describe all 301 as gut/stool inputs to the cross-omics figure. The source has 88 participant-material/region cells, including 15 saliva cells; these are not the 70/60 analysis-cell denominators.

Participants 5 and 7 leave the complete DNA design because eligible colon DNA (5) and duodenal DNA (7) are unavailable. The completed-run exclusion is separately documented: SRR18585253, participant 15/specimen1637, lacked the accepted mate-2 input. It is not the explanation for dropping participants 5/7.

## Output schema for the parent workbook

`ANALYSIS_INDEX.tsv` has one row per stable analysis_id, human purpose, biological unit, inclusion rule, counts, testing family, source file/hash, roster file and limits. `specimens` counts the union of explicit source specimen identifiers; `proteomic_specimens` and `DNA_specimens` retain layer-specific counts. Acquisition counts denote input records; not every feature is detected in every acquisition. `gut_stool_participant_region_cells` and `saliva_participant_region_cells` keep upstream material categories separate.

`ANALYSIS_MEMBERSHIP.tsv` is a long table. `entity_type` is acquisition, specimen (proteomic), DNA_specimen, DNA_run, participant, or participant_region. Filter by analysis_id AND entity_type before counting biological replicates. Each per-analysis TSV has exactly the same schema. `entity_id` is a display key; original acquisition, biology specimen, cross-omics specimen, run, participant and region fields remain separate. The two specimen namespaces were linked using the exact common acquisition IDs and verified as one-to-one; no prefix stripping or guessed joins was used.

`BASE_ACQUISITIONS.tsv` and `IDENTIFIER_CROSSWALK.tsv` preserve experimental identity. `BENCHMARK282_ACQUISITIONS.tsv` retains original A16 keys and counts. `DNA301_RUN_SPECIMEN_ROSTER.tsv` and `DNA301_SOURCE_CATEGORIES.tsv` expose upstream DNA membership/materials. `DNA_CELL_ELIGIBILITY.tsv` gives explicit 14-by-five region availability in all three functional representations. `REGIONAL_PAIRED_LINKS.tsv` lists the 56 participant-specific capsule-to-stool links; it does not imply acquisition pairing.

`EXCLUSIONS.tsv` records stage-specific eligibility, comparator absence, zero quantitative input and the documented missing-mate exclusion. `disposition_type` distinguishes them. `ACQUISITION_NONMEMBERSHIP.tsv` gives every acquisition outside each MS analysis roster relative to the 331-search design; nonmembership is not automatically a processing failure or a causal exclusion. No DNA-to-acquisition Cartesian join is created.

## Evidence and limits

Producer selection rules are in A33/revised/derive_capscan_pc.py:24–42 and139–145; A16/frozen/search.py:40–60; A19/frozen/derive_capscan_extras.py:94–115; A22/frozen/derive_crossomics.py:61–97; capscan_dna/extend_current.py:21–89; native_pyramids.py:47–78 and204–233; and capscan_metag_analysis/run_selection.py:18–40. Exact full paths and hashes are in SOURCE_MANIFEST.tsv.

This is a fresh roster/eligibility reconstruction from retained producing sources. It does not reidentify spectra, re-estimate LFQ, assemble DNA, rerun metabolite preprocessing or validate experimental causality. Keep the local restricted-review designation; root owns the reader workbook and final Methods locator.
''')
if Path(__file__).resolve() != (O/'build_rosters.py').resolve(): shutil.copy2(Path(__file__),O/'build_rosters.py')
print('FINAL PASS',len(CHECKS),'checks;',len(mem),'membership rows;',len(non),'nonmembership rows')
# Final independent crosswalk guards and stable output receipt.
j=active.merge(S[['specimen_id_crossomics','participant_id','region']],left_on='specimen',right_on='specimen_id_crossomics',validate='1:1')
check('all_explicit_DNA_to_proteomic_specimen_links_metadata_match',(j.subject_id==j.participant_id).all() and (j.group==j.region).all())
a16full=read(A16/'results/all_acquisitions.tsv');a16full=a16full[(a16full.cohort=='CAPSCAN')&(a16full.test_arm=='FL_public')&(a16full.reference_arm=='MG')]
check('public_MG5_and9_same_acquisition_roster',set(a16full.loc[a16full.min_length==5,'sample'])==set(a16full.loc[a16full.min_length==9,'sample']))
p=O/'VALIDATION.json';v=json.loads(p.read_text());v['checks']=CHECKS;v['observations']['upstream_DNA_specimens_with_explicit_proteomic_design_link']=len(j);p.write_text(json.dumps(v,indent=2)+'\n')
if Path(__file__).resolve() != (O/'build_rosters.py').resolve(): shutil.copy2(Path(__file__),O/'build_rosters.py')
outputs=[dict(path=str(p),sha256=sha(p),bytes=p.stat().st_size) for p in sorted(O.iterdir()) if p.is_file() and p.name!='OUTPUT_MANIFEST.json']
(O/'OUTPUT_MANIFEST.json').write_text(json.dumps(dict(status='STABLE_RESTRICTED_REVIEW',files=outputs),indent=2)+'\n')
