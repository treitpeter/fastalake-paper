"""sample_rosters_workbook: build histidine scope. Moved from the original analysis as it ran, apart from path scrubbing; see README.md."""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
B=Path('LAB_ROOT')
O=B/'00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M01_sample_rosters/capscan'
A20=B/'MM_review_analysis/A20_4_7_terminal_routes';A27=B/'MM_review_analysis/A27_4_7_histidine_defence';A28=B/'MM_review_analysis/A28_4_7_histidine_inference';F=A28/'results_revision4/full_CAPSCAN';src=[];checks={}
def reg(p):
 p=Path(p)
 if str(p) not in [x['source_file'] for x in src]:src.append(dict(source_file=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest(),bytes=p.stat().st_size))
 return str(p)
def rd(p,**kw):reg(p);return pd.read_csv(p,sep='\t',**kw)
def jj(p):reg(p);return json.loads(p.read_text())
def ck(k,v):checks[k]=bool(v);assert v,k
base=rd(O/'BASE_ACQUISITIONS.tsv',dtype=str);a20=rd(A20/'results/summary_CAPSCAN/acquisition_design.tsv',dtype=str)
x=a20.rename(columns={'sample':'acquisition_id','participant':'histidine_participant','region':'histidine_region','specimen_id':'histidine_specimen_id'}).merge(base,on='acquisition_id',validate='1:1')
ck('A20_exact331_292_15',len(x)==331 and x.histidine_specimen_id.nunique()==292 and x.histidine_participant.nunique()==15)
ck('A20_exact_design_metadata_match',(x.histidine_participant==x.participant_id).all() and (x.histidine_region==x.region).all() and (x.histidine_specimen_id==x.specimen_id_biology).all())
for name,p in [('A27',A27/'inputs/manifest.json'),('A28',A28/'inputs/manifest_revision4.json')]:
 m=jj(p);tasks=[t for t in m['tasks'] if t['cohort']=='CAPSCAN'];ck(name+'_same331_samples',set(t['sample'] for t in tasks)==set(x.acquisition_id))
 lookup=x.set_index('acquisition_id')
 ck(name+'_design_matches',all(str(t['design']['participant'])==lookup.loc[t['sample'],'participant_id'] and t['design']['region']==lookup.loc[t['sample'],'region'] and t['design']['specimen_id']==lookup.loc[t['sample'],'specimen_id_biology'] for t in tasks))
done=jj(F/'COMPLETE.json');ck('A28_receipt331',done['status']=='PASS' and done['acquisitions']==331)
q=rd(O/'CAPSCAN_QUANT_NONEMPTY.tsv',dtype=str);qids=set(q.loc[q.entity_type.eq('acquisition'),'acquisition_id']);x['in_quantitative330_roster']=x.acquisition_id.isin(qids)
x['histidine_scope']='Input acquisition universe; endpoint-specific finite eligibility differs.';x['source_file']=str(A20/'results/summary_CAPSCAN/acquisition_design.tsv')
x[['acquisition_id','histidine_specimen_id','specimen_id_crossomics','participant_id','region','aliquot_x','in_quantitative330_roster','histidine_scope','source_file']].rename(columns={'aliquot_x':'aliquot'}).to_csv(O/'HISTIDINE_ACQUISITION_ROSTER.tsv',sep='\t',index=False)
empty=x.loc[~x.in_quantitative330_roster,'acquisition_id'].tolist();ck('one_quantitative_exception',len(empty)==1)
z=rd(A20/'results/CAPSCAN'/empty[0]/'terminal_counts.tsv');den=z.groupby('route').denominator.max().to_dict();ck('quant_zero_has_histidine_evidence',all(den[r]>0 for r in ['denovo','semi','unspecific']))
t=rd(F/'complete_statistical_families.tsv');summary=t.groupby(['family','endpoint'],dropna=False).agg(saved_tests=('n','size'),minimum_saved_n=('n','min'),maximum_saved_n=('n','max'),saved_estimable_tests=('estimable','sum')).reset_index();summary['scope']='Saved test n, not independently replayed endpoint eligibility or statistical tests';summary['source_file']=str(F/'complete_statistical_families.tsv');summary.to_csv(O/'HISTIDINE_SAVED_FINITE_N_SUMMARY.tsv',sep='\t',index=False)
# Compact exact identities from existing participant tables, limited to H and reciprocal endpoints.
rows=[]
for p,kind,keys,metric in [(F/'participant_region_enrichment.tsv','H_source_difference',['route','side','origin','region'],'difference'),(F/'participant_reciprocal_controls.tsv','complementary_product_difference',['route','origin','region'],'difference')]:
 d=rd(p,dtype={'participant':str})
 if 'residue' in d:d=d[d.residue.eq('H')]
 for labels,g in d.groupby(keys,sort=False):
  finite=np.isfinite(pd.to_numeric(g[metric],errors='coerce'));ids=sorted(g.loc[finite,'participant'].unique(),key=int);row=dict(endpoint_kind=kind,**dict(zip(keys,labels)),finite_participants=len(ids),participant_ids=';'.join(ids),source_rows=len(g),source_file=str(p));rows.append(row)
compact=pd.DataFrame(rows);compact.to_csv(O/'HISTIDINE_FINITE_PARTICIPANT_SETS.tsv',sep='\t',index=False)
ck('all_finite_participants_in_15',all(set(r['participant_ids'].split(';'))<=set(str(x) for x in range(1,16)) for r in rows if r['participant_ids']))
ck('fluid_up_to14_stool_up_to15',compact.loc[compact.region.ne('stool'),'finite_participants'].max()<=14 and compact.loc[compact.region.eq('stool'),'finite_participants'].max()==15)
for p in [A20/'frozen/collector.py',A28/'frozen/collect_revision4.py',A28/'frozen/histidine_helpers.py']:reg(p)
pd.DataFrame(src).to_csv(O/'HISTIDINE_SOURCE_MANIFEST.tsv',sep='\t',index=False)
(O/'HISTIDINE_VALIDATION.json').write_text(json.dumps(dict(status='PASS_BOUNDED_ROSTER_SCOPE',checks=checks,acquisitions=331,specimens=292,participants=15,quantitative_excluded_acquisition=empty[0],histidine_route_max_denominators_for_that_acquisition=den,compact_endpoint_groups=len(compact),family_saved_n_ranges=t.groupby('family').n.agg(['min','max']).to_dict('index'),limits=['No raw scientific analysis rerun','Full 44,010-test eligibility membership was not reconstructed; remains endpoint/statistical audit scope','Compact participant sets cover H source differences and reciprocal controls; spatial/motif test membership remains in producing tables/code','Stable ANALYSIS_INDEX.tsv and ANALYSIS_MEMBERSHIP.tsv were not changed']),indent=2)+'\n')
(O/'HISTIDINE_ROSTER_SCOPE.md').write_text('''# CapScan histidine roster scope — M01 addendum

The existing 14-participant quantitative/regional roster is not the complete histidine acquisition universe. A20, A27 and the actual A28 revision-4 manifest all contain the same 331 acquisitions, 292 specimens and 15 participants as the search design. Their specimen/participant/region identifiers agree exactly. A28 full_CAPSCAN completion records all331 acquisitions. The 14-participant regional set gives the potential capsule–stool paired participants only; it must not supply every histidine denominator.

The zero-LFQ acquisition `20210513_EXPL0_PeTr_SA_T760_1341_20210528032734` (participant5, jejunum, specimen main:1181) is excluded from the quantitative330 roster but has actual A20 terminal evidence: maximum recorded route denominators are1264 for de novo,443 for semi-tryptic and301 for nonspecific. Thus copying the quantitative acquisition/specimen exclusions into histidine would be incorrect.

A20 descriptive participant-region summaries include15 participants in stool and up to14 in each fluid region. Fractions with zero denominators are missing, not zero (`A20/frozen/collector.py:28–47`). Paired contrasts merge the available capsule/stool endpoint cells (`:49–63`). A28 source differences are aggregated acquisition→specimen→participant-region (`A28/frozen/collect_revision4.py:182–205`). Proximal/distal summaries require both component regions, and ordinal slopes require all four fluid regions (`:208–219`). Its statistical helper drops nonfinite differences before counting n and marks n<8 non-estimable (`histidine_helpers.py:81–94`).

Saved full-family n ranges are: spatial H0–14; source enrichment0–15; motifs14–15; reciprocal controls14–15. These ranges come from the saved complete_statistical_families.tsv and do not mean every planned contrast is estimable. They are not a fresh replay of all44,010 tests. The compact finite-participant table directly lists observed participant IDs for H source-difference cells and reciprocal controls only; missing planned source cells and spatial/motif eligibility retain their original source locators.

## Files and schemas

- `HISTIDINE_ACQUISITION_ROSTER.tsv`: exact331 rows; acquisition ID, histidine specimen ID, cross-omics specimen ID, participant, region, aliquot, quantitative330 membership flag and source locator. This is the input universe, not an endpoint-specific finite-data claim.
- `HISTIDINE_SAVED_FINITE_N_SUMMARY.tsv`: one row per saved family/endpoint; test count, minimum/maximum saved n, saved estimable count and source locator. No scientific test was rerun.
- `HISTIDINE_FINITE_PARTICIPANT_SETS.tsv`: compact exact participant lists and finite counts by route/side/origin/region for H source difference, and route/origin/region for complementary-product difference. These are grouped directly from retained participant tables; they do not substitute for spatial or motif test membership.
- `HISTIDINE_SOURCE_MANIFEST.tsv`: exact source paths, SHA-256 and bytes, including producer code and actual revision-4 manifest.
- `HISTIDINE_VALIDATION.json`: bounded crosswalk and finite-count checks, with limits.

Exact sources: `MM_review_analysis/A20_4_7_terminal_routes/results/summary_CAPSCAN/acquisition_design.tsv`; `MM_review_analysis/A27_4_7_histidine_defence/inputs/manifest.json`; `MM_review_analysis/A28_4_7_histidine_inference/inputs/manifest_revision4.json`; and its `results_revision4/full_CAPSCAN/{complete_statistical_families.tsv,participant_region_enrichment.tsv,participant_reciprocal_controls.tsv}`. Detailed acquisition/specimen denominators remain in A20 per_acquisition_residues.tsv.gz/per_specimen_residues.tsv.gz and A28 acquisition_enrichment.tsv.gz. Motif participant membership remains in A28 participant_motif_differences.tsv.gz. Do not duplicate the full opportunity tables to explain cohort membership.

M01 can therefore index CAPSCAN_HISTIDINE_ROUTES as331 acquisitions/292 specimens/15 participants, with endpoint-specific finite n and the paired14-person potential universe stated separately. Full endpoint/statistical reconciliation remains M22; this addendum does not close it. Stable CapScan ANALYSIS_INDEX and ANALYSIS_MEMBERSHIP files are unchanged. Local restricted review only.
''')
print(json.dumps({'status':'PASS','checks':len(checks),'compact_endpoint_groups':len(compact),'quantitative_exception':empty[0]},indent=2))
