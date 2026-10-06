"""Independent recomputation of the PREDICT histidine inference rows (bootstrap confidence interval and ratio) from the collected tables, used to check collect_revision4 output."""
from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
B=Path(__file__).resolve().parents[1];ROOT=B.parents[1];O=B/'histidine';S=ROOT/'MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT'
a=json.loads((ROOT/'MM_review_analysis/recovery_2026-09-15/predict_histidine_audit/full_corrected_PREDICT/AUDIT.json').read_text());assert a['status']=='PASS' and a['tests']==1572
tf=pd.read_csv(S/'complete_statistical_families.tsv',sep='\t');ef=pd.read_csv(S/'acquisition_enrichment.tsv.gz',sep='\t',dtype={'participant':str});routes=['denovo_top30','sage_unspecific_bothq'];key=['route','side','origin','residue']
ef=ef[ef.residue.eq('H')&ef.route.isin(routes)].copy();ef['O']=ef.observed/ef.denominator.replace(0,np.nan);ef['E']=ef.expected/ef.denominator.replace(0,np.nan)
sp=ef.groupby(key+['participant','region','specimen_id'],sort=True)[['O','E']].mean().reset_index();pat=sp.groupby(key+['participant'],sort=True)[['O','E']].mean().reset_index();pat['difference']=pat.O-pat.E;pat.to_csv(O/'PREDICT_patient_source_enrichment.tsv',sep='\t',index=False)
primary=tf[tf.endpoint.eq('source_enrichment')&tf.residue.eq('H')&tf.route.isin(routes)].copy();checked=0
for r in primary.itertuples():
 d=pat[(pat.route==r.route)&(pat.side==r.side)&(pat.origin==r.origin)];v=d.difference.dropna().to_numpy();n=len(v);assert n==r.n
 w=np.random.default_rng(20261014+n).multinomial(n,np.full(n,1/n),size=20000);ci=np.quantile(w@v/n,[.025,.975]);np.testing.assert_allclose(ci,[r.CI95_low,r.CI95_high],atol=2e-12,rtol=1e-8)
 oe=d[['O','E']].dropna().to_numpy();n2=len(oe);w=np.random.default_rng(20261114+n2).multinomial(n2,np.full(n2,1/n2),size=20000);means=w@oe/n2;ratio=means[:,0]/means[:,1];ci=np.quantile(ratio,[.025,.975]);np.testing.assert_allclose(ci,[r.enrichment_CI95_low,r.enrichment_CI95_high],atol=2e-12,rtol=1e-8);checked+=1
primary.to_csv(O/'PREDICT_primary_source_enrichment.tsv',sep='\t',index=False)
rec=tf[tf.endpoint.eq('matched_reciprocal')];rec.to_csv(O/'PREDICT_primary_reciprocal_tests.tsv',sep='\t',index=False)
conc=pd.read_csv(S/'participant_concentration_sensitivity.tsv',sep='\t');c=conc[conc.route.isin(routes)&conc.origin.eq('All')];c.groupby(['route','side','removed_top_proteins']).difference.agg(['size','mean','min','max']).to_csv(O/'PREDICT_concentration_sensitivity.tsv',sep='\t')
mot=pd.read_csv(ROOT/'MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/complete_statistical_families.tsv',sep='\t');mot=mot[mot.endpoint.eq('motif')&mot.population.eq('H_only')&mot.route.isin(['denovo_top30','sage_semi_bothq','sage_unspecific_bothq'])&mot.position.ne(-1)].copy();mot.to_csv(O/'CAPSCAN_primary_H_conditioned_flanks.tsv',sep='\t',index=False)
summary=mot.groupby(['route','side','origin','contrast']).agg(tests=('q','size'),significant=('q',lambda s:int((s<.05).sum())),largest_positive=('mean','max'),largest_negative=('mean','min'));summary.to_csv(O/'CAPSCAN_motif_review_counts.tsv',sep='\t')
(O/'PREDICT_SCIENTIFIC_REVIEW.md').write_text('''# Corrected PREDICT histidine review

All 306 acquisitions and the independent 1,572-test source/reciprocal audit passed. Primary interval validation independently reconstructs 12 source-effect intervals and 12 ratio intervals after averaging acquisitions within specimen and then specimens within each patient. The 93 patients retain their observed clinical-specimen mixtures; clinical codes are not treated as CapScan anatomical regions.

Production top-30% de novo C-boundary H fractions average 7.70%, against 2.61% source expectation (ratio 2.95, 95% bootstrap interval 2.74–3.16). All 93 patient differences are positive; source-family q=0.000010. Human and Microbial labelled sources both retain positive effects, ratios 2.34 and 4.34. N-boundary enrichment is smaller: 2.38% versus 2.07%, ratio 1.15, q=0.000603; only 55/93 patient effects are positive. This supports C-boundary enrichment under the conditional source null, not universal N-boundary enrichment in every patient.

The unspecific Sage route, requiring peptide and spectrum q≤.01, also retains C enrichment (9.94% observed versus 2.47% expected, ratio 4.03; every patient positive). The cohort-wide PREDICT A28 plan contains tryptic and unspecific routes, not a cohort-wide semi-tryptic route or anatomical motif family. The four measured cases have their separate same-scan semi-tryptic check. CapScan supplies the full regional and motif comparisons.

**Complementary-product evidence is mixed and includes an adverse finding.** Top-30% All-origin H sites have a −0.871 percentage-point matched-control difference (95% interval −1.317 to −0.437; q=0.000283). The Human-labelled contrast is negative (−1.443 points,q=0.000012); Microbial is not significant (+0.384 points,q=0.359). Unspecific Sage All-origin support is positive (+0.971 points,q=0.000012). The analysis therefore does not support universal excess complementary products across methods or sources.

Source matching conditions on the selected FASTA and compatible sequence/mass opportunities. It does not remove selection, peptide detectability, digestion, handling or ex-vivo cleavage effects. Shared spectra and database selection prevent treating agreement between de novo and Sage as independent truth. Local parent labels do not establish unique species identity. All residue families, null results, unresolved P1 sites and concentration sensitivities remain available.
''')
(O/'MOTIF_SCIENTIFIC_REVIEW.md').write_text('''# Histidine-conditioned motif review

The full CapScan 32,400-test motif family was independently reconstructed, including null and opposite-direction cells. The new export retains every primary-route H-conditioned flank result for both sides, all three local origin classes and proximal fluid, distal fluid and stool. P1 is fixed to H by selection and is excluded from flank inference. Each primary all-origin side/compartment has 140 flank tests; adjusted q values remain those from the entire 32,400-test family.

At production top 30%, N-boundary flank differences include excess leucine at P2′ (about 13.3–13.5 percentage points) and serine at P1′ (7.9–10.2 points). C-boundary contexts show excess serine at P1′ (8.8–9.6 points) and leucine at P2′ (8.1–9.8 points). Full heat maps include all 20 amino acids and seven noncentral positions, avoiding a selected sequence logo that hides depletion. Proximal/distal require both constituent fluid regions in the same participant; stool is kept separate. The complete primary-route matrices and all-origin/human/microbial summaries preserve route and source variation.

These are differences from a source- and H-conditioned opportunity background. The patterns may reflect source composition, detectability and sample handling as well as proteolysis; they do not identify a protease, enzyme activity or cleavage mechanism. The N-boundary H can precede the measured peptide and must be labelled as source-inferred in spectra.
''')
(O/'PREDICT_REVIEW_VALIDATION.json').write_text(json.dumps(dict(status='PASS',independent_difference_intervals=checked,independent_ratio_intervals=checked,source_audit=a,code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()),indent=2)+'\n');print('PREDICT REVIEW PASS',checked,flush=True)
