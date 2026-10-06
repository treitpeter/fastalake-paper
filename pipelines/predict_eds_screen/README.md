# predict_eds_screen

Internal ID: A70 (previous folder name A70_eds_biochemical_screen)

Which of 20 fixed display features, all measured metabolites, and all KO and KEGG pathway profiles of the PREDICT non-host protein signal correlate with the dysbiosis score (EDS) across 85 patients, and which metabolite-enzyme pairs linked by KEGG reactions are inverse? Used by Fig4 h, SFig10 a, STab4 and the EDS text of the supplementary methods.

Status: code copied as it ran (path edits only), not re-run in this repository. See RECEIPTS.md for original job ids and hashes.

## Inputs

| input (LAB_ROOT-relative) | sha256 |
|---|---|
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/inputs/matrices_repaired/PREDICT_KO_microbial_broad.tsv.gz` | `7184e3da856ef7e3a3eda180d21691dc4a8615f06fb5dba11df3c64238921474` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/inputs/matrices_repaired/PREDICT_Pathway_microbial_broad.tsv.gz` | `9b7a89360209f0803ed84693679873533ac13f4b40d6192b690fccbc2d9dbf26` |
| controlled input `heatmap_85_patient_order.tsv` (not in repository) | `4ee42bf87fdc5b73878ff70d20c2f0d995bacd8e35ea378baf5738a45256384d` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/PREDICT_repaired__original/source_data/Figure4/heatmap_85_values.npz` | `0cf55ba3f6ce26b4051ad33fb7b85fb066ed37c3a0291cb3282f8abfbf9be886` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/PREDICT_repaired__original/source_data/Figure4/heatmap_columns.tsv` | `60a773ff51861574930a4a17da3c55143f2c02a05872010c0eed2521c62490ad` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/PREDICT_repaired__original/source_data/Figure4/heatmap_features.tsv` | `a5acdfec0e270dc04097699c610d01dfae518a33bbb09e932138fd96588d9eb1` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/PREDICT_repaired__original/source_data/Figure4/heatmap_values.npz` | `862037a97826a3e729bb52d37e37d675fe9024e61d84c39c9e5ec269e0d097b5` |
| `analyses/enterotypes_2026-09-07/PREDICT_metabolite_annotations.tsv` | `14e044fbbf5f6dedcd68231a6b25053f38ec4e4b6fbc2d7ad953c38b580f0121` |
| `analyses/enterotypes_2026-09-07/PREDICT_metabolite_intensity.tsv` | `0663b71545e31e677368c1fcbfe7d2c503e987f4e5a4eff677da6d8848418cc9` |
| `denovo_vs_tpm/metab_reverse_2026-09-14/reference/compound_reaction.tsv` | `09084463946eeb27bcccb6871c489248935d676cd7d1d66772e0f10d9d2e5689` |
| `denovo_vs_tpm/metab_reverse_2026-09-14/reference/reaction_ko.tsv` | `539b237200129db56c91ab5308f0e76d7ea58500d4c098b415808ad6cf00a2e1` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/ADDITIONAL_LINKS.json` | `2bae424ab02b718a7a2b592928214be7ca103d7a373fce88c611519ca766374d` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/INPUTS.json` | `1d7d15d2a79bddedb6cacafd8d66c929b253c8fc049f1b805ab2d06a5d45c337` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/all_EDS_associations.tsv` | `a46342612d7f4bf90c538f39cd4854d641c8c8fe1db55f223c48a82cd1b1a4d2` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/all_linked_associations.tsv` | `33264126a194d835877b8c095ce85c9c00b6eeddc159a1e3152919f4c4095f00` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/analyze.py` | `9925977db928a866b73ef424b34de45a99d04cc571c160b83f099e5d496b1a1a` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/augmented_pair_screen.tsv` | `8c81ed1bded410d86ddc450c8e7e1833eed334a20d5d1e81bec27e6ba71d2f30` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/individual_followup.py` | `86ecade6b17e9291a38b57b72cf671af01a7dfa8fb627d5f2da0f52eb85fadfb` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/separate_adjustments/TARGETS.json` | `6270ac79a4ae241989a81df47e920316d3dfcb95b2925779a8b135de2d380a5d` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/separate_adjustments/all_models.tsv` | `8991f713a7597c77ff7cdf2faa7b764630e6355a3c1902b350b2d673d2ad461e` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/separate_adjustments/run_analysis.py` | `e27b62da1d938229a6cea7688d7ecbf335417fe0c3b62c16b533bb08832597a3` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/input_data/KO_profile.tsv.gz` | `164fe6cde20cec1e6cc6803c7b2ad50a3f373f02e411e211769670b0f39cec02` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/input_data/Pathway_profile.tsv.gz` | `8453f2b543726afcdb30a7b29159c7cf02f009a7598dd72e9f777dd2396c40f7` |
| controlled input `heatmap_85_patient_order.tsv` (not in repository) | `43dc72c4b31e4bcf44c2c6bfbcc16a0a713ea47a8f377e86c87280fc5327b9ba` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/input_data/heatmap_85_values.npz` | `11c171206d200b9051c0f785041288bbc621be61695baae1bfa8eaffb88b7b37` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/input_data/heatmap_features.tsv` | `b883e1659743f5437c3d1d3981abd07a36d1ae41ee5d87b6321596fe90968e51` |
| `reviews/ipath3_v6.2_MM_2026-09-13/source_data/Simple_PREDICT_metabolite_id_audit.tsv` | `5f6875b82b2b4baeaa0fa8ab993e026df4216a29c0ed12adca3dafb8e9aec927` |
| `simplified_grouping/final_revision/Fig5/kegg_cache/compound_to_pathway.tsv` | `692cc003e36f7753beae7fa45f33c4ea9cc28d2b4f540bc78a2d118d4500bfb0` |

## Outputs

`eds_refit/all_EDS_associations.tsv`, `eds_refit/augmented_pair_screen.tsv`, `eds_refit/input_data/heatmap_85_values.npz` (85-patient values and within-feature z scores), `heatmap_features.tsv`, `heatmap_85_patient_order.tsv`, `all_linked_associations.tsv`, `joined_evidence.tsv` (hashes in RECEIPTS.md).

## Slurm resources and run time (original logs, `sacct`)

| step | job | cores | requested mem | peak RSS | wall |
|---|---|---|---|---|---|
| baseline replay of the earlier screen (`eds_biochemistry.py baseline`) | 6044262 | 2 | 8G | 0.27 GB | 1 min |
| refit (`eds_biochemistry.py refit`, one of five steps of the downstream job) | 6048399 | 4 | 16G | 1.22 GB for the whole job | 4 min 9 s for the whole job |

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| minimum paired observations | 43 positives; at least 3 distinct values per variable | the original screen's positive-observation threshold, kept | none tried |
| non-positive values | treated as missing, no imputation | zeros are absence of signal, not measurements | none |
| statistic | Spearman rho; adjusted model uses rank residuals after least squares on covariates | same rule as the earlier screen | unadjusted model is reported beside it |
| covariates | rank host fraction, rank log2 protein groups, antibiotic state, clinical group dummies | host depth and clinical state confound EDS | 16 subsets in A78 |
| multiplicity | BH within layer x normalization x model | one family per displayed table | none |
| normalizations | native intensity; annotated-relative (row sum) for KO and pathway | host depth sensitivity | none |
| features | the same 20 fixed display features, no reselection after the profile repair | avoids selecting on the new result | none |
| bootstrap in `individual_followup.py` | seed 2026093031, 2000 resamples, Spearman for K01667 against two metabolites | literature-guided follow-up, stated post hoc | none |
| reaction links | all KEGG compound-reaction-KO links kept, including targets absent from the data; no direction filter | no favourable-direction selection | none |

## Notes

Code edits relative to the original: lab absolute paths replaced by `FASTALAKE_DATA_DIR`; lab-tree roots located by a marker file replaced by `FASTALAKE_DATA_DIR`. No statistics changed. The scripts keep their original integrity gates (`verify()`, `gate()`, hash checks against `*_CODE.json` and `*_INPUTS.json` manifests of the frozen run). Those manifests contain absolute lab paths and are not in the repository; the hashes they held are listed in RECEIPTS.md. Because the code here carries the path edits described below, its own sha256 differs from the original (both are in RECEIPTS.md), so a re-run needs the manifests regenerated for this copy.The earlier driver versions in the original `biology_refit/` folder produced the pre-repair screen and are not used for the published outputs. `analyze.py` and `individual_followup.py` are the copies stored beside the refit outputs; the driver itself copies them from the original screen folder at run time. The refit mode also needs `groups/` (supported-human and other-unresolved group matrices, written by the grouping step `build_groups.py`, which is outside this pipeline group), `$OUT/C/profiles/PREDICT_noclust` (A88) and `$OUT/C/clinical_refit/acquisition_accounting.tsv` (A91, controlled). `common.py` is shared with A84, A87 and A108, which add this folder to `PYTHONPATH`. The separate-adjustment models are in A78.
