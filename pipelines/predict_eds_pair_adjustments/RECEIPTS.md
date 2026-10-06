# RECEIPTS

## Original Slurm jobs

| role | job |
|---|---|
| baseline | 6044262 |
| refit | 6048399 |

## Code

| original (LAB_ROOT-relative) | file here | sha256 original | sha256 here |
|---|---|---|---|
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/separate_adjustments/run_analysis.py` | `run_analysis.py` | `e27b62da1d938229a6cea7688d7ecbf335417fe0c3b62c16b533bb08832597a3` | `fb80effee67693c25749053816b125bbb0e43c0f184d6ed8295a15f73d8f98e3` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/separate_adjustments/TARGETS.json` | `TARGETS.json` | `6270ac79a4ae241989a81df47e920316d3dfcb95b2925779a8b135de2d380a5d` | `6270ac79a4ae241989a81df47e920316d3dfcb95b2925779a8b135de2d380a5d` |

## Inputs (sha256 recorded by the original run)

| input (LAB_ROOT-relative) | sha256 |
|---|---|
| `analyses/enterotypes_2026-09-07/PREDICT_metabolite_annotations.tsv` | `14e044fbbf5f6dedcd68231a6b25053f38ec4e4b6fbc2d7ad953c38b580f0121` |
| `analyses/enterotypes_2026-09-07/PREDICT_metabolite_intensity.tsv` | `0663b71545e31e677368c1fcbfe7d2c503e987f4e5a4eff677da6d8848418cc9` |
| `denovo_vs_tpm/metab_reverse_2026-09-14/reference/compound_reaction.tsv` | `09084463946eeb27bcccb6871c489248935d676cd7d1d66772e0f10d9d2e5689` |
| `denovo_vs_tpm/metab_reverse_2026-09-14/reference/reaction_ko.tsv` | `539b237200129db56c91ab5308f0e76d7ea58500d4c098b415808ad6cf00a2e1` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/analyze.py` | `9925977db928a866b73ef424b34de45a99d04cc571c160b83f099e5d496b1a1a` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/input_data/KO_profile.tsv.gz` | `164fe6cde20cec1e6cc6803c7b2ad50a3f373f02e411e211769670b0f39cec02` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/input_data/Pathway_profile.tsv.gz` | `8453f2b543726afcdb30a7b29159c7cf02f009a7598dd72e9f777dd2396c40f7` |
| controlled input `heatmap_85_patient_order.tsv` (not in repository) | `43dc72c4b31e4bcf44c2c6bfbcc16a0a713ea47a8f377e86c87280fc5327b9ba` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/input_data/heatmap_85_values.npz` | `11c171206d200b9051c0f785041288bbc621be61695baae1bfa8eaffb88b7b37` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/input_data/heatmap_features.tsv` | `b883e1659743f5437c3d1d3981abd07a36d1ae41ee5d87b6321596fe90968e51` |
| `reviews/ipath3_v6.2_MM_2026-09-13/source_data/Simple_PREDICT_metabolite_id_audit.tsv` | `5f6875b82b2b4baeaa0fa8ab993e026df4216a29c0ed12adca3dafb8e9aec927` |
| `simplified_grouping/final_revision/Fig5/kegg_cache/compound_to_pathway.tsv` | `692cc003e36f7753beae7fa45f33c4ea9cc28d2b4f540bc78a2d118d4500bfb0` |

## Outputs (sha256 recorded by the original run)

| output (LAB_ROOT-relative) | sha256 |
|---|---|

## Path aliases

Lab path segments are shortened in code and tables (`jess_integration_2026-09-29` to `integration_2026-09-29`, `v10_mm_philip_comment_audit_2026-09-22` to `v10_comment_audit_2026-09-22`, `coauthor_comment_check_v25_2026-10-05` to `comment_check_v25_2026-10-05`). `FASTALAKE_DATA_DIR` must hold the frozen inputs under the shortened names; the original names are the lab directory names.
