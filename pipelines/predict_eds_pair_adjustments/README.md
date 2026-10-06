# predict_eds_pair_adjustments

Internal ID: A78 (previous folder name A78_eds_five_pair_adjustments)

For five post hoc metabolite-enzyme pairs, how does the partial rank correlation change when host fraction, protein-group depth, antibiotic state and clinical group are added separately and together (all 16 covariate subsets)? Used by STab4 and the counts at v28d line 799.

Status: code copied as it ran (path edits only), not re-run in this repository. See RECEIPTS.md for original job ids and hashes.

## Inputs

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

Also the frozen `FREEZE.json` and `TARGETS.json` of the original run (`TARGETS.json` is shipped here).

## Outputs

`eds_refit/separate_adjustments/all_models.tsv` (saved model table written by `run_analysis.py`; no file named `five_pair_all_models.tsv` was found, see OPEN_ITEMS.md), `category_counts.tsv`, `host_depth_correlation.tsv`, `patient_rosters.tsv` (hashes in RECEIPTS.md).

## Slurm resources and run time (original logs, `sacct`)

Run inside the EDS drivers: baseline job 6044262 (2 cores, 1 min wall, 0.27 GB peak RSS, whole job); refit job 6048399 (4 cores, 4 min 9 s wall, 1.22 GB peak RSS, whole job).

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| covariate subsets | all 16 subsets of host, depth, antibiotics, clinical | each adjustment shown separately, no model preferred by significance | none |
| statistic | Spearman-type correlation of rank residuals, df = n - rank - 1, two-sided t | same as the main screen | none |
| observations | positive finite values of both variables | zeros are missing | none |
| patients | 85 fixed | fixed EDS roster | none |
| pairs | `TARGETS.json` (5 pairs) | post hoc pairs approved after the screen | pairs are stated as post hoc |

## Notes

Code edits relative to the original: lab absolute paths replaced by `FASTALAKE_DATA_DIR`; lab-tree roots located by a marker file replaced by `FASTALAKE_DATA_DIR`. No statistics changed. The scripts keep their original integrity gates (`verify()`, `gate()`, hash checks against `*_CODE.json` and `*_INPUTS.json` manifests of the frozen run). Those manifests contain absolute lab paths and are not in the repository; the hashes they held are listed in RECEIPTS.md. Because the code here carries the path edits described below, its own sha256 differs from the original (both are in RECEIPTS.md), so a re-run needs the manifests regenerated for this copy.Executed by `eds_biochemistry.py` of A70, which copies `run_analysis.py` and `TARGETS.json` into `eds_*/separate_adjustments/` and runs `pilot` then `full`. Needs `A70` outputs (`INPUTS.json`, `augmented_pair_screen.tsv`) one directory up.
