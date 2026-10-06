# predict_eds_focused_pairs

Internal ID: A71 (previous folder name A71_eds_focused_pairs)

Which saved screen rows and covariate-subset estimates belong to the three approved post hoc metabolite-enzyme pairs shown in SFig10 a and STab4 (Focused_* sheets)? The script reads the saved full-family tables of A70 and recomputes the 16-subset adjustment decomposition; it does not rerun the discovery screen.

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

## Outputs

`EDS_approved/source_data/`: `focused_features.tsv`, `focused_values.npz`, `focused_raw_measurements.tsv`, `pair_patient_rosters.tsv`, `selected_pair_summary.tsv`, `separate_adjustments_all_models.tsv`, `*_full_family_results.tsv`, and `patient_order_and_covariates.tsv` (carries clinical covariates, controlled; hash only). The plan names `focused_values.tsv` and `patients.tsv`; the files found have the names listed here.

## Slurm resources and run time (original logs, `sacct`)

Original job id and resources not recorded in the logs found.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| pairs | glucuronic acid with K01812, tryptophan with K01667, glycochenodeoxycholic acid with K01442 | approved post hoc examples | all pairs are in A70 |
| q values | taken from the saved full-family tables | the decomposition is descriptive | none |
| covariate subsets | 16 subsets of host, depth, antibiotics, clinical | as A78 | none |

## Notes

Code edits relative to the original: lab absolute paths replaced by `FASTALAKE_DATA_DIR`; lab-tree roots located by a marker file replaced by `FASTALAKE_DATA_DIR`. No statistics changed. The scripts keep their original integrity gates (`verify()`, `gate()`, hash checks against `*_CODE.json` and `*_INPUTS.json` manifests of the frozen run). Those manifests contain absolute lab paths and are not in the repository; the hashes they held are listed in RECEIPTS.md. Because the code here carries the path edits described below, its own sha256 differs from the original (both are in RECEIPTS.md), so a re-run needs the manifests regenerated for this copy.The script also renders figure files (matplotlib, PyMuPDF) that are not used by the repository's figure units, which never redraw figures; only the `source_data/` tables are inputs here. See OPEN_ITEMS.md.
