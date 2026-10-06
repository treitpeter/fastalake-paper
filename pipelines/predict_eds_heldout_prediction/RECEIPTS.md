# RECEIPTS

## Original Slurm jobs

| role | job |
|---|---|
| baseline replay | 6044249 |
| refit | 6048399 |

## Code

| original (LAB_ROOT-relative) | file here | sha256 original | sha256 here |
|---|---|---|---|
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_eds.py` | `nested_eds.py` | `1551285b6b1164c179682e344edae519c87f36c97e8857354bddecc4fb0d8bf5` | `1551285b6b1164c179682e344edae519c87f36c97e8857354bddecc4fb0d8bf5` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/code/analyze_score.py` | `analyze_score.py` | `bc24a718c5e6913e5b3fbc6105c328faf36261e0d4c085b52cab27752f48f0bb` | `6c8b641acb7e48077e82c0d9a1e6c6ed4c6832d8e0bef189bd1c013eae4310fc` |

## Inputs (sha256 recorded by the original run)

| input (LAB_ROOT-relative) | sha256 |
|---|---|
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/code/analyze_score.py` | `bc24a718c5e6913e5b3fbc6105c328faf36261e0d4c085b52cab27752f48f0bb` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/A19_six_models_reproduced.tsv` | `2587f378c1019b7dd5c7dd11c267a89d98dfd6fadb92792384e65f220eef3c7d` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/OG_reconstruction.tsv` | `dc0cd475d1d9caba05dedad4282d101d3612e9b49afc9be08a7dcf0b7637a5c6` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/distance_reference_folds.tsv` | `f3127cfeaac2575b63149fa4520aafea91ce62d8e240b4a171ac6a6d938a36eb` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/fit_boundary_audit.tsv` | `4c9d5f9f27cdb719e6370b0fa87a34c42dfd3e4dddea4fa96a37c24f2fe628b3` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/fixed_heatmap_review.tsv` | `f2d23b99334258bb84b8fdbafb6e5d4eec17d439678c243c230f2c4f4809a9c0` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/grouping_score_agreement.tsv` | `d9e345101d301b6bf7a29adcfeba4a88dbc9909c64940366b7592dc7279bfaeb` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/heldout_scores.tsv` | `56a0e6fde92da121995eeb38acf0823c574c2a1ac614690e4c45203aec5fb8a6` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/historical_model_reproduction.tsv` | `997f7449ab9f02a83d1aec7ddd443860b4872e3a086747200963d17d0d67ae9f` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/increment_over_baseline.tsv` | `492812b7a8978ab967aed78f0f465fda751c3defc3098d164518a9f4939224cc` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/model_metrics.tsv` | `03f26d11826110e2c6ad87b90c4a72254a47c142e8c54daa497a9b1a534e76c8` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/outer_folds.tsv` | `731ae0bff69de204498168c4cad4e6406fe231840393af178a44774188cbdfe7` |
| controlled input `patient_roster.tsv` (not in repository) | `aa4f4ba19c5f39a9715b7c8d29b808640c1d9b790dc5ec3821630fa6f1bb187c` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/permutation_results.tsv` | `c7db39ace48c9b02aed476a0fbf66323295338d192b87b123fd6ddc030e1c643` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/phenotype_associations.tsv` | `e9f6bba70bd29a985c201adbe3886c6307d3158f628676bc9e53e32c29c94e94` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/score_context.tsv` | `4fc149dda74d0510e3656d242f30375b2abf808d7269a4fa5088b9564609c942` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/selected_alphas.tsv` | `0968e48db95cd10ac5355a341b80c15343716e63702854f355f6c45e428e7adc` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/full/source_data/training_mean_baseline.tsv` | `adaae6059c15c01b16424bd655c00e6ea583ac6e086be23249c107e2530271e3` |
| controlled input `EDS_context_patient_roster.tsv` (not in repository) | `d91cd9be292091f380e7bdeb7285adaf9ef64fc40263ba89c6c4e5422d47c8a1` |
| controlled input `acquisition_accounting.tsv` (not in repository) | `5f7dac66d38cca1fc72d9899771b3bcb65fab9ca02f5328482ee856d142c55b6` |
| controlled input `heatmap_85_patient_order.tsv` (not in repository) | `4ee42bf87fdc5b73878ff70d20c2f0d995bacd8e35ea378baf5738a45256384d` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/inputs/PREDICT__original/source_data/PREDICT_noclust/microbial_OG.tsv` | `d468652ce315291e69a23f57027974a96f4c0d7a23a08ff253ebfc705e6d16c7` |
| controlled input `EDS_context_patient_roster.tsv` (not in repository) | `86460216d07ecc640c981900e17308d4d7f5a706868817a3fc1cb3237e025429` |
| controlled input `acquisition_accounting.tsv` (not in repository) | `1a2dac477320d9d7362b572db683369adab168879df0dab40ce901fe621a43a5` |
| controlled input `heatmap_85_patient_order.tsv` (not in repository) | `fb5abc826fba293df5628898c4cea9296645f6efbdac9edf7aa2763633f2446d` |
| `MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/inputs/PREDICT__simplified/source_data/PREDICT_noclust/microbial_OG.tsv` | `0b7250aa87327e182a7f0a64511569a2a7be1ab0f592b4384115ea998c66f30b` |
| controlled input `specimen_predictions_and_metadata.tsv` (not in repository) | `a8e22468e8b4501b3e2b93d9f8f522b3a8733caa6255e57ea4233e764610bfb6` |

## Outputs (sha256 recorded by the original run)

| output (LAB_ROOT-relative) | sha256 |
|---|---|
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/outer_folds.tsv` | `731ae0bff69de204498168c4cad4e6406fe231840393af178a44774188cbdfe7` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/grouping_score_agreement.tsv` | `e6b9f9eb3b32ca26a328decf4d9d6850d062bdc356dbd966efeb572a2555a5a4` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/phenotype_associations.tsv` | `e9f6bba70bd29a985c201adbe3886c6307d3158f628676bc9e53e32c29c94e94` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/patient_roster.tsv` | `aa4f4ba19c5f39a9715b7c8d29b808640c1d9b790dc5ec3821630fa6f1bb187c` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/fit_boundary_audit.tsv` | `4c9d5f9f27cdb719e6370b0fa87a34c42dfd3e4dddea4fa96a37c24f2fe628b3` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/INPUT_PROVENANCE.json` | `a54ea5fe46a66130e80718006f733f567d6cd9ac01dacfdd99e9a6b249ed7701` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/COMPUTATION_VALIDATION.json` | `dd3c623b5899d2dae577680752b37f0b7b448d5dbf247b42d4331b59000f10df` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/score_context.tsv` | `4fc149dda74d0510e3656d242f30375b2abf808d7269a4fa5088b9564609c942` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/selected_alphas.tsv` | `0968e48db95cd10ac5355a341b80c15343716e63702854f355f6c45e428e7adc` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/historical_model_reproduction.tsv` | `997f7449ab9f02a83d1aec7ddd443860b4872e3a086747200963d17d0d67ae9f` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/heldout_scores.tsv` | `e3c1a9b23d22cd6ff52e910d6f57a3e54932b920aabf52f983a027446157075c` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/model_metrics.tsv` | `9acf62a34aabc6b462fbea5cee218a1b21f3d89f4282da34ddd2b14dc02f472a` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/increment_over_baseline.tsv` | `2122ac3e4a0c6110f51cbc3b78af4969f53d8ca50ff95f5058c56281bd8a4151` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/distance_reference_folds.tsv` | `f3127cfeaac2575b63149fa4520aafea91ce62d8e240b4a171ac6a6d938a36eb` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_baseline/permutation_results.tsv` | `c7db39ace48c9b02aed476a0fbf66323295338d192b87b123fd6ddc030e1c643` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/outer_folds.tsv` | `731ae0bff69de204498168c4cad4e6406fe231840393af178a44774188cbdfe7` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/grouping_score_agreement.tsv` | `f0755e572fa1ff817b3bb785291e4ccbaa2e001244e42777fe0ef6f7cecc293c` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/phenotype_associations.tsv` | `15a8ba2125cfa91edbc070128e9dd93db6575956c5d53932ac1d0fa6e0f0b167` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/patient_roster.tsv` | `6664b7dbfdb84c84370ab3d09cde2e57f23b8df237b8302414f0b1eab6b161de` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/fit_boundary_audit.tsv` | `cf2115b8c1c8325f560a4364c9c40856a665be25d821d44b2898cf8268c7a752` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/INPUT_PROVENANCE.json` | `a54ea5fe46a66130e80718006f733f567d6cd9ac01dacfdd99e9a6b249ed7701` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/COMPUTATION_VALIDATION.json` | `49f0e10bce6e511638e5a0566f6ab9b4b01c5d0a8baf71cf590bc0ff65f1e245` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/score_context.tsv` | `4b38a4863fe860f9eb694d87a8902a20604b96d1eb181898b2a45e22fd9302eb` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/selected_alphas.tsv` | `127cbda331150bd5986beee763fb6049b4b237a6923eac5e9b65d5a46d8e19bc` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/historical_model_reproduction.tsv` | `997f7449ab9f02a83d1aec7ddd443860b4872e3a086747200963d17d0d67ae9f` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/heldout_scores.tsv` | `e4d6d384d74be1268606a19e6cdc052c05b2666177e0110b7d6d3e2121214483` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/method_comparison.tsv` | `c9fcc705106d092a4b48665633b55427d7492a5e340d2e6e9f9651ef7237c06e` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/model_metrics.tsv` | `df5b688227c6083df06dcf98c5f394a0413b8399cf4db552961333380a53ce7f` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/increment_over_baseline.tsv` | `c5ba12b277ccbbac6e4b1a9b7aa347866fff74110915ba28dbab6871f2932d31` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/distance_reference_folds.tsv` | `f3127cfeaac2575b63149fa4520aafea91ce62d8e240b4a171ac6a6d938a36eb` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/permutation_results.tsv` | `16e32d66a9d4fd63fb57a9ebf8893f5bf4fc978440e1b1cce2482e125c2c11a0` |

## Path aliases

Lab path segments are shortened in code and tables (`jess_integration_2026-09-29` to `integration_2026-09-29`, `v10_mm_philip_comment_audit_2026-09-22` to `v10_comment_audit_2026-09-22`, `coauthor_comment_check_v25_2026-10-05` to `comment_check_v25_2026-10-05`). `FASTALAKE_DATA_DIR` must hold the frozen inputs under the shortened names; the original names are the lab directory names.

## Four-test BH recompute (v30)

Slurm job 6074462 (p.hpcl92, 2 cores, 2G). Code: `eds_bh_direct4.py` = `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/grouping_recount/eds_bh.py` (sha256 `fef8ac935a68d923a725b9e2622302456eb54099f93df289cc542b6a2b2c171d`) with a docstring added.

| file | role | sha256 |
|---|---|---|
| LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/nested_refit/model_metrics.tsv | input | df5b688227c6083df06dcf98c5f394a0413b8399cf4db552961333380a53ce7f |
| LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/grouping_recount/eds_bh_direct4.tsv | output | 36263b74f2b515dc86ef02e75c61bd5534ef7f50851fb854ed735a1bb19d4d5f |
