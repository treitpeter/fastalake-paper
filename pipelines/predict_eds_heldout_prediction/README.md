# predict_eds_heldout_prediction

Internal ID: A108 (previous folder name A108_eds_heldout_prediction)

How well do OG function profiles and protein-group depth predict the dysbiosis score in 85 patients under nested cross-validation, against label permutations and a depth-only baseline? Used by the held-out prediction text of the supplementary methods (`eds_prediction_*.tsv`: model metrics, held-out scores, permutations, folds).

Status: code copied as it ran (path edits only), not re-run in this repository. See RECEIPTS.md for original job ids and hashes.

## Inputs

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

## Outputs

`nested_refit/{model_metrics,heldout_scores,permutation_results,selected_alphas,outer_folds,increment_over_baseline,phenotype_associations}.tsv`, `method_comparison.tsv` (hashes in RECEIPTS.md).

## Slurm resources and run time (original logs, `sacct`)

| step | job | cores | requested mem | peak RSS | wall |
|---|---|---|---|---|---|
| baseline replay | 6044249 | 2 | 8G | 2.4 MB | 16 s |
| refit (one of five steps) | 6048399 | 4 | 16G | 1.22 GB whole job | 4 min 9 s whole job |

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| outer / inner folds | 5 / 3 | protocol set before the run | none |
| split seeds | 3 | stability across splits | none |
| ridge alphas | 0.01, 0.1, 1, 10, 100 | grid fixed in advance | none |
| label permutations | 199 | protocol | none |
| fixed-prediction resamples | 2000 | protocol | none |
| patients | 85 | fixed EDS roster | none |
| covariates | log10 protein groups (original and simplified grouping) and human-assigned signal | both depth definitions kept | none |
| function profile | OG other-unresolved profile of the verified-orthologue rerun (A88) | replaces the earlier profile | baseline mode replays the earlier run to tolerance 1e-7 |

## Notes

Code edits relative to the original: lab absolute paths replaced by `FASTALAKE_DATA_DIR`; lab-tree roots located by a marker file replaced by `FASTALAKE_DATA_DIR`. No statistics changed. The scripts keep their original integrity gates (`verify()`, `gate()`, hash checks against `*_CODE.json` and `*_INPUTS.json` manifests of the frozen run). Those manifests contain absolute lab paths and are not in the repository; the hashes they held are listed in RECEIPTS.md. Because the code here carries the path edits described below, its own sha256 differs from the original (both are in RECEIPTS.md), so a re-run needs the manifests regenerated for this copy.`nested_eds.py` imports `analyze_score.py` from the original eds_phenotype folder; here a copy sits beside it, and the driver loads it from `N/'code'` as in the original (path `eds_phenotype/code/analyze_score.py` under the data directory). Uses `common.py` from A70.

## Four-test multiple-testing family (v30)

With the original protein grouping removed, the EDS held-out prediction family is the four
direct-assignment model tests instead of eight grouping-by-model tests.
`eds_bh_direct4.py` reads `model_metrics.tsv` of the nested refit, recomputes the stored
eight-test Benjamini–Hochberg q values as a control (maximum absolute difference 7e-17) and
computes q over the four `grouping == simplified` rows (seed 20260915; plus-one
permutation P, 199 permutations; unit: patient, n = 85).

| model | rho | permutation P | q (8 tests) | q (4 tests) |
|---|---|---|---|---|
| protein | 0.644 | 0.005 | 0.0067 | 0.0067 |
| baseline (depth + human share) | 0.351 | 0.050 | 0.050 | 0.050 |
| combined | 0.646 | 0.005 | 0.0067 | 0.0067 |
| distance | 0.367 | 0.005 | 0.0067 | 0.0067 |

No q value changes; only the stated family size changes. Launcher: `slurm/eds_bh_direct4.sh`
(`FASTALAKE_NESTED_REFIT_DIR`, `FASTALAKE_OUT_DIR`, `FASTALAKE_PYTHON`). Resources: 2 cores,
2 GB requested, seconds. Job id and hashes in RECEIPTS.md.
