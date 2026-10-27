# Source data: Supplementary Methods 8-11

Tables read by `text_statistics/supplementary_methods_8_11/supplementary_methods_8_11.ipynb`. They
were made from the original analysis outputs by
`text_statistics/supplementary_methods_8_11/prepare_source_data.py` (Slurm, 9 s, 1.4 GB peak);
origins and hashes are in that unit's `PROVENANCE.md`. The unit also reads tables owned by `Fig5`,
`Fig6`, `SFig10`, `SFig11`, `SFig14`, `STab3` and `STab5`; they are not duplicated here. Read floats
with `float_precision="round_trip"`. PREDICT patients appear only as Source Data pseudonyms
(`PREDICT_Pnnn`); no clinical covariate is in any table.

## own_stool_distance_matrices.tsv
sha256: 5c77ba103194697cd2c3db3bfcdf81a7f1eb058ea8826c45cc47a9517f8e9140
rows: 675 (one row per OG prevalence filter, capsule participant and stool participant; 3 x 15 x 15)
produced by: capscan_functional_profiles (distances recomputed with the producer's own code from its OG profile matrix)

| column | unit | meaning |
|---|---|---|
| prevalence | fraction of profiles | minimum fraction of the 75 participant-region profiles in which an OG is present (0, 0.1, 0.2) |
| capsule_participant | id | participant whose four capsule-region profiles are compared |
| stool_participant | id | participant whose stool profile is compared |
| mean_distance | Euclidean distance (Hellinger profiles) | mean over the four capsule regions of the distance to that stool profile |

## capscan_dna_runs.tsv
sha256: 31bcbe66f1c27a368e004a26207fd0e965c9220eb2f912a83476271ef2400eb5
rows: 302 (one row per CapScan metagenome run, public SRA accession)
produced by: capscan_metagenome_profiles; native_crossomics_profiles (native membership)

| column | unit | meaning |
|---|---|---|
| run_accession | id | SRA run accession |
| participant | id | participant pseudonym (1-15) |
| region | id | sampling group (duodenum, jejunum, ileum, colon = regions 1-4; stool) |
| functional_aggregate | boolean | run has positive signal in the functional DNA aggregate (OG, KO and pathway tables agree) |
| native_profile | boolean | run is one of the gene-sequence inputs of the native CapScan DNA profile (13 participants with all five groups) |

## predict_native_feature_summary.tsv
sha256: f2ff941609446bad7e5764b15c75ae064566b6302a594c84966d4496f62dafbc
rows: 4 (one row per feature filter and DNA or RNA layer)
produced by: native_crossomics_profiles

| column | unit | meaning |
|---|---|---|
| cohort, analysis, layer | id | PREDICT; sensitivity_all_features (all observed) or primary_prevalence10 (≥10% of profiles); DNA or RNA |
| input_features | count of sequences | encoded-sequence features before filtering |
| retained_features | count of sequences | features retained (all observed: positive in at least one of 77 profiles) |
| median_signal_fraction, minimum_signal_fraction, maximum_signal_fraction | fraction of TPM | retained share of each profile's sequence-resolved abundance |

## eds_prediction_model_metrics.tsv
sha256: df5b688227c6083df06dcf98c5f394a0413b8399cf4db552961333380a53ce7f
rows: 24 (one row per grouping rule, split seed and model)
produced by: predict_eds_heldout_prediction

| column | unit | meaning |
|---|---|---|
| grouping, seed, model | id | protein-grouping rule (original, simplified), outer split seed, model (protein, baseline, combined, distance) |
| n | count of patients | 85 |
| rho, MAE, r2 | dimensionless; EDS units; dimensionless | held-out Spearman rho, mean absolute error, R² (MAE and R² empty for the distance score) |
| minimum, maximum, outside_unit_interval | EDS units; count | range of held-out predictions; predictions outside 0-1 |
| permutation_p, permutation_q_8 | dimensionless | label-permutation P (first seed only) and BH q over eight tests |
| rho_ci_low, rho_ci_high | dimensionless | 95% bootstrap interval of rho (first seed only) |

## eds_prediction_permutations.tsv
sha256: 16e32d66a9d4fd63fb57a9ebf8893f5bf4fc978440e1b1cce2482e125c2c11a0
rows: 1592 (one row per grouping, model and label permutation; 8 x 199)
produced by: predict_eds_heldout_prediction

| column | unit | meaning |
|---|---|---|
| grouping, model | id | as above |
| permutation | id | permutation index (1-199) |
| rho | dimensionless | held-out Spearman rho with permuted labels |

## eds_prediction_heldout.tsv
sha256: 3a8082a5a5e080e85f404a9b71c14d80dd2181cc4222573148de040314626795
rows: 2040 (one row per patient, grouping, seed and model; 85 x 24)
produced by: predict_eds_heldout_prediction

| column | unit | meaning |
|---|---|---|
| row | id | patient position (0-84) in the analysis order; the bootstrap indexes this order |
| patient_pseudonym | id | Source Data pseudonym |
| dysbiosis_score | dimensionless | Enterotype Dysbiosis Score (target) |
| grouping, seed, model | id | as above |
| score | EDS units (distance model: distance) | held-out prediction |

## eds_prediction_increment.tsv
sha256: c5ba12b277ccbbac6e4b1a9b7aa347866fff74110915ba28dbab6871f2932d31
rows: 4 (one row per grouping and model, protein or combined)
produced by: predict_eds_heldout_prediction

| column | unit | meaning |
|---|---|---|
| grouping, model | id | as above |
| MAE_baseline_minus_model | EDS units | mean over patients of baseline absolute error minus model absolute error |
| ci_low, ci_high | EDS units | 95% bootstrap interval |
| bootstrap | text | resampling description |

## eds_prediction_selected_alphas.tsv
sha256: 127cbda331150bd5986beee763fb6049b4b237a6923eac5e9b65d5a46d8e19bc
rows: 90 (one row per grouping, seed, ridge model and outer fold)
produced by: predict_eds_heldout_prediction

| column | unit | meaning |
|---|---|---|
| grouping, seed, model, fold | id | as above; outer fold 0-4 |
| alpha | dimensionless | ridge penalty selected by inner cross-validation |

## eds_prediction_outer_folds.tsv
sha256: 9ea8807e33d63cff8e15092f4ca7521377a6a2c03a72be0acce4c2cc1b21e5b8
rows: 255 (one row per seed and patient; 3 x 85)
produced by: predict_eds_heldout_prediction

| column | unit | meaning |
|---|---|---|
| seed, fold | id | outer split seed and the outer fold in which the patient was held out |
| row | id | patient position (as in eds_prediction_heldout.tsv) |

## eds_prediction_fit_audit.tsv
sha256: b3cb9564e65bb4ec21422a81d43564b13e60b1da27a23f08246a0a6d97d63e0c
rows: 360 (one row per fit; both grouping rules are logged without a grouping column)
produced by: predict_eds_heldout_prediction

| column | unit | meaning |
|---|---|---|
| seed, model, outer_fold, stage | id | split seed, ridge model, outer fold, inner or outer fit |
| n_train, n_test | count of patients | training and evaluation patients of the fit |
| features | count of OGs | OG features retained in the training data |

## eds_prediction_distance_reference.tsv
sha256: 3693ae167a1ec1efa7ebb24adec98921de74ca489a015b481c965273497895cd
rows: 30 (one row per grouping, seed and outer fold)
produced by: predict_eds_heldout_prediction

| column | unit | meaning |
|---|---|---|
| grouping, seed, fold | id | as above |
| training_EDS_cutoff | dimensionless | lower-quartile EDS of the training patients |
| n_reference | count of patients | training patients at or below the cutoff (reference set of the distance score) |

## nonspecific_search_summary.tsv
sha256: 96f06a076849ce5d0b54cb68bd5bf5cdb3a30a1ad09bb230b2244cf97180b0bf
rows: 2 (one row per cohort)
produced by: capscan_nonspecific_search (CapScan, including participant 15), predict_nonspecific_search (PREDICT)

| column | unit | meaning |
|---|---|---|
| cohort | id | CAPSCAN or PREDICT |
| configured_searches, completed_searches | count of acquisitions | search configurations and searches with a passing completion receipt |
| cleave_at, semi_enzymatic | text; boolean | digestion settings (empty cleave_at = non-specific) |
| min_len, max_len | residues | peptide length range |
| peptide_min_mass, peptide_max_mass | Da | peptide mass range |
| precursor_tol_ppm, fragment_tol_ppm | ppm | symmetric tolerances |
| static_mod_C, variable_mod_M | Da | carbamidomethylation, methionine oxidation |
| max_variable_mods | count | variable modifications per peptide |
| lfq | boolean | label-free quantification |

Every value was checked to be identical in all acquisitions' configurations.

## hstool_grouping_receipt.tsv
sha256: d088b6f1401a962cc16fb5aaa57606b7d3183aa27f8eda5713f6f029f515cfe7
rows: 6 (one row per recorded count)
produced by: healthy_stool_study_grouping

| column | unit | meaning |
|---|---|---|
| key | id | acquisitions; quantified protein groups under direct peptide assignment and under original grouping; eligible LFQ features; the KO refit's quantified features and positive acquisitions |
| value | count | value recorded in the receipt |

## Corrected own-stool distances

`own_stool_distance_matrices_checked.tsv` is rebuilt from the corrected OG
profiles used by Figure 5. The original matrix is retained separately.
`text_statistics/supplementary_methods_8_11/prepare_checked_distances.py`
checks the input hashes, reconstructs all 675 distances, and verifies the own
and other-stool means against `source_data/Fig5/own_stool_participants.tsv`.
The adjacent JSON records the input and output hashes. The current statistics
use the checked matrix.
