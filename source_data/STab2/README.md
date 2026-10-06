# Source data: Supplementary Table 2

Inputs of Supplementary Table 2 (full-route shuffled-control acceptance and control
exposure) and of the Methods counts at v28d line 743. Written once by
`tables/STab2/prepare_source_data.py` from the original producer outputs; origins and
hashes are in `tables/STab2/PROVENANCE.md`.

Roster: 710 acquisitions (CAMPI 15, healthy stool 58, PREDICT 306, CapScan 331). This is
the roster of the full-route experiment as run on 2026-09-06; it does not include the
16 CapScan acquisitions of participant 15 that the v28f roster (347 CapScan, 726 FASTAs)
adds (see `tables/STab2/OPEN_ITEMS.md`).

Identifiers: `sample` is the raw-file stem of the acquisition, as in the public data
deposits and the Source Data package. No clinical covariate is in any file.

Comparison with the Source Data package: the published workbook
`Supplementary_Tables/Supplementary_Table_2.xlsx` (sha256 `4768d39f...`) is rebuilt from
these files by `tables/STab2/STab2.ipynb`; all 18 sheets agree cell by cell
(`tables/STab2/checks/table_compare.json`). `hstool_pooling_per_acquisition.tsv` is
byte-identical to package file `Figure3/per_acquisition.tsv` (map row RES-53).

Common columns of the full-route tables:

| column | unit | meaning |
|---|---|---|
| cohort | id | CAMPI_MOCK, HSTOOL (healthy stool), PREDICT or CAPSCAN |
| sample | id | acquisition (raw-file stem) |
| q_threshold | dimensionless | peptide q-value threshold (0.001, 0.005, 0.01, 0.02, 0.05) |
| min_residues | count of residues | lower length bound of the accepted peptides (5 or 9; upper bound 50) |

## per_acquisition.tsv
sha256: b473dab7ef4d992961ae9850f8c4f91668434249340984bbcf5c1b3ff7b5c9ad
rows: 7100 (one row per acquisition, q threshold and length floor; 710 x 5 x 2)
produced by: cohort_full_route_controls (collect_cohort.py, per cohort; concatenated)

| column | unit | meaning |
|---|---|---|
| n_canonical | count of peptides | accepted target-labelled canonical peptides (modification-stripped, I/L-collapsed) |
| n_target_only | count of peptides | accepted peptides assigned only to ordinary target proteins |
| n_shared | count of peptides | accepted peptides assigned to both target and control proteins |
| n_shuffle_only, n_plant_only, n_thermophilic_only | count of peptides | accepted peptides assigned only to shuffled, A. thaliana or P. furiosus control proteins |
| shuffle_only_fraction, plant_only_fraction, thermophilic_only_fraction | dimensionless | the count above divided by n_canonical; empty when n_canonical is 0 |
| shuffle_only_not_in_original_reservoir (and plant_, thermophilic_) | count of peptides | control-only peptides with no exact I/L-normalised substring match in the original stage-0 reservoir |
| ..._not_in_original_reservoir_fraction | dimensionless | that count divided by n_canonical |

## stage_control_counts.tsv
sha256: e2f7fe8bb4d32088c5d7680f04fa74d40ed315b5f286e084ba33c57c8265c65d
rows: 8520 (one row per acquisition, stage and control class; 710 x 4 x 3)
produced by: cohort_full_route_controls (collect_cohort.py)

| column | unit | meaning |
|---|---|---|
| stage | id | 0_control_bank, 1_cohort (cohort evidence selection), 3_acquisition (acquisition selection), 4_razor (razor FASTA) |
| control_class | id | shuffle, thermophilic or plant |
| n_proteins | count of sequences | control protein sequences present at that stage |

## baseline_comparison.tsv
sha256: f2f7298688a7a7057f3d49f187bce638c93ad1ed5c793d71c95b6ae452b604eb
rows: 710 (one row per acquisition)
produced by: cohort_full_route_controls (collect_cohort.py)

| column | unit | meaning |
|---|---|---|
| baseline_canonical | count of peptides | canonical target peptides at peptide q <= 0.01 in the ordinary search without controls |
| full_route_target_supported | count of peptides | canonical peptides at q <= 0.01 with a target assignment in the full-route search |
| intersection | count of peptides | peptides in both sets |
| lost | count of peptides | baseline peptides not in the full-route set |
| gained | count of peptides | full-route peptides not in the baseline set |
| count_ratio | dimensionless | full_route_target_supported / baseline_canonical |
| own_set_recovery | dimensionless | intersection / baseline_canonical |

## quantile_summary.tsv
sha256: 5be454201c9622a6b6fea724737d1523b94eae06f30270d527901b8736d22de2
rows: 240 (one row per cohort, q threshold, length floor, control class and metric)
produced by: cohort_full_route_controls (collect_cohort.py)

| column | unit | meaning |
|---|---|---|
| control_class | id | shuffle, thermophilic or plant |
| metric | id | name of the per-acquisition fraction summarised |
| n_attempted | count of acquisitions | acquisitions in the cohort |
| n_defined | count of acquisitions | acquisitions with a defined fraction (n_canonical > 0) |
| median, q1, q3, max | dimensionless | median, first and third quartile (Python `statistics.quantiles`, method exclusive) and maximum across acquisitions |

## control_discoveries_CAMPI_MOCK.tsv, _HSTOOL.tsv, _PREDICT.tsv, _CAPSCAN.tsv
sha256: 155c43af8f931cb0d15b98c9c385a4ab3cd9200e981eebb3739dacf0650cbe42 (CAMPI_MOCK); 62110ad8b152b282fdbb0cfca3a74cf0a470d6622f6db2d5860265fe360a9209 (HSTOOL); 7dafc42ccc5d25e6740c294da9fcce13b97207844ffd040dde8d40256f8b36dc (PREDICT); 8987ef72a5173d40abf3321e6bc30d0eceb1bb6602f93e83ba02231603138c21 (CAPSCAN)
rows: 175; 36742; 186062; 143766 (one row per accepted peptide with any control assignment, per acquisition and q threshold)
produced by: cohort_full_route_controls (collect_cohort.py). Kept per cohort because the concatenation exceeds 20 MB; the workbook sheet Control_discoveries is their concatenation in this order.

| column | unit | meaning |
|---|---|---|
| peptide | sequence | canonical peptide (modification-stripped, I/L-collapsed) |
| assignment | id | shuffle, plant, thermophilic (control only) or shared (also a target) |
| classes | id | semicolon-separated protein classes the peptide maps to |
| present_in_original_reservoir | 0/1 | 1 if the peptide occurs as an exact I/L-normalised substring in the original stage-0 reservoir |

## fusion_*.tsv and fusion_validation.json
sha256: fusion_cohort_q001 3c562fb726dd3a42c9455625279612181abfa68d2ccd0b6d798e492244d8f4d6; fusion_cohort_threshold_summary bb656c79a97c46fb88350c1a8f4fd8f0e745ff7a790ea7d32b50857c2709717f; fusion_matched_seed_summary bb656c79a97c46fb88350c1a8f4fd8f0e745ff7a790ea7d32b50857c2709717f (byte-identical to the previous file); fusion_short_5_to_8_summary 6711214fd4621e4c9a374a6d51d3bad6db9a884368df7e2480622d536d2239f2; fusion_acquisition_outliers_q001 439cb2763bd3f5d04c6f8343acfb955bdee2ebe0aecfc0990492fd27ce931e6d; fusion_completeness bfff2e36a17978ea94fd9634756f24a08f338a418aeb2eca1eb8dd01b82f5131; fusion_all_acquisitions c1fffd53310b355e6c69237d120eb7d76da62bc6dfb3278a963754193e718956; fusion_validation.json fe6c73b3a0aa177a28038066b6d4bb633f5aa355ac05cf640d2a1d6696fe6096
rows: 16; 80; 80; 40; 2840; 8; 14200 (fusion_all_acquisitions: one row per acquisition, seed, q threshold and length floor; 710 x 2 x 5 x 2)
produced by: cohort_fusion_controls (summaries written by audit_completed.py). Fusion-linked control experiment of 2026-09-08, seeds 20260907 and 20260908, 1,420 searches. These are the workbook sheets prefixed `fusion_`.

| column | unit | meaning |
|---|---|---|
| seed | id | fusion control permutation seed |
| n_canonical | count of peptides | accepted canonical peptides (N) |
| n_target_only, n_entrapment_only, n_shared | count of peptides | accepted peptides assigned to targets only, controls only, or both; they sum to n_canonical |
| n_unpaired_accepted | count of peptides | accepted peptides outside disjoint one-to-one target-control pairs (U) |
| n_paired_entrapment | count of peptides | accepted paired-control peptides (E) |
| paired_numerator | count of peptides | paired numerator without U (E + E-only + 2B in the producer's notation) |
| paired_upper_diagnostic_with_unpaired_worst_case | dimensionless | (paired_numerator + n_unpaired_accepted) / n_canonical |
| combined_upper_diagnostic_with_unpaired_worst_case | dimensionless | (2 n_paired_entrapment + n_unpaired_accepted) / n_canonical |
| paired_subset_fraction | dimensionless | as named by the producer |
| whole_procedure_fdr_control_established | boolean | always False |
| acquisitions | count of acquisitions | acquisitions pooled in a summary row |
| paired_pct, combined_pct, unknown_pct, entrapment_only_pct | % | pooled ratios of summed acquisition counts (100 (paired_numerator + U)/N, 100 (2E + U)/N, 100 U/N, 100 E-only/N) |
| paired_median_pct, paired_q25_pct, paired_q75_pct | % | median and 25th/75th percentiles (pandas `quantile`, linear) of the per-acquisition paired diagnostic |
| acquisitions_paired_above_1pct | count of acquisitions | acquisitions whose paired diagnostic exceeds 1% |
| completed, expected | count of searches | completed and expected searches per cohort and seed |

## indep_ordinary_counts.tsv, indep_fusion_all_thresholds.tsv, indep_validation.json
sha256: f0aa932eb0f90785f840eda9ddc2537d20602abf403bfd0d3cdb1c56b07adc44; 6c912a08e2ff8554d5ee807f3a9b46f033572902ee348b80c03dde41ab205a12; d7ed061abf46f4af38415ca0bb4911c202ecd8a8b1a2891b3347745ba4f1d9e2
rows: 18; 120
produced by: campi_independent_selection (collect_campi.py). CAMPI S01 and S02 against a 29,673-sequence reference; these are the workbook sheets prefixed `indep_`. Hashes of the two tables equal those recorded in indep_validation.json.

| column | unit | meaning |
|---|---|---|
| evaluated | id | evaluated CAMPI acquisition (S01 or S02) |
| arm | id | fixed (fixed database), self (selection from the evaluated acquisition), independent (selection only from the other acquisition) |
| max_residues | count of residues | upper length bound |
| canonical_target_peptides | count of peptides | accepted canonical target peptides at q_threshold |
| recovered_fixed | count of peptides | of those, also accepted in the fixed arm |
| lost_fixed | count of peptides | fixed-arm peptides not accepted in this arm |
| additional_vs_fixed | count of peptides | peptides accepted in this arm but not in the fixed arm |
| unit, interpretation | id | constant labels written by the producer |
| other columns | | as in the fusion tables |

## hstool_pooling_per_acquisition.tsv
sha256: 474de2e78bbfc37fcf132930b688069df06909699b6f7cb9a0c12bf0c42222f8
rows: 360 (one row per acquisition, arm, q threshold and length floor; 12 x 3 x 5 x 2)
produced by: healthy_stool_pooling_controls (reduce.py). Twelve healthy-stool technical-replicate acquisitions searched with full-route controls against their own database (`own`, reused from A91), the union of all twelve (`pool_all12`) and the union of the other eleven (`pool_other11`). Used for v28d line 743.

| column | unit | meaning |
|---|---|---|
| arm | id | own, pool_all12 or pool_other11 |
| other columns | | as in per_acquisition.tsv |
