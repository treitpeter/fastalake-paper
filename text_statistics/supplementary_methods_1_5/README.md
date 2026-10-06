# Supplementary Methods 1-5

Supplementary Methods 1-5 (v28d lines 493-659) describe the datasets, sequence selection, reference comparisons, database searching, protein inference and annotation. They print counts, percentages and the results of the functional-allocation sensitivity analysis. This unit recomputes those numbers and tests from `source_data/`. The text has no figure panels.

Scope: the statistics listed in the table below. Configuration constants (for example tolerances and length windows), roster sizes and reservoir sizes printed in the same lines are covered by the units that own them (`methods`, `Fig1`, `Fig2`, `SFig01`, `STab6`); see `OPEN_ITEMS.md`.

Run from the repository root: `python paper.py build supplementary_methods_1_5`. To rebuild the tables in `source_data/supplementary_methods_1_5/` from the original producer outputs (lab access needed): `FASTALAKE_LAB_ROOT=<lab project directory> python text_statistics/supplementary_methods_1_5/prepare_source_data.py`. To regenerate `NUMBERS.tsv` and `STATS.tsv`: `python text_statistics/supplementary_methods_1_5/stats.py --write`.

## Statistics

| statistic | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| AlphaNovo permuted control | exact-match rate after permuting predictions across spectra | "reduced the match rate to 0.006%" (553) | 801,613 reference identifications (peptide) in 15 acquisitions | none (negative control) | Statistic: AlphaNovo permuted-prediction control | source_data/Fig2/alphanovo_accuracy_per_acquisition.tsv | campi_alphanovo_accuracy |
| bag-aware lookup | gains and losses against exact lookup | "added 6,961 ... removed 5,941 (net +1,020, positive in 38 of 48)" and following (561, 579) | 48 acquisitions; 336 arm-acquisition rows | none (counts) | Statistic: bag-aware and mass-equivalent lookup | source_data/SFig02/bags_bounded48_per_acquisition.tsv, bags_bounded48_orderings.tsv | bag_mass_lookup_search |
| study-supported space | proteins retrieved by pooled study queries | "18,007,876 proteins" ... "26,881,038" (567) | reservoir proteins; 58 and 97 pooled acquisitions | none | Statistic: healthy-stool study-supported protein space | source_data/Fig1/stage1_stats.json, hstool97_stage1_stats.tsv, hstool97_pooled_prediction_files.tsv | healthy_stool_stage1_lookup, healthy_stool_97_study_lookup |
| parser audit | archived and current parser agree | "all 710 retained prediction files (112,835,982 rows)" (573) | 710 (v28f 726) prediction files | none | Statistic: parser audit | parser_audit_acquisitions_726.tsv | prediction_parser_audit |
| CAMPI target databases | records and contaminants | "29,673 reference and 19,435 MG records, each including 116 contaminants" (609) | FASTA records | none | Statistic: CAMPI target databases and member presence | campi_target_fasta_records.tsv | none (supplied files) |
| member presence | exact presence in the reservoir | "eight ... (29,557 member proteins)" (611) | 29,557 proteins in 8 proteomes | none | Statistic: CAMPI target databases and member presence | source_data/Fig2/reservoir_protein_presence_by_proteome.tsv | campi_reservoir_containment |
| input completeness | closing tag and scan coverage of every search input | "All passed (lowest coverage 99.6%)" (629) | 637 (v28f 653) + 58 searches (acquisition) | none (check) | Statistic: search input completeness | search_input_completeness.tsv | sage_input_completeness |
| partial-origin human share | human-assigned intensity with partly unannotated proteins | "0.12% ... 0.60% ... 1.96%" (651) | LFQ intensity pooled over 58, 331 (v28f 347), 306 acquisitions | none | Statistic: human-assigned share with partly unannotated compatible proteins | origin_rows_*.tsv | origin_annotation_share |
| annotation checks | seed-ortholog check, MGnify coverage, CAMPI supplement | "Of 9,368,708 MGnify annotations, 372,277 (4.0%) failed ..." (657) | gene annotations (protein) | none | Statistic: annotation sources and the seed-ortholog check | mgnify_seed_ortholog_check_by_genome.tsv, annotation_mgnify_coverage_by_arm.tsv, campi_supplemental_annotation_queries.tsv | mgnify_seed_ortholog_check, annotation_repair_coverage, campi_supplemental_eggnog |
| allocated share | median allocated intensity per rule | "7-61% instead of 62-81%" (657) | acquisitions (58, 306, 346) | none | Statistic: functional-allocation sensitivity, allocated share | allocation_coverage_*.tsv | functional_allocation_sensitivity |
| own stool | CapScan intestinal samples closer to own stool | "All fourteen CapScan participants ... (P≤1 × 10⁻⁵)" (657) | 14 (v28f 15) participants | permutation test, two-sided | Statistic: functional-allocation sensitivity, CapScan own stool | allocation_capscan_own_stool_distances_347.tsv | functional_allocation_sensitivity |
| PREDICT cross-omics | protein-profile distance correlations | "ρ=0.16-0.64 against 0.16-0.51 ... all P≤0.026"; "0.43-0.51 to 0.28-0.35 (P≤0.002)" (657) | 77 patients | Mantel test (Spearman), two-sided | Statistic: functional-allocation sensitivity, PREDICT cross-omics | allocation_predict_distances.tsv | functional_allocation_sensitivity, predict_crossomics_refit |
| healthy-stool KO r | median pairwise correlation by design | "(r=0.89-0.91)" ... "0.81 to 0.69-0.73" (657) | 66 and 45 acquisition pairs | none (descriptive) | Statistic: functional-allocation sensitivity, healthy-stool KO correlation | allocation_hstool_ko_pairs.tsv | functional_allocation_sensitivity, healthy_stool_ko_precision_refit |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| ACCURACY_WINDOW | 9-50 | length window of the accuracy and control statement (line 553) | 5-50 in source_data/Fig2/alphanovo_accuracy_per_acquisition.tsv |
| BAG_BOUNDED_ARM, ORDERING_BOUND | bags_le8, 8 | the setting named in line 561 | bags_le2, bags_le4, bags_all, mass_all, comb_all in the same table |
| PARSER_ORIGINAL_FILES | 710 | the v28d audit set; indices 710-725 are the 16 participant-15 files added for v28f | both sets reported in NUMBERS.tsv |
| COMPLETENESS_THRESHOLD_TEXT | 0.95 | criterion stated in line 629 | 0.90 used by the original check (COMPLETENESS_THRESHOLD_CODE); every search passes both |
| SEED_ORTHOLOG_MIN_RATIO | 0.5 | "at least half of the reported bitscore" (line 657) | 0.2 and 0.8 in the producer summary (not shipped) |
| RESTRICTED_RULES | unique, complete, unique_complete | the three rules of line 657 | equal allocation (production rule) |
| OWN_STOOL_SEED, OWN_STOOL_DRAWS | 20260910, 99,999 | producer settings (alloc_capscan347.py) | none |
| PREVALENCE_FILTERS | 0, 0.1, 0.2 | producer settings; text pools all three | per-filter rows in STATS.tsv |
| MANTEL_SEED, MANTEL_PERMUTATIONS | 20260907, 1,999 | producer settings (check_prevalence.py) | none |

## Known limitations

- The allocation-sensitivity tests are re-run here from intermediate distance matrices and per-pair correlations that `prepare_source_data.py` derived with the producers' own functions; the profile matrices themselves (up to 42 MB each) are not in the repository.
- The P value of the own-stool test is at the smallest value the permutation scheme can return (1 / 100,000).
- The healthy-stool KO correlation ranges pool the profile with and without human-assigned signal; the text does not say so.
- The 26,881,038-protein space is the output of the original production run, read from its statistics file; it was not recomputed.
