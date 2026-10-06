# Supplementary Table 5

Individual-function correlations between DNA, RNA and protein layers in the PREDICT
cohort (77 patients). For each functional representation (OG, KO, pathway) and each pair
of microbial layers (metagenome DNA, metatranscriptome RNA, non-host protein), one
Spearman correlation across patients per function: 30,581 tests in nine
Benjamini-Hochberg families. The table also lists the eligibility of all 184,386
candidate function-by-layer-pair rows and compares the revised functional allocation with
the previous one for functions eligible under both. The notebook rebuilds the workbook,
compares it cell by cell with the published table, and recomputes every test from the
patient-level profiles.

## Sheets

| sheet | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| Read_me | sheet descriptions | 483 | NA | NA | Table | NA (text in panels.PACKAGE_README) | STab5_workbook_build |
| All_30581_tests | per function and layer pair: patients, rho, P, q | 141, 829 | 77 patients; 35 to 77 per test | Spearman, two-sided; BH within family | Table; Recomputation of the 30,581 tests | crossomics_same_function_correlations.tsv; profiles_*.tsv.gz | predict_function_correlations |
| Nine_family_summary | tests, median rho, tests with q < 0.05 per family | 141 | 9 families | summary of the above | Table; Recomputation of the 30,581 tests | same_function_summary.tsv; profiles_*.tsv.gz | predict_function_correlations |
| All_eligibility | every candidate function and layer pair, patients positive in both layers | 829 | 184,386 rows | eligibility rule (at least 35 patients) | Table; Recomputation of the 30,581 tests | eligibility.tsv; profiles_*.tsv.gz | predict_function_correlations |
| Eligibility_common | common, newly eligible and lost functions per family; median rho revised and previous | 141, 829 | 9 families | none (descriptive) | Table | crossomics_same_function_correlations.tsv; previous_allocation_same_function_correlations.tsv | STab5_workbook_build |

The two Source Data support files (`Supplementary_Tables/Table5_support/`
`eligibility_matched_correlations.tsv`, `common_eligibility_summary.tsv`) are written by
the same notebook cell and checked byte for byte.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| N_PATIENTS | 77 | patients with complete DNA, RNA, protein and metabolite data, one earliest complete specimen each (Supplementary Methods 8, v28d lines 809, 829) | none |
| MIN_PATIENTS | 35 | minimum patients with positive values in both layers (Supplementary Methods 8, v28d line 829) | none reported |
| LEVELS | OG, KO, Pathway | the three functional representations (Supplementary Methods 5 and 8, v28d lines 657, 829) | none |
| LAYERS | metaG, metaT, microbial_metaP | all three pairs tested (Supplementary Methods 8, v28d line 829) | none |
| CORRELATION | spearman | rank correlation across patients, relative signal over the layer pair's shared function vocabulary (Supplementary Methods 8, v28d line 829) | none |
| CORRECTION | fdr_bh, per level x layer pair | nine families (Supplementary Methods 8, v28d line 829) | none |
| Q_SIGNIFICANT | 0.05 | threshold of the `significant` column (producer) | none |
| MATCH_KEYS | level, layer_a, layer_b, feature | matching of revised and previous allocation (Supplementary Methods 8, v28d line 829) | none |
| FLOAT_RTOL | 1e-12 | tolerance of the producer's own workbook comparison | none; all 1,351,138 cells are exactly equal |

## Known limitations

- These tests were not repeated under the allocation-sensitivity rules of Supplementary
  Methods 5; the manuscript states this (v28d line 657).
- Profile normalization, eligibility and families follow the original code exactly; the
  P value is the t approximation of `scipy.stats.spearmanr` (n - 2 degrees of freedom),
  not an exact or permutation P.
- The previous-allocation test table (17,914 tests) is shipped as produced and is not
  recomputed here; its protein profiles are not in this repository.
- `prepare_source_data.py` documents how `source_data/STab5/` was made from the lab files
  and two controlled inputs; reviewers do not need to run it.
