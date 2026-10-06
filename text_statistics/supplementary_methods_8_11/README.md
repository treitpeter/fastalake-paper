# Supplementary Methods 8-11

Supplementary Methods 8-11 (v28d lines 779-911) describe the clinical, regional and cross-omics
analyses, the non-specific searches and histidine-boundary analyses, the statistical framework with
its testing-family map, and computational provenance. This unit recomputes the statistics and counts
these sections print, and checks every Benjamini-Hochberg family of the testing-family map (v28d line
893) by recomputing each adjusted q value from the raw P values within exactly the stated family. It
has no figure panels. Run with `python paper.py build supplementary_methods_8_11`.

Files: `supplementary_methods_8_11.ipynb` (the notebook), `stats.py` (all computations, constants and
the NUMBERS/STATS row builders), `NUMBERS.tsv`, `STATS.tsv`, `write_tables.py` (writes the two tables
once; the notebook asserts that it reproduces them), `prepare_source_data.py` (how this unit's
`source_data/` files were made from the original outputs; not on the reviewer path).

## Sections

| section | what is checked | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| CapScan regional human share | sign-flip P, BH q, bootstrap intervals recomputed; equal to producer to 0 | four regional contrasts, 32,768 patterns, 2,000 resamples (843) | 15 participants, paired | exact sign-flip, BH over 4 | CapScan regional human-assigned share | `source_data/Fig5/origin_fraction_cells.tsv`, `four_region_tests.tsv` | capscan_functional_profiles |
| CapScan own-stool test | permutation P, own-closer count, mean advantage, interval recomputed from distance matrices | "every participant closer ... P≤1 × 10⁻⁵" (847), 99,999 permutations (849) | 15 participants, 75 profiles | stool-label permutation, plus-one | CapScan own-stool permutation test | `own_stool_distance_matrices.tsv` (this unit), `source_data/Fig5/similarity_summary.tsv` | capscan_functional_profiles |
| Cross-omics families | family sizes 28/16 and 10/6, all q recomputed, populations, permutation count, classification stability | 809-839, 893 | PREDICT 77 patients; CapScan 13 participants (65 profiles), sensitivity 15 (75) | Spearman of distance matrices, permutation, BH | Cross-omics profile concordance: families | `source_data/SFig11/*_PREDICT.tsv`, `source_data/Fig5/{functional,native}_correlations.tsv` | predict_crossomics_refit, capscan_crossomics |
| Controlled distance sensitivity | 176 estimates, 24 families of 10 or 6, q recomputed | 819-821, 893 | 77 patients; 13 participants | permutation, BH per variant | Controlled distance sensitivity | `source_data/SFig11/controlled_176_comparisons.tsv` | crossomics_distance_sensitivity |
| Individual-function correlations | nine families, q recomputed, n ≥ 35 | 829, 893 | 35-77 patients per test | Spearman (t approximation), BH | Individual-function correlations | `source_data/STab5/crossomics_same_function_correlations.tsv` | predict_function_correlations |
| EDS screen | ten families per layer/normalization/model, q recomputed, n ≥ 43 | 799 | 43-85 patients per feature | Spearman or partial rank correlation, BH | EDS screen | `source_data/SFig10/all_EDS_associations.tsv` | predict_eds_screen |
| Held-out EDS prediction | held-out rho, permutation P, BH over 8, rho and MAE intervals, outer folds regenerated | 803-805 | 85 patients | label permutation (199), bootstrap (2,000) | Held-out EDS prediction | `eds_prediction_*.tsv` (this unit) | predict_eds_heldout_prediction |
| PREDICT antibiotic endpoints | two-test families, q recomputed, Monte Carlo P grid | 793-795, 893 | 42 patients (primary) | paired sign-flip, BH over 2 | PREDICT antibiotic endpoints | `source_data/STab3/endpoints_by_definition.tsv` | predict_rifaximin_sensitivity |
| Leading species and persistence | eligibility minimum (12 of 15), 24 functions, 111 of 360 profiles | 785-787 | 15 participants | none (descriptive) | CapScan leading species and Fig. 5e persistence | `source_data/Fig5/contributor_cells.tsv`, `mosaic_donor_regions.tsv` | peptide_centric_downstream |
| DNA runs and native features | 301 functional runs, 245 native runs in 13 complete participants; 14,822,078 / 5,391,526 features | 809, 827, 833 | runs; sequences | none | DNA runs and native sequence features | `capscan_dna_runs.tsv`, `predict_native_feature_summary.tsv` (this unit) | capscan_metagenome_profiles, native_crossomics_profiles |
| Non-specific searches, histidine families | completed searches and settings; 13 family sizes, q recomputed, exact and Monte Carlo P grids | 855, 869, 873, 893 | 347 + 306 acquisitions; 15 participants; 93 patients | sign-flip, BH | Supplementary Methods 9 | `nonspecific_search_summary.tsv` (this unit), `source_data/Fig6/*_statistical_families.tsv` | capscan_nonspecific_search, predict_nonspecific_search, terminal_participant_inference, terminal_histidine_inference |
| Healthy-stool grouping counts | 58 acquisitions, 59,428 groups, 1,909,973 features | 909 | acquisitions | none | Supplementary Methods 11 | `hstool_grouping_receipt.tsv` (this unit) | healthy_stool_study_grouping |

Recomputation depth. Regional tests, own-stool test, held-out EDS statistics and every BH q value are
recomputed from unit-level values. Family sizes, populations and counts are counted from row-level
tables. The 59,428, 1,909,973, 14,822,078 and 5,391,526 counts are read from producer summary tables
(the gene-level and peptide-level data are not in the repository). Correlations and permutation P values
of the cross-omics, controlled-distance and histidine tests are recomputed by the figure units named in
`STATS.tsv`, not here.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `BOOT_SEED_REGIONAL` | 20261003 | producer seed of regional and own-stool bootstrap (`capscan_fig5bcd.py` `boot()`); not stated in the text | none |
| `BOOT_DRAWS_REGIONAL` | 2000 | v28d line 843 | none |
| `OWN_STOOL_SEED` | 20260910 | v28d line 849 | none |
| `OWN_STOOL_PERMUTATIONS` | 99999 | v28d line 849 | none |
| `OG_PREVALENCES` | 0, 0.1, 0.2 | v28d line 847 (all observed is primary) | 0.1 and 0.2 are sensitivities, same table |
| `Q_THRESHOLD` | 0.05 | v28d line 889 | none |
| `CROSSOMICS_SEED` | 20260907 | v28d line 815 (native); same seed in the producer for functional profiles, not stated at line 827 | none |
| `CROSSOMICS_PERMUTATIONS` | 1999 | v28d lines 815, 827, 835 | none |
| `EDS_SPLIT_SEEDS` | 20260915, 20260916, 20260917 | v28d line 803 | seeds 2 and 3 are sensitivities |
| `EDS_LABEL_PERMUTATION_SEED` | 20260915 | producer (`analyze_score.py` `main()`); not stated | none |
| `EDS_BOOT_SEED` | 20260918 | producer (`analyze_score.py` `main()`); not stated | none |
| `EDS_BOOT_DRAWS` | 2000 | v28d line 805 | none |
| `EDS_OUTER_FOLDS`, `EDS_INNER_FOLDS` | 5, 3 | v28d line 803 | none |
| `EDS_PATIENTS` | 85 | v28d line 803 | none |
| `EDS_ALPHAS` | 0.01, 0.1, 1, 10, 100 | v28d line 803 | none |
| `ANTIBIOTIC_SIGNFLIP_DRAWS` | 200000 | v28d line 795 | exact enumeration (Codex reimplementation, see PROVENANCE.md) |
| `HISTIDINE_MC_DRAWS` | 200000 | v28d line 869 (n > 16) | exact enumeration for n ≤ 16 |
| `REG` | duodenum, jejunum, ileum, colon, stool | producer region ids of regions 1-4 and stool | none |

## Numbers on these lines recomputed by other units

These printed numbers are recomputed in the unit named (snapshot of their `NUMBERS.tsv`, 6 October
2026) and are not repeated here.

| v28d line | printed | unit | status there |
|---|---|---|---|
| 785 | 1,467; 21,168 | SFig12 | MATCH |
| 793 | 12, 12, 35, 38, 40, 42, seven | STab3 | MATCH (controlled inputs) |
| 799 | 1,000; 16; 192; 20; 21; 43; 84; 85; <0.05; five; three | SFig10, STab4 | MATCH |
| 799 | 50; six; 95% | SFig10, STab4 | OPEN |
| 801 | 2; 20; 259; 85; Nineteen; six; three | SFig10 | MATCH |
| 813 | 26.5%; 21.0% | SFig11 | MATCH |
| 823 | 11 of 24; three; four; two; two; 24 | SFig11 | v28f_change: 9 of 24, one; others MATCH |
| 829 | 0.10-0.12; 30,581; 35; nine | STab5 | MATCH |
| 835 | 1,999; 16; thirteen | SFig11 | MATCH |
| 847 | 73,316; 53,023; 44,330 | Fig5 | MATCH |
| 849 | 10.6784%; 4.1584%; 5.2590%; 14.8369%; 60 | Fig5 | MATCH |
| 863 | fourteen | Fig6 | v28f_change: fifteen |
| 869 | 2,000; 20,000; 200,000; 50,000; fifteen | Fig6 | MATCH |
| 869 | fourteen | Fig6 | v28f_change: fifteen |
| 869 | 16 | Fig6 | OPEN |
| 903 | 15; 58; 306; 368; 97 | SFig15 | MATCH |
| 903 | 331 | SFig15 | v28f_change: 347 |

Software versions and reference numbers on line 907 are not statistics and have no rows.

## Known limitations

- 25 printed numbers are OPEN (`OPEN_ITEMS.md`): mostly analysis parameters (thresholds, seeds, the
  ridge grid) that the producers use but that no repository code exercises, plus three counts whose
  source was not located.
- The permutation seed 20260910 cannot be confirmed from the output: no permuted value was as extreme
  as observed, so P equals 1/(99,999 + 1) for any seed.
- The antibiotic sign-flip tests and bootstrap intervals need controlled clinical inputs; they are
  recomputed in `tables/STab3`, and only the family structure is checked here.
- Tables owned by other units are read in place; if those units change a table, this notebook's
  assertions will fail rather than silently use different data.
