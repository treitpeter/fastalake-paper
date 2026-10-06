# Supplementary Figure 11

Shows whether samples that differ in one molecular layer also differ in another, at the level of
KEGG pathway profiles, and how the native comparison depends on feature filtering and on the
distance and correlation method. For every pair of layers (DNA, RNA, protein excluding
human-assigned signal, human-assigned protein, metabolites), the Euclidean distances between all
pairs of PREDICT patients (a) or CapScan participant-by-group profiles (b) are correlated. Panels c
to f repeat the native comparison with a 10% prevalence filter and with other distance and
correlation methods. In v30 the native, eggNOG OG and KO matrices of the former panels a and b were
removed because they show the values of Fig. 4i and Fig. 5g; their tests are still recomputed here
(STATS.tsv rows labelled "Source Data") and kept in Source Data.

In this unit the column `v28d_line` of NUMBERS.tsv and STATS.tsv holds the paragraph index of the
v30 text (`visible_v30f.txt`, first column), not a v28d line. Every quote was checked against that
text with `tools/check_numbers.py --manuscript-text` (one paragraph per line).

## Panels

All six panels are drawn on one canvas by the v30 build; the notebook renders it once under
Panel a.

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | PREDICT KEGG pathway distance correlations; all observed features | pathway RNA-protein ρ=0.16, q=0.024 (v30 paragraph 74); legend (241) | 77 patients | Spearman correlation of distance-matrix upper triangles, 1,999 patient-label permutations, two-sided; BH within the 28 OG/KO/pathway comparisons | Panel a | `source_data/SFig11/functional_correlations_PREDICT.tsv`, `distances_PREDICT_functional.tsv` | predict_crossomics_refit (not in `pipelines/`) |
| b | CapScan KEGG pathway distance correlations, region-adjusted | legend (241); BH family of 16, pathways in Supplementary Fig. 11b (427) | 13 participants (65 participant-by-group profiles) | as a, 1,999 whole-participant block permutations; BH within 16 OG/KO/pathway comparisons | Panel b | `source_data/Fig5/functional_correlations.tsv`, `source_data/Fig5/functional_distances_long.tsv` | capscan_crossomics (not in `pipelines/`) |
| c | PREDICT native correlations, all observed vs ≥10% prevalence | legend (241); retained signal 26.5% and 21.0% (416) | 77 patients | as a; BH within each analysis (10) | Panel c | `native_correlations_PREDICT.tsv`, `distances_PREDICT_native.tsv`, `predict_native_retained_signal.tsv` | predict_crossomics_refit, native_crossomics_profiles |
| d | CapScan native correlations, all observed vs ≥10% prevalence | legend (241) | 13 participants | as b; BH within each analysis (6) | Panel d | `source_data/Fig5/native_correlations.tsv`, `distances_CAPSCAN_native_prevalence10.tsv`, `source_data/Fig5/native_distances_long.tsv` | capscan_crossomics |
| e | PREDICT controlled grid M1-M8 | legend (241) | 77 patients | outer Pearson or Spearman of distance upper triangles, 1,999 patient-label permutations; BH within each method family (10) | Panel e | `controlled_176_comparisons.tsv`, `distances_PREDICT_controlled.tsv` | crossomics_distance_sensitivity |
| f | CapScan controlled grid M1-M8, region adjusted | Hellinger vs rank, 9 of 24 (83); 9 of 24 / one of four / all 24 (421) | 13 participants | as e with whole-participant block permutations; BH within each method family (6) | Panel f | `controlled_176_comparisons.tsv`, `distances_CAPSCAN_controlled.tsv` | crossomics_distance_sensitivity |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `SEED` | 20260907 | permutation seed of all three upstream tests (Supplementary Methods 8, v28d line 815) | none |
| `N_PERMUTATIONS` | 1999 | v28d lines 453, 815, 835; smallest possible P 0.0005 | none |
| `P_TOLERANCE_NATIVE` | 1e-12 | native and controlled tests count a permuted value as extreme when \|null\| >= \|rho\| - 1e-12; functional tests use no tolerance (code only) | none |
| `FUNCTIONAL_CHUNK` | 200 | permutation chunking in the functional test; no effect on P | none |
| `Q_CUTOFF` | 0.05 | asterisks (v28d line 453) | none |
| `FUNCTIONAL_PREVALENCE` | 0.0 | all observed features in a and b (v30 paragraph 241) | ≥10% prevalence: Source Data rows with prevalence 0.1 |
| `PREDICT_SCOPE` | displayed77 | 77 patients with all layers (v28d line 453) | none |
| `CAPSCAN_SCOPE` | matched12 | producer label of the DNA-matched CapScan family, 13 participants at 347 acquisitions (v28d line 833) | 15-participant non-DNA family (`displayed_mixed12_14`): Source Data, v28d line 837 |
| `CAPSCAN_DONORS` | 1, 2, 3, 4, 6, 8-15 | participants with DNA in all five sampling groups (upstream COMPLETE.json) | none |
| `REGIONS` | duodenum, jejunum, ileum, colon, stool | producer labels of sampling regions 1-4 and stool; order defines the permutation blocks | none |
| `FIG_MM`, `FONT_PT`, `MATRIX_CELL_MM`, `PNG_DPI` | 180 x 190 mm, 6 pt, 7.6 mm, 300 dpi | v30 build: authored at the placed size in fl_style | none |

## Known limitations

- The distance matrices are recomputed into statistics here; the profiles behind them (gene TPM,
  protein-group intensities, metabolite intensities) are not in this repository. The upstream
  analyses that built the profiles and distances (A82, A84, A85, A86) are not in `pipelines/`.
- No interval is computed for any correlation; the upstream code computes only coefficients,
  permutation P values and BH q values.
- Pairwise distances are not independent observations; inference relies on permuting whole
  patients or whole five-group participant blocks.
- In CapScan, DNA, protein and metabolite cells come from the same participant and sampling group,
  not necessarily from the same physical specimen (v28d line 833).
- The v30 build asserted that the removed native, OG and KO matrices equal Fig. 4i and Fig. 5g
  against v29 build tables that are not in this repository; that duplicate test is not repeated
  here (OPEN_ITEMS.md).
