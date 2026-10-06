# Source data: Supplementary Figure 11

All files were written by `figures/SFig11/prepare_source_data.py` from the original producer
inputs and the distance matrices saved by the upstream analyses. No patient or participant
identifier other than the CapScan participant number appears; PREDICT patients are numbered by
their position in the 77-patient analysis roster. No clinical covariate is present.

The figure also reads four tables written by the Figure 5 unit (not duplicated here):
`source_data/Fig5/functional_correlations.tsv` (CapScan rows of the producer's functional input,
byte-identical), `source_data/Fig5/native_correlations.tsv` (CapScan rows of the producer's native
input, byte-identical), `source_data/Fig5/functional_distances_long.tsv` (CapScan functional
distances, scope `matched12`, prevalence 0, used for the recomputation) and
`source_data/Fig5/native_distances_long.tsv` (CapScan native distances, all observed features).
`prepare_source_data.py` asserts that these equal the upstream outputs exactly.

Correlation tables use the producer's column names. Layer names in functional tables: `metaG`
(DNA), `metaT` (RNA), `microbial_metaP` (protein excluding human-assigned signal), `host_metaP`
(human-assigned protein), `metabolites`. Native tables: `DNA`, `RNA`, `Non-host protein`,
`Host protein`, `Metabolites`.

## functional_correlations_PREDICT.tsv
sha256: e4719d73c2629b400a39ccf9adc3e3d678790bde413845a654408b8dfc9a0f00
rows: 56 (one row per PREDICT functional-profile comparison; 28 per prevalence setting)
produced by: predict_crossomics_refit

| column | unit | meaning |
|---|---|---|
| cohort | id | `PREDICT` |
| scope | id | `displayed77`: the 77-patient population |
| level | id | `OG`, `KO` or `Pathway` |
| layer_a, layer_b | id | the two molecular layers compared |
| prevalence | fraction | 0 = all observed features (plotted); 0.1 = features positive in at least 10% of profiles (sensitivity, not plotted) |
| n_units | count of patients | profiles in the comparison |
| n_donors | count of patients | patients in the comparison |
| rho | dimensionless | Spearman correlation between the upper triangles of the two Euclidean distance matrices |
| p | probability | two-sided permutation P, (1 + extreme) / 2,000, 1,999 permutations |
| q | probability | Benjamini-Hochberg q within the 28-comparison family of the same prevalence setting |

## native_correlations_PREDICT.tsv
sha256: c7b9a1b7291510bbfb4788eb5c59062cb30a0c696b158d49d2d32dd1bbc6436e
rows: 20 (one row per PREDICT native-profile comparison; 10 per analysis)
produced by: predict_crossomics_refit

| column | unit | meaning |
|---|---|---|
| cohort | id | `PREDICT` |
| analysis | id | `sensitivity_all_features` (all observed features; plotted in a and c) or `primary_prevalence10` (at least 10% prevalence; plotted in c) |
| layer_a, layer_b | id | the two layers compared |
| rho | dimensionless | Spearman correlation between distance-matrix upper triangles |
| p | probability | two-sided permutation P |
| permutations | count | 1,999 |
| n_profiles | count of patients | 77 |
| n_participants | count of patients | 77 |
| q | probability | Benjamini-Hochberg q within the 10 comparisons of the analysis |

## controlled_176_comparisons.tsv
sha256: 1705248887b442a01f3c52b36c584b3107db38d5f8e5c53067c32081d506ffca
rows: 176 (one row per controlled distance comparison: 80 PREDICT, 96 CapScan)
produced by: crossomics_distance_sensitivity

| column | unit | meaning |
|---|---|---|
| cohort | id | `PREDICT` or `CAPSCAN` |
| features | id | `prevalence10` (M1-M4) or `all_observed` (M5-M8) |
| distance | id | `hellinger` (Euclidean distance of Hellinger-transformed profiles) or `rank` (rank-profile distance) |
| adjustment | id | `unadjusted` or `region_adjusted` (CapScan only; sampling-group means removed) |
| outer_correlation | id | `pearson` or `spearman` correlation between distance-matrix upper triangles |
| layer_a, layer_b | id | layers compared |
| coefficient | dimensionless | outer correlation |
| p | probability | two-sided permutation P, 1,999 permutations |
| p_positive_tail | probability | one-sided (greater) permutation P |
| n_profiles | count of profiles | 77 (PREDICT) or 65 (CapScan, 13 participants x 5 sampling groups) |
| n_participants | count of patients or participants | 77 or 13 |
| q | probability | Benjamini-Hochberg q of `p` within the family of the same cohort, features, distance, adjustment and outer correlation (10 PREDICT or 6 CapScan comparisons) |

## distances_PREDICT_functional.tsv
sha256: 8f1d903f2fc48d18d77ee8d26a6537bb8684f51353f0d5912d4711a38446c72f
rows: 43,890 (one row per patient pair per matrix; 15 matrices x 2,926 pairs)
produced by: predict_crossomics_refit

| column | unit | meaning |
|---|---|---|
| level | id | `OG`, `KO` or `Pathway` |
| layer | id | functional layer name (see above) |
| profile_i, profile_j | id | 1-based positions of the two patients in the 77-patient roster, i < j, row-major |
| distance | dimensionless | Euclidean distance between the transformed profiles (all observed features) |

## distances_PREDICT_native.tsv
sha256: b1adee7ce195ba29d983c7a5a893a40d5e75af497cbc3d6ace486c448350122f
rows: 29,260 (10 matrices x 2,926 pairs)
produced by: predict_crossomics_refit

| column | unit | meaning |
|---|---|---|
| analysis | id | `primary_prevalence10` or `sensitivity_all_features` |
| layer | id | native layer name |
| profile_i, profile_j | id | as above |
| distance | dimensionless | Euclidean (Hellinger) distance; metabolites: Euclidean distance of standardized log10 profiles |

## distances_PREDICT_controlled.tsv
sha256: a6d4b5ce0fe5da490788d80f51cd2ca185fc10a20a266d68912641a0ad20e9e3
rows: 58,520 (20 matrices x 2,926 pairs)
produced by: crossomics_distance_sensitivity

| column | unit | meaning |
|---|---|---|
| features | id | `all_observed` or `prevalence10` |
| distance_kind | id | `hellinger` or `rank` |
| adjustment | id | `unadjusted` |
| layer | id | native layer name |
| profile_i, profile_j | id | as above |
| distance | dimensionless | distance of that kind |

## distances_CAPSCAN_native_prevalence10.tsv
sha256: 22fb297d74d1e7bf542258feee55d960e7bf8115c325253c9e148984ceed7bb2
rows: 8,320 (4 matrices x 2,080 pairs)
produced by: capscan_crossomics

| column | unit | meaning |
|---|---|---|
| analysis | id | `primary_prevalence10` |
| layer | id | native layer name |
| participant_i, participant_j | id | CapScan participant numbers (13 DNA-complete participants) |
| region_i, region_j | id | sampling group (`duodenum`, `jejunum`, `ileum`, `colon` are the producer labels for intestinal sampling regions 1-4; `stool`) |
| distance | dimensionless | Euclidean distance of sampling-group-adjusted profiles |

## distances_CAPSCAN_controlled.tsv
sha256: d4e18fe35ebd39650d86971f9a83f36d1d7295b9a3feb9558fb6f5a60c48e217
rows: 66,560 (32 matrices x 2,080 pairs)
produced by: crossomics_distance_sensitivity

| column | unit | meaning |
|---|---|---|
| features, distance_kind | id | as for PREDICT |
| adjustment | id | `unadjusted` or `region_adjusted` |
| layer | id | native layer name |
| participant_i, region_i, participant_j, region_j | id | as above |
| distance | dimensionless | distance of that kind |

## predict_native_retained_signal.tsv
sha256: c87ff05fe7eccec989fe5a0c08c32d0fb26e821829d38af348b91c370dec2848
rows: 154 (one row per patient and layer; 77 DNA, 77 RNA)
produced by: native_crossomics_profiles

| column | unit | meaning |
|---|---|---|
| cohort | id | `PREDICT` |
| analysis | id | `primary_prevalence10` |
| layer | id | `DNA` or `RNA` |
| profile_index | id | 1-based row position within the layer (original order; the original patient label is removed) |
| retained_signal_fraction | fraction | share of sequence-resolved TPM kept by the at-least-10%-prevalence filter |

## Correspondence with the Source Data package

| this file | package sheet | agreement |
|---|---|---|
| functional_correlations_PREDICT.tsv | `Biology_refit/cross_refit/functional_correlations.tsv`, PREDICT rows | identical (the package sheet is the same file, sha256 4c05b200...) |
| native_correlations_PREDICT.tsv | `Biology_refit/cross_refit/native_correlations.tsv`, PREDICT rows | identical |
| CapScan rows (source_data/Fig5) | same two package sheets, CAPSCAN rows | differ: the package still holds the 12-participant (331-acquisition) CapScan rows |
| controlled_176_comparisons.tsv | `Biology_refit/distance_refit/controlled_176_comparisons.tsv` | PREDICT rows identical; the 96 CapScan rows differ (package holds 12 participants) |
| distances_* | `Biology_refit/cross_refit/*_distances.npz`, `Biology_refit/distance_refit/*_distances.npz` | PREDICT matrices identical; CapScan matrices in the package are the 12-participant ones |
| predict_native_retained_signal.tsv | no package sheet; from the Crossomics companion `native/SOURCE_DATA_retained_signal_by_profile.tsv` | values identical, identifier column removed |
