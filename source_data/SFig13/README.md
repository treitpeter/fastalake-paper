# Source data: Supplementary Figure 13

Tables read by `figures/SFig13/SFig13.ipynb`. The CapScan complete family table is shared with
Figure 6 and lives at `source_data/Fig6/capscan_histidine_statistical_families.tsv` (sha256
`3c24a42eec725e0b43bdeecca3928f971925c25395dde7478d92a1db9e1c28d7`); it is not duplicated here.
How each file was made: `figures/SFig13/prepare_source_data.py`.

PREDICT identifiers: raw patient identifiers, visit codes, specimen identifiers and raw file
names are replaced by pseudonyms (`P001`..`P093`, `P###_S##`, `P###_S##_A##`, `U##`). The patient
pseudonyms are the same as in `source_data/Fig6/`. No clinical covariate is in any file.

## capscan_histidine_motif_display.tsv
sha256: 269fa990cd52029a6dd4961e7bf21e347436c374c89bee4fa5c9631dd87364ea
rows: 840 (one row per displayed motif cell: compartment x boundary side x flank position x residue)
produced by: terminal_histidine_inference (row selection by `figures/SFig13/panels.py:derive_display_table`)

Rows of the CapScan family table with endpoint `motif`, population `H_only`, route `denovo_top30`,
origin `All`, position not -1. Byte-identical to the v28f replacement Source Data file
`Supplementary_Figure_13/CapScan_histidine_motif_display.tsv`.

| column | unit | meaning |
|---|---|---|
| endpoint | id | `motif` |
| family | id | correction family, `all_motif_positions_residues` (32,400 tests) |
| route | id | evidence route, `denovo_top30` (de novo top-30% predictions) |
| side | id | boundary side, `N` or `C` |
| origin | id | peptide origin class, `All` |
| residue | id | amino acid at the flanking position |
| contrast | id | compartment: `proximal` (duodenum and jejunum), `distal` (ileum and colon), `stool` |
| n | count of participants | participants with a value in every region of the compartment |
| mean | fraction | participant mean of observed minus source-expected residue frequency (x100 = percentage points, the colour) |
| p | dimensionless | two-sided sign-flip P (exact enumeration for n <= 16) |
| estimable | boolean | n >= 8 |
| CI95_low, CI95_high | fraction | 95% percentile bootstrap interval of the mean (20,000 participant resamples) |
| n_positive, n_negative | count of participants | participants with a positive or negative value |
| observed_fraction_mean, expected_fraction_mean, enrichment_ratio, ratio_bootstrap_defined, enrichment_CI95_low, enrichment_CI95_high | | empty for motif rows (used by source-enrichment rows only) |
| population | id | `H_only`: boundaries whose P1 residue is histidine |
| position | residue positions | flank position relative to the boundary (-4..-2, +1..+4) |
| q | dimensionless | Benjamini-Hochberg q within the 32,400-test family; dots at q < 0.05 |

## capscan_participant_motif_differences.tsv.gz
sha256: df84cc2720b0ed63448048ef0676ad33c44ef9fe46eb05414a1d98a3ecb5ae77
rows: 864,000 (one row per route x side x origin x population x position x residue x participant x region)
produced by: terminal_histidine_inference (`participant_motif_differences.tsv.gz`, 347 CapScan acquisitions, 15 participants)

Input to the recomputation of the complete motif family (`panels.stat_motif_family`).

| column | unit | meaning |
|---|---|---|
| route | id | evidence route (six routes) |
| side | id | `N` or `C` boundary |
| origin | id | `All`, `Human` or `Microbial` |
| population | id | `all_nonKR` (all non-K/R boundaries) or `H_only` (P1 = histidine) |
| position | residue positions | flank position -4..-1, +1..+4 |
| residue | id | amino acid |
| cohort | id | `CAPSCAN` |
| participant | id | CapScan participant number 1-15 (study numbering, no clinical data) |
| region | id | duodenum, jejunum, ileum, colon or stool |
| difference | fraction | observed minus source-expected residue frequency, averaged over acquisitions within specimen and specimens within participant and region |

## predict_histidine_acquisition_enrichment.tsv.gz
sha256: b40ea94ad1ba55f9931bce83035728cd7429f4d4c5d20733c08733beb3bfb951
rows: 344,760 (one row per acquisition x route x side x origin x population x residue; 306 acquisitions, 302 with a patient)
produced by: terminal_histidine_inference (`acquisition_enrichment.tsv.gz` of the PREDICT run), pseudonymised by `figures/SFig13/prepare_source_data.py`

Input to panels g,h (`panels.predict_patient_means`) and to the recomputation of the PREDICT
source family (`panels.stat_predict_source_family`). The analysis uses `O = observed / denominator`
and `E = expected / denominator`.

| column | unit | meaning |
|---|---|---|
| route | id | evidence route (de novo score bins, `denovo_top30`, `sage_unspecific`, `sage_unspecific_bothq`, ...) |
| side | id | `N` or `C` boundary |
| origin | id | `All`, `Human` or `Microbial` |
| population | id | `all_nonKR` |
| residue | id | P1 amino acid |
| observed | count of boundaries | observed internal non-K/R boundaries with this P1 residue |
| expected | count of boundaries | sum over the same boundaries of the source-protein expected frequency of this residue |
| denominator | count of boundaries | observed internal non-K/R boundaries |
| sites | count of boundaries | as written by the A28 per-acquisition step; not used here |
| undefined_background | count of boundaries | as written by the A28 per-acquisition step; not used here |
| enrichment | dimensionless | observed / expected; not used here |
| cohort | id | `PREDICT` |
| sample | id | acquisition pseudonym `P###_S##_A##` (or `U##` without a patient) |
| participant | id | patient pseudonym `P001`..`P093`; empty for 4 acquisitions without a patient (excluded) |
| region | id | `withheld` (visit code withheld); empty without a patient |
| aliquot | id | `stool` |
| specimen_id | id | specimen pseudonym `P###_S##` (300 specimens) |

## predict_histidine_statistical_families.tsv
sha256: 85c774cb604a1d9e6c9526b04750b10d6601a07855fe1552357f66be016d28d7
rows: 1,572 (1,560 source-enrichment tests and 12 complementary-product tests)
produced by: terminal_histidine_inference (`complete_statistical_families.tsv` of the PREDICT run; byte copy)

Source of the plotted means and 95% intervals in panels g,h; compared with the recomputation.

| column | unit | meaning |
|---|---|---|
| endpoint | id | `source_enrichment` or `matched_reciprocal` |
| family | id | `source_P1_all_residues` (1,560) or `all_reciprocal_controls` (12) |
| route, side, origin, residue | id | as above (side and residue empty for reciprocal rows) |
| contrast | id | `ALL_observed_specimens` (patient mean over all specimens) |
| n | count of patients | patients with a value |
| mean | fraction | patient mean of observed minus expected (x100 = percentage points) |
| p | dimensionless | two-sided sign-flip P (200,000 random assignments, plus-one) |
| estimable | boolean | n >= 8 |
| CI95_low, CI95_high | fraction | 95% percentile bootstrap interval (20,000 patient resamples) |
| n_positive, n_negative | count of patients | patients with a positive or negative value |
| observed_fraction_mean, expected_fraction_mean | fraction | patient means of O and E |
| enrichment_ratio | dimensionless | ratio of the two means |
| ratio_bootstrap_defined | count of resamples | bootstrap resamples with a defined ratio |
| enrichment_CI95_low, enrichment_CI95_high | dimensionless | 95% bootstrap interval of the ratio |
| q | dimensionless | Benjamini-Hochberg q within the family |
