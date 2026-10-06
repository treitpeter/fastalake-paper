# Source data: Figure 6

Tables read by `figures/Fig6/Fig6.ipynb` and `figures/Fig6/panels.py`. Roster: CapScan 347
acquisitions from 15 participants (fixed roster, participant 15 included); PREDICT 93
patients (unchanged). `figures/Fig6/prepare_source_data.py` documents how each file was
copied or derived; `figures/Fig6/PROVENANCE.md` gives the origin of each file and its
comparison with the submission Source Data package.

Common identifiers: `participant` is the CapScan study participant number (1-15, not a
clinical identifier); `patient` is a PREDICT pseudonym (`P001`-`P093`). Regions are
`duodenum`, `jejunum`, `ileum`, `colon` (printed as regions 1-4) and `stool`.
Origins: `All`, `Human`, `Microbial` (conservative rule), plus `Mixed`, `Unknown` and
`NotInSelectedFASTA`, which are retained here but not drawn.

## capscan_terminal_participant_region.tsv
sha256: dd5097208a5d37e163a77d28e6ead4e9642b9db0b7fc12a5a61de03760d920e1
rows: 44532 (one row per route x unit x origin x terminus x metric x participant x region, peptide length 9-50)
produced by: terminal_participant_inference (fixed 347-acquisition roster)

| column | unit | meaning |
|---|---|---|
| cohort | id | always CAPSCAN |
| route | id | `denovo` (top 30% de novo), `semi` and `unspecific` (Sage semi-tryptic and non-specific, peptide q <= 0.01), `tryptic` (conventional search), `denovo_highconfidence*` (selection sensitivities) |
| min_length | residues | minimum peptide length; 9 for every row here (maximum 50 applied upstream) |
| unit | id | `canonical_peptide` (distinct I/L-equivalent sequences) or `accepted_psm_or_prediction_row` |
| origin | id | peptide origin class |
| terminus | id | `C_last` (C-terminal residue), `N_first` (first residue), `N_cleavage_P1` (residue before the N terminus) |
| metric | id | `nonKR_fraction_all` (non-K/R termini / all), `H_fraction_all` (H termini / all), `H_fraction_nonKR` (H / non-K/R termini) |
| participant | id | CapScan participant |
| region | id | sampling region |
| value | fraction (0-1) | metric averaged within acquisition, then specimen, then participant-region |
| specimens | count of specimens | specimens of this participant-region |
| defined_specimens | count of specimens | specimens with a defined value |

## capscan_residue_participant_region.tsv
sha256: 4ee46ef1cdab6e9d94ae0a61324c650f0522af70fc0ac2c3dac73543bd8fdc2e
rows: 163284 (as above, unit `assigned_bare_sequence` only)
produced by: terminal_participant_inference

Columns as `capscan_terminal_participant_region.tsv`. `unit` is `assigned_bare_sequence`
(assigned unmodified strings, I and L distinct). Additional `metric` values:
`residue_<X>_fraction_nonKR` (share of residue X among non-K/R termini, fraction 0-1, for
the 18 non-K/R residues) and `other_residue_fraction_nonKR`.

## capscan_terminal_statistical_families.tsv
sha256: 919360392724871c2cef8e99e5a0172cdb7d8ad8a364a476778aac7c8f10f70b
rows: 22896 (one row per planned test; all A24 families, including lengths 0 and 5 and the selection-threshold families)
produced by: terminal_participant_inference

| column | unit | meaning |
|---|---|---|
| cohort, route, min_length, unit, origin, terminus, metric | id | as above |
| contrast | id | `<a>__minus__<b>`; `proximal` = mean of regions 1-2, `distal` = mean of regions 3-4, `capsule` = mean of proximal and distal; selection families use their own labels |
| family | id | Benjamini-Hochberg family (Fig. 6 uses `terminal_fraction__length_9`, 4,536 tests, and `all18_nonKR_composition__length_9`, 13,608 tests) |
| n_participants | count of participants | paired participants with both values |
| mean_difference, median_difference | fraction (0-1) | of the paired participant differences |
| n_positive, n_negative | count of participants | participants with a positive or negative difference |
| p | dimensionless | two-sided sign-flip P (exact for n <= 16); 1 when n < 8 |
| estimable | boolean | n >= 8 |
| CI95_low, CI95_high | fraction (0-1) | 95% percentile participant-bootstrap interval of the mean difference |
| test | id | `exact_sign_flip` |
| region | id | region label of the selection families; empty otherwise |
| q | dimensionless | Benjamini-Hochberg q within `family` |

## capscan_acquisition_design.tsv
sha256: 4638d7894cdd3423e33a7416bfe2de11d9e761e6bbde4749d29e0ed4a03c7fc4
rows: 347 (one row per CapScan acquisition)
produced by: terminal_participant_inference

| column | unit | meaning |
|---|---|---|
| cohort | id | CAPSCAN |
| sample | id | raw-file name of the acquisition |
| participant | id | CapScan participant |
| region | id | sampling region |
| aliquot | id | `main` or `pellet` |
| specimen_id | id | canonical specimen (acquisitions of one specimen are averaged first) |

## capscan_histidine_statistical_families.tsv
sha256: 3c24a42eec725e0b43bdeecca3928f971925c25395dde7478d92a1db9e1c28d7
rows: 44010 (one row per planned test in four families: `source_P1_all_residues` 9,000, `all_reciprocal_controls` 90, `H_spatial_all_routes` 2,520, `all_motif_positions_residues` 32,400)
produced by: terminal_histidine_inference (fixed 347-acquisition roster)

| column | unit | meaning |
|---|---|---|
| endpoint | id | `source_enrichment`, `matched_reciprocal`, `spatial_H_*`, `ordinal_fluid_H_*`, `motif` |
| family | id | Benjamini-Hochberg family |
| route | id | analysis route (`denovo_top30`, `sage_semi_bothq`, `sage_unspecific_bothq` are drawn; `_bothq` = peptide and spectrum q <= 0.01) |
| side | id | `C` or `N` boundary |
| origin | id | All, Human, Microbial |
| residue | id | P1 residue (source tests) or motif residue |
| contrast | id | region (source and reciprocal tests), region pair, or compartment (motif) |
| n | count of participants | participants with a defined value |
| mean | fraction (0-1) | mean participant value (observed minus expected for source tests; histidine minus matched control for reciprocal tests) |
| p | dimensionless | two-sided sign-flip P; 1 when n < 8 |
| estimable | boolean | n >= 8 |
| CI95_low, CI95_high | fraction (0-1) | 95% percentile bootstrap interval of the mean, 20,000 whole-participant resamples |
| n_positive, n_negative | count of participants | sign counts |
| observed_fraction_mean, expected_fraction_mean | fraction (0-1) | mean observed residue fraction at boundaries and mean source-protein expectation |
| enrichment_ratio | dimensionless | observed_fraction_mean / expected_fraction_mean |
| ratio_bootstrap_defined | count of resamples | bootstrap resamples with a defined ratio |
| enrichment_CI95_low, enrichment_CI95_high | dimensionless | 95% bootstrap interval of the ratio |
| population, position | id, residues | motif population and position relative to the boundary (motif tests) |
| q | dimensionless | Benjamini-Hochberg q within `family` |

## capscan_histidine_primary_source_enrichment.tsv
sha256: aa6587cd2d51c19a1ed0d38d510708d420ade2552387435a4d7a31f432f9845b
rows: 90 (source-enrichment rows of `capscan_histidine_statistical_families.tsv` with residue H and the three drawn routes)
produced by: `panels.derive_histidine_inputs` from terminal_histidine_inference

Columns as `capscan_histidine_statistical_families.tsv`. Panels d and e draw the `origin == All` rows.

## capscan_histidine_primary_reciprocal_tests.tsv
sha256: 6b12f31b2e7360dd3675ab37c6e2e2c33689ce97c582fc84728d5ca356fe4109
rows: 45 (matched-reciprocal rows for the three drawn routes)
produced by: `panels.derive_histidine_inputs` from terminal_histidine_inference

Columns as above; panel f draws the `origin == All` rows.

## capscan_histidine_participant_region_enrichment.tsv
sha256: a830d42712865d23f6642060c782b2dd00b49db19b551dc3575cdb359f619722
rows: 106160 (one row per route x side x origin x residue x participant x region)
produced by: terminal_histidine_inference

| column | unit | meaning |
|---|---|---|
| route, side, origin, residue, cohort, participant, region | id | as above |
| O | fraction (0-1) | observed fraction of uniquely localized internal non-K/R boundaries with this P1 residue, averaged within specimen then participant-region |
| E | fraction (0-1) | source-protein expectation over eligible internal non-K/R positions, averaged the same way |
| difference | fraction (0-1) | O - E |
| enrichment | dimensionless | O / E |

## capscan_histidine_participant_reciprocal.tsv
sha256: 74613f6da70994b5307cf95660f40f1eb3a22cc9757175c48d5b28161980263c
rows: 1350 (one row per route x origin x participant x region)
produced by: terminal_histidine_inference

| column | unit | meaning |
|---|---|---|
| route, origin, cohort, participant, region | id | as above |
| difference | fraction (0-1) | complementary-product rate at histidine boundaries minus matched non-histidine boundaries |

## predict_histidine_primary_source_enrichment.tsv
sha256: 8f2cb3185f336832a59d0ba0be519e76eda1134275ff8febb98c604c887de3cf
rows: 12 (residue H, routes `denovo_top30` and `sage_unspecific_bothq`, two sides, three origins)
produced by: terminal_histidine_inference (PREDICT run), selected by the end-of-day histidine review

Columns as `capscan_histidine_statistical_families.tsv`; `contrast` is
`ALL_observed_specimens`, the unit is the patient (n = 93; specimen means averaged within
patient), sign-flip P from 200,000 random sign vectors with plus-one correction.

## predict_histidine_primary_reciprocal_tests.tsv
sha256: 61e8490827b12dcd29720e3750706e4900184b2e88fda92ecdafd243db258ee8
rows: 12 (all PREDICT matched-reciprocal tests, four routes x three origins)
produced by: terminal_histidine_inference (PREDICT run)

Columns as above.

## predict_histidine_statistical_families.tsv
sha256: 85c774cb604a1d9e6c9526b04750b10d6601a07855fe1552357f66be016d28d7
rows: 1572 (complete PREDICT families: `source_P1_all_residues` 1,560, `all_reciprocal_controls` 12)
produced by: terminal_histidine_inference (PREDICT run)

Columns as `capscan_histidine_statistical_families.tsv` (no motif columns). No patient-level values.

## predict_histidine_patient_enrichment.tsv
sha256: aec3d14f13c53e75d177375300e5e94c4f777fb97d5cabaffb60f5620e7b9f0e
rows: 1116 (one row per route x side x origin x patient, residue H)
produced by: terminal_histidine_inference (PREDICT run), patient means by the end-of-day histidine review; raw patient identifiers replaced by pseudonyms

| column | unit | meaning |
|---|---|---|
| route, side, origin, residue | id | as above |
| patient | id | pseudonym P001-P093, assigned in ascending order of the raw identifier (the order used by the original bootstrap); the map is a controlled input |
| O | fraction (0-1) | observed histidine fraction at internal boundaries, specimen means averaged within patient |
| E | fraction (0-1) | source-protein expectation, averaged the same way |
| difference | fraction (0-1) | O - E |
