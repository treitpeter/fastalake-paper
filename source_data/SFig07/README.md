# Source data: Supplementary Figure 7

All tables are tab-separated with a header row. No clinical covariates are included. PREDICT
participants carry pseudonymised identifiers (`PREDICT_P001` ...); CapScan participants carry
study numbers 1-15. Acquisition names are the public raw-file names of the healthy-stool study.
Origins, hashes and agreement with the Source Data package are in
`figures/SFig07/PROVENANCE.md`.

Arm labels (see `arm_definitions.tsv`): `original` acquisition-specific searches;
`union_individual` each acquisition searched separately against the common 153,267-protein
FASTA; `joint_corrected` joint search of all 58 acquisitions (corrected Sage build), peptide q
estimated once across the study (study-level confidence); `joint_original` the same joint
search with the uncorrected Sage build; `joint_fileq` joint-search PSMs with peptide q
re-estimated within each acquisition (per-acquisition confidence); `joint_sageq_repro`
`joint_corrected` re-extracted with study-level q (positive control); `original_fileq` the
acquisition-specific searches re-scored with the per-acquisition estimator (sensitivity);
`repeat_individual` the acquisition-specific searches repeated. Tiers: `TechnicalReplicates`
(12 acquisitions), `WorkflowReplicates` (15 independently prepared samples), `OnePiece`,
`Intra`, `Interindividual`.

## identity_comparison.tsv (panel a)

One row per acquisition, arm and minimum peptide length (58 x 4 arms x 2 lengths = 464 rows).

| column | meaning |
|---|---|
| sample | acquisition |
| tier | experiment tier |
| test | arm compared with `original` |
| reference | always `original` |
| min_length | 5 (5-50 residues) or 9 (9-50 residues) |
| original | accepted canonical peptides in the original search |
| test_count | accepted canonical peptides in the test arm |
| overlap | peptides accepted in both |
| gained, lost | peptides only in the test arm, only in the original |
| count_ratio | test_count / original |
| original_identity_recovery | overlap / original (fraction) |

## identity_per_acquisition.tsv (panel a, Results)

One row per arm and acquisition (6 arms x 58), peptides of 9-50 residues (the accepted
5-50-residue identities with at least 9 residues; recount `hstool_recount.py`, see
PROVENANCE.md).

| column | meaning |
|---|---|
| arm | arm label |
| sample, tier | acquisition and tier |
| original | accepted 9-50-residue peptides in the original search |
| test | accepted 9-50-residue peptides in the arm |
| overlap | accepted in both |
| count_ratio | test / original |
| original_identity_recovery | overlap / original (fraction) |

## identity_summary_by_length.tsv (Results, Methods)

One row per minimum length (5 or 9) and arm, written by the 9-50 recount.

| column | meaning |
|---|---|
| min_length | 5 (5-50 residues) or 9 (9-50 residues) |
| arm | arm label |
| median_count_ratio | median over 58 acquisitions of count_ratio |
| n_higher | acquisitions with count_ratio > 1 |
| median_identity_recovery | median over 58 acquisitions of original_identity_recovery |
| study_distinct | distinct accepted canonical peptides across all 58 acquisitions |
| ids_in_all12_technical | peptides accepted in all 12 technical replicates |

## completeness_by_length.tsv (Results)

| column | meaning |
|---|---|
| min_length | 5 or 9; features kept if the canonical peptide has at least this many residues |
| space | `features` (peptide features, modifications kept, charges combined) |
| arm, tier | arm and tier |
| acquisitions | acquisitions in the tier (12 or 15) |
| rows_complete_all | features positive with an accepted same-acquisition identification in every acquisition of the tier |

## common_search_participants.tsv (panel b, Methods)

One row per participant (93 PREDICT, 15 CapScan). Peptide columns count peptides of 9-50
residues (recount `clin_recount.py`); `gained_peptides` and `lost_peptides` were not
recounted at 9-50 and are empty. Feature and protein-group columns are unchanged.
Acquisition values were averaged within specimen, then specimens within participant. Prefix `original_` is the participant's own
FASTA, `common_` the common study FASTA.

| column | meaning |
|---|---|
| participant | pseudonymised PREDICT identifier or CapScan participant number |
| original_peptides, common_peptides | accepted peptides |
| gained_peptides, lost_peptides | peptides only under common, only under original |
| original_positive_features, common_positive_features, gained_positive_features, lost_positive_features | positive quantified features, same pattern |
| original_sum_groups, common_sum_groups | quantified protein groups by summed intensity |
| original_directlfq_groups, common_directlfq_groups | quantified protein groups after directLFQ |
| original_directlfq_group_cells_lost, common_directlfq_group_cells_lost | sum groups minus directLFQ groups |
| `<metric>_difference` | common minus original |
| `<metric>_percent_difference` | 100 x (common / original - 1); plotted: peptides (x), directlfq_groups (y) |
| cohort | `PREDICT` or `CAPSCAN` |

## microbial_group_pairs.tsv (panels c, d, line 769)

One row per definition, comparison, method, tier and acquisition pair.

| column | meaning |
|---|---|
| definition | `own_rows` (positive in both acquisitions under the arm scored) or `matched` (positive under both arms in both acquisitions) |
| comparison | arm (own_rows) or `<reference>_vs_<test>` (matched) |
| method | arm whose intensities are scored |
| tier | `TechnicalReplicates` (66 pairs) or `WorkflowReplicates` (105 pairs) |
| sample_a, sample_b | the two acquisitions |
| common_rows | microbial protein-group rows used |
| log2_pearson | Pearson r of log2 intensities over those rows |
| median_pair_cv_pct | median over rows of 100 x population SD / mean of the two raw intensities |

## microbial_group_pair_summary.tsv (line 769)

Upstream median over pairs of `microbial_group_pairs.tsv`: columns definition, comparison,
method, tier, pairs, common_rows (median), log2_pearson (median), median_pair_cv_pct (median).

## common_search_pair_CV.tsv (producer consistency check)

Earlier pair table for study-level arms, used by the producer to check that the per-pair CVs
of the study-level columns are unchanged. Columns: protein_class, arm, definition
(`A_arm_alone` = own rows, `B_matched_both_arms` = matched), tier, sample_a, sample_b,
n_rows, median_pair_cv_pct, log2_pearson.

## technical_agreement_summary.tsv (line 769)

Function-averaged pairwise agreement of the leading annotated species, study-level confidence.

| column | meaning |
|---|---|
| method | `fixed_<arm>` (fixed annotation dictionary) or `native_*` |
| tier | `TechnicalReplicates` or `WorkflowReplicates` |
| threshold | `both_resolved` (leading species resolved in both acquisitions; printed values) or `both_strict` |
| original_agreement, test_agreement | agreement (fraction) under the original and the test arm |

## arm_definitions.tsv

Two columns, `arm` and `definition`: the arm labels used above.
