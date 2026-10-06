# Source data: Supplementary Table 1

CAMPI candidate-clustering sensitivity. Acquisition ids `S01`..`S12` and
`S11_Fraction1`..`S11_Fraction4` are the 15 CAMPI acquisitions. Arms: `unclustered`, and
`cluster_0XX` for MMseqs2/Linclust at XX% minimum sequence identity (99, 97, 95, 90, 75,
50). All ratios, losses and gains are relative to the unclustered search of the same
acquisition. Peptides are distinct accepted canonical peptides (peptide q <= 0.01, 5-50
residues, I and L merged). No clinical or personal data.

The two `per_acquisition_*` tables are the inputs of the table. The `check_*` and
`independent_recount_*` tables are used only by the consistency checks in
`tables/STab1/panels.py:checks`. The `support_97_99_*` tables are the eight earlier
97%/99%-only analyses that the Source Data package carries beside Supplementary Table 1
(`Supplementary_Tables/Table1_support/`); they are kept here unchanged and are not used to
compute any printed number.

## per_acquisition_identifications.tsv
sha256: 0187448fe338050276ab9b8748e15ef820e1888ba62019556e08487788511c06
rows: 105 (one row per acquisition and arm)
produced by: campi_candidate_clustering (collector `collect.py`)

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition |
| arm | id | `unclustered` or `cluster_0XX` |
| identity | dimensionless | nominal Linclust --min-seq-id (1.0 for unclustered) |
| protein_entries | count of protein entries | entries in the searched FASTA |
| protein_entry_ratio | dimensionless | protein_entries / unclustered protein_entries |
| peptides | count of peptides | accepted peptides in this arm |
| baseline_peptides | count of peptides | accepted peptides in the unclustered arm |
| peptide_count_ratio | dimensionless | peptides / baseline_peptides |
| peptide_reduction_pct | % | 100 x (1 - peptide_count_ratio) |
| baseline_peptide_identities_retained | count of peptides | unclustered peptides also accepted in this arm |
| baseline_identity_recovery | fraction | baseline_peptide_identities_retained / baseline_peptides |
| REF_peptides | count of peptides | peptides accepted in the reference-database search of the acquisition |
| REF_recovered | count of peptides | REF_peptides also accepted in this arm |
| REF_recovery | fraction | REF_recovered / REF_peptides |
| lost | count of peptides | accepted in unclustered, not in this arm |
| gained | count of peptides | accepted in this arm, not in unclustered |
| lost_sequence_absent | count of peptides | lost peptides whose exact sequence is absent from this arm's searched FASTA |
| lost_sequence_available | count of peptides | lost peptides whose exact sequence is still present in this arm's searched FASTA |

## per_acquisition_KO_coverage.tsv
sha256: 4dd098b16a55e0cd208f7f6b23e5eb0e190e6cc236d8a9be2fc1c4473f6c124f
rows: 105 (one row per acquisition and arm)
produced by: campi_candidate_clustering (collector `collect.py`, one fixed peptide-to-KO mapping)

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition |
| arm | id | arm |
| compatible_KOs | count of KOs | union of KEGG Orthology ids compatible with the arm's accepted peptides |
| baseline_compatible_KOs | count of KOs | same for the unclustered arm |
| baseline_KOs_retained | count of KOs | KOs in both |
| KO_recovery | fraction | baseline_KOs_retained / baseline_compatible_KOs |
| KO_lost | count of KOs | in unclustered, not in this arm |
| KO_gained | count of KOs | in this arm, not in unclustered |
| single_KO_evidence_functions | count of KOs | KOs supported by peptides that map to exactly one KO |
| baseline_single_KOs | count of KOs | same for the unclustered arm |
| single_KO_recovery | fraction | single-KO functions retained / baseline_single_KOs |
| peptides | count of peptides | accepted peptides in this arm |
| unmapped_peptides | count of peptides | peptides not found in the fixed sequence pool |
| mapped_without_KO | count of peptides | found in the pool, no KO |
| single_KO_peptides | count of peptides | peptides compatible with exactly one KO |
| multi_KO_peptides | count of peptides | peptides compatible with more than one KO |

## independent_recount_identifications.tsv
sha256: 2aaf6ffe20d9cfadf8cc7f87b1ed4bb2c4606b2c18f46c2896fa4ed5d301de56
rows: 105 (one row per acquisition and arm)
produced by: independent recount from the Sage result files (own parser), see tables/STab1/PROVENANCE.md

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition |
| arm | id | arm |
| peptides | count of peptides | recounted accepted peptides |
| baseline_peptides | count of peptides | recounted unclustered accepted peptides |
| ratio | dimensionless | peptides / baseline_peptides |
| reduction_pct | % | 100 x (1 - ratio) |
| lost | count of peptides | in unclustered, not in this arm |
| gained | count of peptides | in this arm, not in unclustered |
| baseline_equals_saved_H08_list | boolean | recounted unclustered set equals the saved accepted list of the earlier 97/99% run |
| entries | count of protein entries | entries in the searched FASTA |
| entry_ratio | dimensionless | entries / unclustered entries |

## independent_recount_KO.tsv
sha256: 7a552d9dac7eea7bba7abd15b651dda2f970fc369f44dd5166c805cebcd1ac9e
rows: 105 (one row per acquisition and arm)
produced by: independent recount, see tables/STab1/PROVENANCE.md

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition |
| arm | id | arm |
| KOs | count of KOs | recounted compatible KOs |
| baseline_KOs | count of KOs | recounted unclustered compatible KOs |
| KO_recovery | fraction | retained baseline KOs / baseline_KOs |
| peptides_absent_from_mapping | count of peptides | accepted peptides missing from the fixed peptide-to-KO mapping (0 everywhere) |

## check_new_rows_95_90_75_50.tsv
sha256: 689dc72952f0abf0bfe9542c20171415a19921ffd99872d7ca0f6a2a09b558cc
rows: 60 (one row per acquisition and arm, arms 95/90/75/50% only)
produced by: earlier extraction of the four new arms in the 45-row table's column names

Columns as in `support_97_99_identification_comparison.tsv` below.

## support_97_99_identification_comparison.tsv
sha256: 43c8a578785dd6781f67c28bc84e5d2baf956faa6d4b675410e770192c550fe4
rows: 45 (one row per acquisition and arm; unclustered, 97%, 99%)
produced by: campi_candidate_clustering, 97/99% baseline searches (H08)

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition |
| arm | id | arm |
| accepted | count of peptides | accepted peptides |
| reference_denominator | count of peptides | peptides accepted in the reference-database search |
| reference_recovered | count of peptides | of those, accepted in this arm |
| lost_vs_unclustered | count of peptides | in unclustered, not in this arm |
| gained_vs_unclustered | count of peptides | in this arm, not in unclustered |
| count_ratio | dimensionless | accepted / unclustered accepted |
| reference_recovery_fraction | fraction | reference_recovered / reference_denominator |
| baseline_identity_retained | fraction | unclustered peptides also accepted in this arm / unclustered accepted |

## support_97_99_all_search_counts.tsv
sha256: cfdca8b39a5216a82f6c551d7a2d065921052cf07fd9b0934e33d6c360257831
rows: 45 (one row per acquisition and arm; unclustered, 97%, 99%)
produced by: campi_candidate_clustering, H08 searches

Columns `sample` to `gained_vs_unclustered` as in `support_97_99_identification_comparison.tsv`.

## support_97_99_paired_count_ratios.tsv
sha256: bd6067eb5080bc6ec77f739a408b43ba20157219b769daf030a8ee1d342bfbf3
rows: 15 (one row per acquisition)
produced by: campi_candidate_clustering, H08 searches

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition |
| cluster_097 | dimensionless | accepted peptides at 97% / unclustered |
| cluster_099 | dimensionless | accepted peptides at 99% / unclustered |

## support_97_99_configuration_and_identity_checks.tsv
sha256: e072078cc313abf21a5a2333723548bcd55e109cd345b6894fbcced5b2289993
rows: 15 (one row per acquisition)
produced by: campi_candidate_clustering, H08 review

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition |
| all_three_identification_settings_equal | boolean | unclustered, 97% and 99% searches used the same identification settings |
| accepted_identity_counts_recounted | boolean | accepted peptide counts were recounted and agree |

## support_97_99_S02_all_baseline_sensitivities.tsv
sha256: fdb568fc0de07b744f4425509fa11d5479d1243182f69155df2b5d4bafa94f1f
rows: 12 (acquisition S02; one row per baseline and arm)
produced by: campi_candidate_clustering, H08 S02 recovery

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition (S02 only) |
| baseline | id | which unclustered S02 search is the denominator: `archived`, `baseline`, `repeat_1`, `repeat_2`, `saved_contemporaneous` |
| arm | id | arm |
| baseline_accepted | count of peptides | accepted peptides in that baseline |
| accepted | count of peptides | accepted peptides in the arm |
| gained | count of peptides | in the arm, not in the baseline |
| lost | count of peptides | in the baseline, not in the arm |
| baseline_retained | fraction | baseline peptides also accepted in the arm / baseline_accepted |

## support_97_99_availability_loss_gain_attribution.tsv
sha256: 419bba008bee26a202a3177b034e97a4f7afbe9cac71c3061e26aa23e9ce85b1
rows: 120 (one row per acquisition, arm (97%, 99%) and category)
produced by: campi_candidate_clustering, H08 searches

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition |
| arm | id | `cluster_097` or `cluster_099` |
| category | id | `lost_sequence_absent_from_clustered_selection`, `lost_sequence_still_available`, `gained_sequence_absent_from_original_selection`, `gained_sequence_already_available` |
| count | count of peptides | peptides in the category |
| reference_peptides | count of peptides | of those, also accepted in the reference-database search |
| scope | text | fixed scope statement written by the producer |

## support_97_99_all_H04_availability.tsv
sha256: 50bdabf4640ffb346a3a81338f03ce84075bfde5a563f40eb050bf337de750c7
rows: 180 (one row per acquisition, peptide universe and selection stage)
produced by: campi_candidate_clustering, H04 retention (exact substring availability, no search)

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition |
| universe | id | fixed peptide set queried: `FL` (FastaLake search) or `REF` (reference-database search) |
| stage | id | sequence set: candidates, cluster representatives or selected sequences, unclustered or at 97/99% |
| n_proteins | count of protein entries | sequences in that stage |
| n_query_peptides | count of peptides | peptides in the universe |
| available | count of peptides | universe peptides present as exact substrings (I=L) in the stage |
| lost_relative_to_candidates | count of peptides | available in the unclustered candidates but not in this stage |

## support_97_99_member_specific_recovery.tsv
sha256: bc2d1d9101182683aa1d054d9db13b72e84d320d53d31f3b20b2fbbf19c45e35
rows: 360 (one row per acquisition, arm and CAMPI community member)
produced by: campi_candidate_clustering, H08 searches

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition |
| arm | id | `unclustered`, `cluster_097`, `cluster_099` |
| member | id | community member species |
| member_specific_REF_denominator | count of peptides | reference-search peptides assigned to this member only |
| recovered | count of peptides | of those, accepted in the arm |
| lost_vs_unclustered | count of peptides | member-specific peptides in unclustered, not in the arm |
| gained_vs_unclustered | count of peptides | member-specific peptides in the arm, not in unclustered |
