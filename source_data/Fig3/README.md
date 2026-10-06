# Source data: Figure 3

Four tables, each a byte-identical copy of the table the original producer read. Their
origins and the comparison with the Source Data package are in
`figures/Fig3/PROVENANCE.md`. The healthy-stool roster (58 acquisitions) is not affected by
the CapScan roster change for v28f, so no table here was replaced.

## transfer_pool_comparisons.tsv
sha256: a7618438c577a20b6037bf801b4daa89c181dcafb724066cfada229b06035947
rows: 420 (one row per target acquisition, experiment, database arm and peptide length floor:
58 targets x 3 transfer arms x 2 floors, plus 12 targets x 3 pool arms x 2 floors)
produced by: healthy_stool_database_transfer

| column | unit | meaning |
|---|---|---|
| sample | id | target acquisition (raw file name) |
| kind | id | `transfer` (database selected from one other acquisition) or `pool` (union of replicate databases) |
| group | id | sampling series: TechnicalReplicates, WorkflowReplicates, OnePiece (spatial sampling), Intra (intra-individual), Interindividual |
| arm | id | `own`; `T1_within_group` (donor from the same series); `T4_different_donor` (donor from the ten-participant inter-individual series); `pool_all12`; `pool_other11` |
| min_length | residues | lower peptide length bound: 5 (5-50 residues; v28d and Source Data) or 9 (9-50 residues; plotted in Fig. 3b-e from v30, the accepted 5-50 set filtered to >= 9 residues) |
| own | count of peptides | canonical target peptides at peptide q <= 0.01 with the target's own database |
| peptides | count of peptides | the same with this arm's database |
| overlap | count of peptides | peptides found with both databases |
| lost | count of peptides | own-database peptides not found with this arm's database |
| gained | count of peptides | peptides found only with this arm's database |
| count_retention_percent | % | peptides / own x 100 |
| own_set_recovery_percent | % | overlap / own x 100 |
| own_sequence_present | count of peptides | own-database peptides whose sequence occurs in this arm's FASTA |
| own_digest_compatible | count of peptides | of those, peptides that are also a compatible tryptic (K/R, not before P) digest product |
| own_digest_compatible_unidentified | count of peptides | digest-compatible own-database peptides not identified with this arm's database |
| database_accessions | count of proteins | protein entries in the searched FASTA |
| donor | id | acquisition whose predictions selected the database; the target itself for `own` transfer rows; `pooled` for pool rows |

## protein_group_pairwise.tsv
sha256: b45b6398514eebef4d92ac080f5df139f75fb2efcde79aebccfc475cab787317
rows: 1896 (one row per space, class and acquisition pair: 316 pairs x 3 classes x 2 spaces)
produced by: peptide_centric_downstream

Figure 3f uses only `space == groups` with `cls` host and microbial. The `space == KO` rows
are an earlier KO allocation that the KO refit superseded; they are kept because the file is
copied unchanged, and no figure or text number is computed from them.

| column | unit | meaning |
|---|---|---|
| space | id | `groups` (protein groups) or `KO` (superseded, unused) |
| design | id | sampling series, as `group` above |
| cls | id | `all`, `host` or `microbial` group-level class |
| i, j | id | the two acquisitions of the pair |
| r_per_pair | dimensionless | Pearson correlation of log2 intensities over features positive in both acquisitions |
| cv_per_pair_pct | % | median over those features of population SD / mean x 100 of the two raw intensities |
| n_features_per_pair | count of protein groups | features positive in both acquisitions |

## ko_pairwise.tsv
sha256: 2cf176c091f92f89162acfb091c3e3da6ca04886380f564c158040ac521a01f3
rows: 1264 (one row per profile and acquisition pair: 316 pairs x 4 profiles)
produced by: healthy_stool_ko_precision_refit

| column | unit | meaning |
|---|---|---|
| design | id | sampling series |
| i, j | id | the two acquisitions of the pair |
| r_per_pair | dimensionless | Pearson correlation of log2 KO intensities over KOs positive in both acquisitions |
| cv_per_pair_pct | % | median over those KOs of population SD / mean x 100 of the two raw intensities |
| n_features_per_pair | count of KOs | KOs positive in both acquisitions |
| profile | id | `new_all` (revised allocation, all origins; Fig. 3f upper triangles), `new_other_unresolved` (revised allocation, human-assigned signal excluded; Fig. 3f lower triangles), `old_all` and `old_microbial` (previous allocation; Methods line 747 comparison, not plotted) |

## design.tsv
sha256: 01ce64a542cfd664eba560aa99c03fc04546045951d74986b44527e55bc6c2c3
rows: 58 (one row per healthy-stool acquisition, in the order used for the heatmap axes)
produced by: peptide_centric_downstream (copied from the cohort design table)

| column | unit | meaning |
|---|---|---|
| sample_id | id | acquisition (raw file name) |
| cohort | id | `HSTOOL` |
| subject_id | id | participant label where the file name records it (inter-individual series only) |
| group | id | sampling series |
| group_ordinal | id | empty |
| aliquot | id | `stool` |
| design_class | id | descriptive class of the series |
| instrument | id | mass spectrometer label from the file name |
| operator | id | operator initials from the file name |
| acquisition_date | date | acquisition date from the file name |
| reacquisition_stamp | id | re-acquisition time stamp from the file name, where present |
| group_source, group_rule, subject_source, metadata_built | id | how the row was derived |
