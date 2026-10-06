# Source data: Supplementary Figure 6

All four tables were written by `figures/SFig06/prepare_source_data.py`, which copies the
original tables, removes the columns named below and changes nothing else (values are
copied as text), and counts participants and specimens before removing them
(`split_half_structure.tsv`). PROVENANCE.md of the unit lists the origin of each file and its sha256.

## split_half_acquisitions_726.tsv
sha256: 21a0dea9ce7045e3835ada95627e4d693d843fdd7fa308438391c7cff0cf1007
rows: 1452 (726 acquisitions x min_length 5 and 9; 15 CAMPI, 58 healthy stool, 306 PREDICT, 347 CapScan)
produced by: cohort_spectrum_split_selection, with the CapScan re-run on 347 acquisitions (half caches pooled over 347)
used by: panels b, c (rows with min_length 9) and the spectrum-split numbers of the v30 text

Columns `participant` and `specimen` were removed (raw PREDICT patient and specimen
identifiers). Not in the Source Data package.

| column | unit | meaning |
|---|---|---|
| index | id | task index of the spectrum-split run |
| cohort | id | CAMPI_MOCK, HSTOOL (healthy stool), PREDICT or CAPSCAN |
| sample | id | acquisition (raw-file stem) |
| min_length | residues | lower peptide-length bound of the view (5 or 9; upper bound 50) |
| full | count of peptides | accepted peptides, database selected from all predictions (FULL) |
| self | count of peptides | accepted peptides, database selected from odd-scan predictions (JOINT, same half as the scored spectra) |
| disjoint | count of peptides | accepted peptides, database selected from even-scan predictions (DISJOINT) |
| D_over_S | dimensionless | disjoint / self (panel b) |
| D_over_F | dimensionless | disjoint / full (panel c) |
| S_over_F | dimensionless | self / full |
| lost_from_self | count of peptides | peptides in self not in disjoint |
| gained_vs_self | count of peptides | peptides in disjoint not in self |
| shared_with_self | count of peptides | peptides in both self and disjoint |
| fraction_self_identities_retained | dimensionless | shared_with_self / self (identity retention, v28d line 101) |
| disjoint_fraction_in_full | dimensionless | fraction of disjoint peptides also in full |
| even_predictions | count of predictions | de novo predictions from even scans |
| odd_predictions | count of predictions | de novo predictions from odd scans |
| prediction_ratio | dimensionless | even_predictions / odd_predictions |
| even_selected_proteins | count of sequences | proteins selected from even-scan predictions |
| odd_selected_proteins | count of sequences | proteins selected from odd-scan predictions |
| selected_protein_ratio | dimensionless | even_selected_proteins / odd_selected_proteins |
| conversion_recovered | boolean | one of five PREDICT inputs reconverted from raw (Supplementary Methods 7) |

## same_spectrum_benchmark_0_predict.tsv
sha256: bf59b19aa057f36621559267b74baf70d8923ac8739dce996ae722c40da56455
rows: 1796 (one row per acquisition, comparison and length view; PREDICT)
produced by: clinical_public_only_comparison (UPSTREAM.json)
used by: panels d, e; v28d line 123 (rows with min_length 5, test_arm FL_public, reference_arm MG: 296)

## same_spectrum_benchmark_1_capscan.tsv
sha256: 37b1cd553132ab1338efdf3dd8b93f5514078fe8d057dd4064e3666c8cc5355c
rows: 1790 (as above; CapScan)
produced by: clinical_public_only_comparison (CapScan searches from capscan_public_only_search)
used by: panels d, e; v28d line 123 (rows with min_length 5, FL_public vs MG: 282)

Removed columns: `specimen_id` in both tables (pseudonymised tube identifier, not used),
and `group` in the PREDICT table only. In the PREDICT sheet `group` is the study visit
code (S = scheduled visit, A = admission, R = readmission), a clinical event that may not
sit next to a participant id (CONVENTIONS.md section 8); no panel or number uses it.
Otherwise equal to Source Data package sheets `Figure4/benchmark_0_per_acquisition.tsv`
and `Figure4/benchmark_1_per_acquisition.tsv`.

Participants: 94 PREDICT and 15 CapScan `subject_id` values in the rows used by d,e
(9-50 residues, FL_public vs MG).

| column | unit | meaning |
|---|---|---|
| cohort | id | PREDICT or CAPSCAN |
| sample | id | acquisition (raw-file stem) |
| subject_id | id | pseudonymised participant id (as in the Source Data package) |
| group | id | CapScan table only: sampling region (duodenum, jejunum, ileum, colon, stool) |
| aliquot | id | sample fraction: `stool` for every PREDICT row; for CapScan `pellet` (intestinal device samples) or `main` (stool) |
| test_arm | id | FL_public (public-sequence FastaLake) or FL_cohortMG (FastaLake on the specimen metagenome) |
| reference_arm | id | MG (matched metagenome search) or FL_cohortMG |
| min_length | residues | lower peptide-length bound (5 or 9; upper bound 50) |
| n_test | count of peptides | accepted peptides of the test arm |
| n_reference | count of peptides | accepted peptides of the reference arm |
| overlap | count of peptides | peptides accepted by both arms |
| gained | count of peptides | in test only |
| lost | count of peptides | in reference only |
| count_ratio | dimensionless | n_test / n_reference |
| reference_recovery | dimensionless | overlap / n_reference (panel e, x 100) |
| test_recovery | dimensionless | overlap / n_test |

## split_half_structure.tsv
sha256: 1107c62d83d18c8f5c049909a4d77264fbe10b81b247fbe927a73328df191f5f
rows: 4 (one per cohort; 726 roster, 9-50 residues)
produced by: `prepare_source_data.py`, function `structure`, from the `participant` and
`specimen` columns of the 726 table before their removal
used by: STATS.tsv (`distinct_or_repeated`), notebook section "Panel b"

| column | unit | meaning |
|---|---|---|
| cohort | id | CAMPI_MOCK, HSTOOL, PREDICT or CAPSCAN |
| n_acquisitions | count of acquisitions | acquisitions in panels b,c |
| n_participants | count of participants | distinct non-empty participant identifiers |
| n_specimens | count of specimens | distinct non-empty specimen identifiers |
| n_acquisitions_without_participant | count of acquisitions | acquisitions with an empty participant identifier |
| n_acquisitions_without_specimen | count of acquisitions | acquisitions with an empty specimen identifier |
