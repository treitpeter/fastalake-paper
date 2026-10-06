# Source data: Supplementary Figure 3

CAMPI mock-community tables for Supplementary Fig. 3 and the CAMPI statements on v28d lines 83, 89, 253, 417 and 685. Peptides are I/L-canonical (I written as L), unmodified sequences accepted at peptide q <= 0.01. All files are complete; nothing was subsampled.

## CAMPI_member_support_recount.tsv
sha256: a230cf25e366bf279891b80875be88ddc8a7838673e4faac74fa8c61a9d6d107
rows: 80 (one row per arm, minimum length and community member (5 x 2 x 8))
produced by: A14_campi_member_gap (recounted by `panels.py:stat_recount_membership`)

| column | unit | meaning |
|---|---|---|
| arm | id | `FL` FastaLake database; `FL_R` FL with the eight reference proteomes appended after selection; `FL_R_M` FL with reference proteomes and MG appended; `REF` the eight reference proteomes; `MG` the metagenome-derived protein database supplied with the CAMPI SIHUMIx data (DNA-only input; see figures/SFig03/OPEN_ITEMS.md) |
| min_length | residues | lower peptide length bound, 5 or 9 (upper bound 50 for both) |
| member | id | community member as named in the reference proteome set |
| pooled_specific_peptides | count of peptides | distinct peptides assigned to this member only, pooled over the 15 acquisitions |
| acquisitions_with_support | count of acquisitions | acquisitions with at least one such peptide |
| total_acquisitions | count of acquisitions | acquisitions searched (15) |

## CAMPI_member_acquisition_counts.tsv
sha256: 4e995436231bf70ab90cfd44e68b56b3b3d0757f015a21967470b96eca138486
rows: 1200 (one row per acquisition, arm, minimum length and member (15 x 5 x 2 x 8), zeros included)
produced by: A14_campi_member_gap (recounted by `panels.py:stat_recount_membership`)

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition (S11 was acquired in four fractions, each a separate acquisition) |
| arm | id | `FL` FastaLake database; `FL_R` FL with the eight reference proteomes appended after selection; `FL_R_M` FL with reference proteomes and MG appended; `REF` the eight reference proteomes; `MG` the metagenome-derived protein database supplied with the CAMPI SIHUMIx data (DNA-only input; see figures/SFig03/OPEN_ITEMS.md) |
| min_length | residues | 5 or 9; upper bound 50 |
| member | id | community member |
| specific_peptides | count of peptides | distinct accepted peptides assigned to this member only |

## CAMPI_reference_recovery.tsv
sha256: 7f2d7a7586c57b47332ba62d7348aceac33aea79dbed3218c46d2dc9d4cb8551
rows: 150 (one row per acquisition, Sage arm and minimum length (15 x 5 x 2))
produced by: A14_campi_member_gap (recounted by `panels.py:stat_recount_membership`)

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition |
| arm | id | `FL` FastaLake database; `FL_R` FL with the eight reference proteomes appended after selection; `FL_R_M` FL with reference proteomes and MG appended; `REF` the eight reference proteomes; `MG` the metagenome-derived protein database supplied with the CAMPI SIHUMIx data (DNA-only input; see figures/SFig03/OPEN_ITEMS.md) |
| min_length | residues | 5 or 9; upper bound 50 |
| peptides | count of peptides | accepted distinct I/L-canonical peptides in this arm |
| reference | count of peptides | accepted peptides in the REF arm of the same acquisition |
| overlap | count of peptides | peptides accepted in both |
| gained | count of peptides | in this arm, not in REF |
| lost | count of peptides | in REF, not in this arm |
| count_ratio | dimensionless | peptides / reference |
| reference_recovery | fraction | overlap / reference (panel c, times 100) |

## CAMPI_within_MSFragger.tsv
sha256: 13490970bb7873ff2dd349f8de1256e656329bc2c8433fc448702c8c98ae5c47
rows: 60 (one row per acquisition, comparator and minimum length (15 x 2 x 2))
produced by: campi_second_engine

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition |
| comparator | id | `REF` or `MG` |
| min_length | residues | 5 or 9; upper bound 50 |
| FL | count of peptides | peptides accepted by MSFragger + Percolator with the FL database |
| reference | count of peptides | peptides accepted with the comparator database (REF or MG; the column name is the producer's) |
| intersection | count of peptides | accepted in both |
| FL_only | count of peptides | FL only |
| reference_only | count of peptides | comparator only |
| reference_recovery | fraction | intersection / reference (panel d, times 100) |
| count_ratio | dimensionless | FL / reference |

## second_engine_between_engine_identities.tsv
sha256: d455994adac00eb735a6955bd13ab265de2ca2bdbc2c04d5f42a4dd6c6274501
rows: 90 (one row per acquisition, database and minimum length (15 x 3 x 2); Source Data only, not drawn)
produced by: campi_second_engine

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition |
| arm | id | `FL_UNIPROT` (the FL database), `REF`, `MG` |
| min_length | residues | 5 or 9; upper bound 50 |
| MSFragger | count of peptides | accepted by MSFragger + Percolator |
| Sage | count of peptides | accepted by Sage with the same database |
| intersection | count of peptides | accepted by both engines |
| MSFragger_only | count of peptides | MSFragger only |
| Sage_only | count of peptides | Sage only |
| Sage_identity_recovered | fraction | intersection / Sage |

## CAMPI_MSFragger_search_receipts.tsv
sha256: 03b0c93ada38888117f940780613b475be98a98e36068cb03e4306ec3c6aead4
rows: 45 (one row per MSFragger + Percolator search (15 acquisitions x 3 databases))
produced by: campi_second_engine (COMPLETE.json of each search, exported field by field without paths)

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI acquisition |
| arm | id | `FL_UNIPROT` (the FL database), `REF`, `MG` |
| status | id | completion status written by the search (`PASS`) |
| engine_and_confidence | text | engine and confidence procedure as recorded by the search |
| accepted_5_50 | count of peptides | accepted I/L-canonical peptides, 5-50 residues, peptide q <= 0.01 |
| accepted_9_50 | count of peptides | the nested 9-50-residue subset |
| target_entries | count of sequences | target protein entries in the searched FASTA |
| raw_sha256 | id | sha256 of the input mzML |
| fasta_sha256 | id | sha256 of the searched target FASTA |
| slurm_job_id | id | Slurm job that ran the search |

## CAMPI_membership_<sample>.tsv (15 files)
One file per acquisition: every peptide accepted by Sage in at least one of the five arms, with its membership in the eight reference proteomes.

| file | sha256 | rows |
|---|---|---|
| CAMPI_membership_S01.tsv | fe3f27ac42f772792572e8400c30e2665eb7ca8b3e9d075f0b5f6e0735035fb2 | 28989 |
| CAMPI_membership_S02.tsv | 4a64e2930e4d27bd74b913b6863e742a8663195a5364f01c86422ce39b55c78c | 25894 |
| CAMPI_membership_S03.tsv | c065d0f2714bfd2e7999666185c1f299ca49ad65a7acb7beaa576fec9f865205 | 33864 |
| CAMPI_membership_S04.tsv | 2831b3686f6ad449fba4e94657827f719993b2fb91e7eb0e76e5d5770fe1ef17 | 32890 |
| CAMPI_membership_S05.tsv | cce302fcbfbd71bc4bb9677863bfc5aa68ad5b5d730cfeab6e822e8ee4f1beb5 | 38600 |
| CAMPI_membership_S06.tsv | 4b81dbb3d5c4ea81e250e274971fd3501f41811d366b6fde8e96845ef18b32ee | 39187 |
| CAMPI_membership_S07.tsv | ec883caeb75013b3ca778ee8d4255b633fde058d99b5cf11e5bc2e38355979a7 | 20010 |
| CAMPI_membership_S08.tsv | 5f1f92c20b3746aec6927465b56f17755545abf88ac8de64afce1a4be9ddd26f | 26239 |
| CAMPI_membership_S09.tsv | dc2aef1c5d7ac110fc202d6b65d8287f0aaa6d28f0542c289d3491ca356feea5 | 24011 |
| CAMPI_membership_S10.tsv | 7520e41bca7aecc7468f7606359c1f80816c51f95396a67e44d377c95031aee0 | 32763 |
| CAMPI_membership_S11_Fraction1.tsv | 55cbde9f72cff8bb5aa72836abcee99047ba581792203d32c5aeb79a7052ee65 | 20251 |
| CAMPI_membership_S11_Fraction2.tsv | 54be84a4f676d0d4d81a902c0c739c20c91338eeeec997405fb13a1c22af7d1b | 20609 |
| CAMPI_membership_S11_Fraction3.tsv | 0620cde139e1acf8eb4cd47c997619dbdbccdeba6917a8d507302bb073df3b32 | 21701 |
| CAMPI_membership_S11_Fraction4.tsv | e93e188f04caf0aa359bb5158b2ced259286c12bddd59caa872b32fa9e807d05 | 12221 |
| CAMPI_membership_S12.tsv | d68c1485743dfc9283263e3c32f8392bd4712a556181c37acede4b00cd2f16a6 | 22410 |

rows: one row per distinct accepted peptide of that acquisition
produced by: A14_campi_member_gap (`results/<sample>/peptide_membership.tsv`, byte-identical)

| column | unit | meaning |
|---|---|---|
| peptide | sequence | accepted peptide, I/L-canonical, 5-50 residues |
| n_members | count of members | number of the eight reference proteomes that contain the peptide (0 = none) |
| members | id | those members, separated by `;` (empty when n_members is 0) |
| arms | id | Sage arms that accepted the peptide in this acquisition, separated by `;` (FL, FL_R, FL_R_M, REF, MG) |
