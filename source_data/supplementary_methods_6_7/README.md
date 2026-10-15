# Source data: Supplementary Methods 6-7

Tables read by `text_statistics/supplementary_methods_6_7/`. Acquisition names are raw-file stems. Raw PREDICT patient identifiers, specimen and tube identifiers are not present; PREDICT participants appear only as the Source Data package pseudonyms (`PREDICT_P###`), in the split-half and cross-fit tables alike. PREDICT visit codes are not present: the `group` column is blank for PREDICT rows. `cluster_order` in crossfit_per_acquisition.tsv is the rank of the original participant identifier within its cohort; it is needed to reproduce the original bootstrap draws and discloses only that rank. No clinical column is present.

Origins, origin sha256 and transformations are in `text_statistics/supplementary_methods_6_7/PROVENANCE.md`. "compared, agrees" means the notebook recomputes the values and asserts agreement within 1e-9.

## acquisition_sizes_726.tsv

sha256: d90de4bfc11efbb9584864a5220d6dfd84b96d291a095962982f611c8549e4cf  
rows: 726  
produced by: cohort_database_sizes  
used for v28d lines: 665, 673-675, 713  
Source Data package: Review_v16/Database_sizes/database_sizes_710.tsv holds 710 rows; this 726-row table (with 16 CapScan participant-15 acquisitions) is not yet in the package

Entries of every selected FASTA (726 acquisitions) and the reservoir each was selected from.

| column | meaning |
|---|---|
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| reservoir_entries | entries of the reservoir (count of sequences) |
| below_human_reference | entries < 20,431 |
| bytes | FASTA size in bytes |
| sha256 | FASTA sha256 |
| entries | FASTA entries |
| residues | residues in the FASTA |
| unique_exact_sequences | distinct exact sequences |
| duplicate_accession_count | repeated accessions |
| sequence_set_sha256 | hash of the sorted sequence set |

## campi_member_summary.tsv

sha256: bdf6d5f0ded8bc8455f95bd1ac8c4ab52af38294ac7030959a269715873591b8  
rows: 80  
produced by: campi_reference_gap  
used for v28d lines: 683  
Source Data package: Figure2/member_summary.tsv (byte copy)

Member-specific peptides per arm and community member.

| column | meaning |
|---|---|
| arm | search or control arm |
| min_length | lower peptide length bound (residues) |
| member | community member |
| pooled_specific_peptides | member-specific peptides, union over 15 acquisitions |
| acquisitions_ge1 | acquisitions with at least one |
| acquisitions_ge10 | acquisitions with at least ten |
| median_specific_peptides | median per acquisition |

## campi_recovery_summary.tsv

sha256: 74658e932aea9f40a62245b6f3153f578289343182a02581ecca83736c060b13  
rows: 10  
produced by: campi_reference_gap  
used for v28d lines: 687  
Source Data package: Figure2/recovery_summary.tsv (byte copy; pooled_reference compared, agrees)

Pooled recovery of reference peptides per arm.

| column | meaning |
|---|---|
| arm | search or control arm |
| min_length | lower peptide length bound (residues) |
| acquisitions | acquisitions pooled |
| pooled_peptides | union of accepted peptides |
| pooled_reference | union of reference peptides |
| pooled_overlap | intersection |
| pooled_reference_recovery | overlap / reference |
| median_acquisition_count_ratio | median per-acquisition count ratio |
| median_acquisition_reference_recovery | median per-acquisition recovery |
| known_member_specific | peptides specific to one member |
| known_member_shared | peptides shared by members |
| absent_eight_member_reference | peptides absent from all eight member proteomes |

## campi_reference_fasta_records.tsv

sha256: f7b9e2df6ce6822cc5a64ecc34b53055ecba723b13ded06e9ab84cada693f1c8  
rows: 29,673  
produced by: none (supplied CAMPI reference FASTA)  
used for v28d lines: 683  
Source Data package: not in the package; reference rows of source_data/supplementary_methods_1_5/campi_target_fasta_records.tsv (unit supplementary_methods_1_5)

One row per record of the supplied CAMPI reference FASTA (29,673).

| column | meaning |
|---|---|
| database | reference |
| record | 0-based position |
| header_token | first header token |
| header_kind | uniprot, assembly_contig or other (contaminant) |
| header_has_trailing_whitespace | header ends in whitespace |
| length | residues |
| sequence_sha256 | sha256 of the sequence |

## campi_reference_peptides.tsv

sha256: 3af55b567347ecce79a1f48f4fdbec4c5970a3d7a4914439689483256ccf4283  
rows: 331,674  
produced by: campi_variability  
used for v28d lines: 687  
Source Data package: concatenation of CAMPI_variability/results/<acquisition>/REF_5_frequencies.tsv (values unchanged)

Accepted reference-search peptides (5-50 residues) per CAMPI acquisition, with accepted spectra per peptide.

| column | meaning |
|---|---|
| sample | acquisition (raw-file stem) |
| peptide | accepted canonical peptide |
| accepted_psms | accepted spectra of the peptide |

## campi_spectral_inventory.tsv

sha256: 41c3c0f554b319bff64e66b141bae2fdfd5ff0a4ed273e9d571cfe67c7bbe65f  
rows: 15  
produced by: campi_variability  
used for v28d lines: 683-687, 699  
Source Data package: CAMPI_variability/summary/spectral_inventory.tsv (byte copy)

MS2 spectra per CAMPI acquisition.

| column | meaning |
|---|---|
| sample | acquisition (raw-file stem) |
| ms2 | MS2 spectra |
| rt_min | first retention time (min) |
| rt_max | last retention time (min) |
| mzml_sha256 | mzML sha256 |

## clustering_per_acquisition_identifications.tsv

sha256: 0187448fe338050276ab9b8748e15ef820e1888ba62019556e08487788511c06  
rows: 105  
produced by: campi_candidate_clustering  
used for v28d lines: 695  
Source Data package: not yet in the package (Supplementary_Tables/Supplementary_Table_1.xlsx holds the 97/99% arms only)

One row per CAMPI acquisition and clustering arm (15 x 7).

| column | meaning |
|---|---|
| sample | acquisition (raw-file stem) |
| arm | search or control arm |
| identity | Linclust --min-seq-id (1.0 for unclustered) |
| protein_entries | entries in the searched FASTA |
| protein_entry_ratio | entries / unclustered entries |
| peptides | accepted peptides |
| baseline_peptides | accepted peptides, unclustered |
| peptide_count_ratio | peptides / baseline_peptides |
| peptide_reduction_pct | 100 x (1 - ratio) |
| baseline_peptide_identities_retained | unclustered peptides also accepted |
| baseline_identity_recovery | retained / baseline |
| REF_peptides | reference-search peptides |
| REF_recovered | reference peptides accepted in this arm |
| REF_recovery | REF_recovered / REF_peptides |
| lost | accepted only unclustered |
| gained | accepted only in this arm |
| lost_sequence_absent | lost and absent from this arm's FASTA |
| lost_sequence_available | lost but present in this arm's FASTA |

## common_fasta_union_summary.tsv

sha256: 144898005afcaeb575f7dc0d0f5430eb92921e284eefd8ab9188bf8dc0e28889  
rows: 1  
produced by: healthy_stool_common_fasta  
used for v28d lines: 751  
Source Data package: Figure3/common_fasta_union_summary.json (flattened)

Exact-sequence union of the 58 selected healthy-stool FASTAs.

| column | meaning |
|---|---|
| source_fastas | FASTAs merged |
| input_records | records before deduplication |
| duplicate_records_removed | duplicates |
| unique_sequences | union size |
| total_residues | residues |
| union_sha256 | union FASTA sha256 |
| identity | deduplication rule |
| status | producer status |

## confidence_fusion_all.tsv

sha256: 19feda566021def7e2321d8a838b636db07f7c59060de88a32f351600720206b  
rows: 1,080  
produced by: campi_fusion_375m  
used for v28d lines: 633  
Source Data package: Review_v16/Confidence/fusion_all.tsv (byte copy)

CAMPI fusion searches under added confidence gates.

| column | meaning |
|---|---|
| index | task index of the producing run |
| sample | acquisition (raw-file stem) |
| seed | random seed of the control |
| gate | acceptance gate: q_only, or q plus PSM PEP <= 0.01 (pep01) or <= 0.001 (pep001) |
| q | peptide q threshold |
| min_length | lower peptide length bound (residues) |
| max_length | upper peptide length bound (residues) |
| n_canonical | as in the package or producer table named above |
| paired_numerator | as in the package or producer table named above |
| n_unpaired_accepted | as in the package or producer table named above |
| n_target_only | as in the package or producer table named above |
| n_entrapment_only | as in the package or producer table named above |
| n_shared | as in the package or producer table named above |
| diagnostic | (paired_numerator + U) / N |

## confidence_ordinary_all_pairs.tsv

sha256: bfb591186d9aace576461f7e884a149a1058b2452f23bc7c0f51ac3b217beabc  
rows: 21,276  
produced by: clinical_public_only_comparison  
used for v28d lines: 633, 777  
Source Data package: Review_v16/Confidence/ordinary_all_pairs.tsv minus subject_id, specimen_id, specimen_basis, aliquot; PREDICT group values blanked

Matched FastaLake (public-sequence) and metagenome counts per acquisition, gate, q and length.

| column | meaning |
|---|---|
| index | task index of the producing run |
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| gate | acceptance gate: q_only, or q plus PSM PEP <= 0.01 (pep01) or <= 0.001 (pep001) |
| q | peptide q threshold |
| min_length | lower peptide length bound (residues) |
| max_length | upper peptide length bound (residues) |
| n_fl | FastaLake peptides |
| n_mg | metagenome peptides |
| shared | both |
| fl_only | FastaLake only |
| mg_only | metagenome only |
| removed_fl | removed by the added gate |
| removed_mg | removed |
| group | sampling group: healthy-stool design or CapScan region or stool; blank for PREDICT (visit codes removed) |
| gain_pct | 100 (n_fl/n_mg - 1) |
| fl_retained | fraction kept by the gate |
| mg_retained | fraction kept |

## confidence_ordinary_summary.tsv

sha256: 5108ef1a8bc0742045240a3e0381f29656939fd0ba259b109f26ebe9d86767b9  
rows: 72  
produced by: clinical_public_only_comparison  
used for v28d lines: 777  
Source Data package: Review_v16/Confidence/ordinary_summary.tsv (byte copy; compared, agrees)

Producer's medians.

| column | meaning |
|---|---|
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| gate | acceptance gate: q_only, or q plus PSM PEP <= 0.01 (pep01) or <= 0.001 (pep001) |
| q | peptide q threshold |
| min_length | lower peptide length bound (residues) |
| max_length | upper peptide length bound (residues) |
| n | as in the package or producer table named above |
| defined_ratios | as in the package or producer table named above |
| median_gain_pct | as in the package or producer table named above |
| median_fl_retained | as in the package or producer table named above |
| median_mg_retained | as in the package or producer table named above |
| median_fl_only | as in the package or producer table named above |
| median_mg_only | as in the package or producer table named above |
| zero_mg | as in the package or producer table named above |

## confidence_pep01_zero_yield_cases.tsv

sha256: 9abb6f353099da3728201a2399ab57d197c2ec1880b9d50e896479bd5cb591dd  
rows: 3  
produced by: clinical_public_only_comparison  
used for v28d lines: 633  
Source Data package: Review_v16/Confidence/pep01_zero_yield_cases.tsv minus specimen columns; PREDICT group values blanked

PREDICT FastaLake acquisitions with zero yield under the PEP <= 0.01 gate.

| column | meaning |
|---|---|
| index | task index of the producing run |
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| gate | acceptance gate: q_only, or q plus PSM PEP <= 0.01 (pep01) or <= 0.001 (pep001) |
| q | peptide q threshold |
| min_length | lower peptide length bound (residues) |
| max_length | upper peptide length bound (residues) |
| n_fl | as in the package or producer table named above |
| n_mg | as in the package or producer table named above |
| shared | as in the package or producer table named above |
| fl_only | as in the package or producer table named above |
| mg_only | as in the package or producer table named above |
| removed_fl | as in the package or producer table named above |
| removed_mg | as in the package or producer table named above |
| group | sampling group: healthy-stool design or CapScan region or stool; blank for PREDICT (visit codes removed) |
| gain_pct | as in the package or producer table named above |
| fl_retained | as in the package or producer table named above |
| mg_retained | as in the package or producer table named above |

## crossfit_fold_caches.tsv

sha256: d85dc40c72d144c2d2bd0dddf739a618118434483ef0f41a56929372f588d800  
rows: 6,370  
produced by: capscan_crossfit_10fold  
used for v28d lines: 773  
Source Data package: not in the package

Out-of-fold prediction tables written for each cohort and fold cache (2 cohorts x 10 folds; 306 PREDICT and 331 CapScan acquisitions per fold).

| column | meaning |
|---|---|
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| fold | fold k (0-9) whose spectra were withheld |
| K | number of folds |
| fold_seed | prefix of the fold key (sha256 of '<fold_seed>:<native scan>', first eight bytes, modulo K) |
| n_acquisitions | acquisitions recorded by the cache |
| sample | acquisition (raw-file stem) |
| oof_rows | rows of the acquisition's out-of-fold prediction table |
| oof_csv_sha256 | sha256 of that table |

## crossfit_independent_lda_check.tsv

sha256: 2748f4e5eab34f5cbe624d962b9d2ab77ffbc24cfd3ba748dfab1c7bb846eae7  
rows: 24  
produced by: capscan_crossfit_10fold  
used for v28d lines: 777  
Source Data package: not in the package

Independent verifier of 8 cross-fit acquisitions: producer pooled-rescoring counts against an independent closed-form Fisher discriminant with the same picked-peptide q, per arm.

| column | meaning |
|---|---|
| i | roster index |
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| structural_checks_ok | fold, leakage and out-of-fold hash checks all passed |
| arm | search or control arm |
| producer_rescoring | lda (fitted discriminant) or fallback (Sage heuristic score) |
| producer_accepted_5 | producer accepted peptides, 5-50 |
| producer_accepted_9 | 9-50 |
| producer_score_independent_q_equal | producer scores with the independent q give the producer's accepted set |
| closed_form_accepted_5 | closed-form discriminant accepted peptides, 5-50 |
| closed_form_accepted_9 | 9-50 |
| closed_form_jaccard | Jaccard index of closed-form and producer accepted sets |

## crossfit_k1_control.tsv

sha256: a13580e4fb8ed76ddd9adc4fd1cb020dd8f14c74fa68cb1330c5af7479ee6e7b  
rows: 24  
produced by: capscan_crossfit_10fold  
used for v28d lines: 775, 777  
Source Data package: not in the package

One-fold (no cross-fitting) control on 6 pilots.

| column | meaning |
|---|---|
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| arm | search or control arm |
| min_length | lower peptide length bound (residues) |
| selection_equals_archived_fasta | selected FASTA equals the archived one |
| sage_equals_frozen | Sage accepted set equals the original |
| sage_k1 | accepted peptides |
| frozen | original accepted peptides |
| pooled_rescoring_equals_frozen | equal through the pooled rescoring |
| pooled_rescoring_k1 | accepted through pooled rescoring |

## crossfit_lda_fallback_acquisitions.tsv

sha256: b4fb01d1633b1ae81781dbed1210549185231339e96bdf69a8ab59103ec4328a  
rows: 3  
produced by: capscan_crossfit_10fold  
used for v28d lines: 777  
Source Data package: Supplementary_Figure_06/metagenome_crossfit/lda_fallback_acquisitions.tsv (PREDICT group values blanked; otherwise unchanged)

The three acquisitions with a singular pooled fit, with producer and independent-refit counts.

| column | meaning |
|---|---|
| i | roster index |
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| group | sampling group: healthy-stool design or CapScan region or stool; blank for PREDICT (visit codes removed) |
| subject_id | as in the package or producer table named above |
| arms_with_fallback | as in the package or producer table named above |
| FL_pooled_5 | as in the package or producer table named above |
| MG_pooled_5 | as in the package or producer table named above |
| gain_cf_pooled_5 | as in the package or producer table named above |
| FL_independent_refit_5 | as in the package or producer table named above |
| MG_independent_refit_5 | as in the package or producer table named above |
| gain_independent_refit_5 | as in the package or producer table named above |
| FL_pooled_9 | as in the package or producer table named above |
| MG_pooled_9 | as in the package or producer table named above |
| gain_cf_pooled_9 | as in the package or producer table named above |
| FL_independent_refit_9 | as in the package or producer table named above |
| MG_independent_refit_9 | as in the package or producer table named above |
| gain_independent_refit_9 | as in the package or producer table named above |

## crossfit_per_acquisition.tsv

sha256: 08f0913375b0b0cc0de3eb64407abc3bc042cea0325c55e75e4993a6545c2adc  
rows: 591  
produced by: capscan_crossfit_10fold  
used for v28d lines: 775-777  
Source Data package: Supplementary_Figure_06/metagenome_crossfit/per_acquisition.tsv plus column cluster_order; PREDICT group values blanked

Cross-fit counts and gains per matched acquisition (578). Column meanings are in the package README (Supplementary_Figure_06/metagenome_crossfit/README.md). cluster_order is the rank of the original participant identifier within its cohort, which fixes the cluster order of the original bootstrap; it discloses that rank only, not the identifiers. PREDICT visit codes in group are blanked.

| column | meaning |
|---|---|
| i | roster index |
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| group | sampling group: healthy-stool design or CapScan region or stool; blank for PREDICT (visit codes removed) |
| subject_id | PREDICT_P### pseudonym, CapScan participant number, or NOID_<acquisition> when no participant was recorded |
| fold_ms2_total | as in the package or producer table named above |
| fl_fasta_sequences_median | as in the package or producer table named above |
| lda_fallback | as in the package or producer table named above |
| FL_pooled_5 | as in the package or producer table named above |
| FL_union_5 | as in the package or producer table named above |
| FL_pert1_5 | as in the package or producer table named above |
| FL_pert2_5 | as in the package or producer table named above |
| FL_pooled_naive_5 | as in the package or producer table named above |
| FULL_pooled_5 | as in the package or producer table named above |
| FULL_union_5 | as in the package or producer table named above |
| FULL_pert1_5 | as in the package or producer table named above |
| FULL_pert2_5 | as in the package or producer table named above |
| FULL_pooled_naive_5 | as in the package or producer table named above |
| MG_pooled_5 | as in the package or producer table named above |
| MG_union_5 | as in the package or producer table named above |
| MG_pert1_5 | as in the package or producer table named above |
| MG_pert2_5 | as in the package or producer table named above |
| MG_pooled_naive_5 | as in the package or producer table named above |
| frozen_fl_5 | as in the package or producer table named above |
| frozen_mg_5 | as in the package or producer table named above |
| gain_cf_pooled_5 | as in the package or producer table named above |
| gain_cf_union_5 | as in the package or producer table named above |
| gain_cf_pert1_5 | as in the package or producer table named above |
| gain_cf_pert2_5 | as in the package or producer table named above |
| gain_full_pooled_5 | as in the package or producer table named above |
| gain_full_union_5 | as in the package or producer table named above |
| gain_frozen_5 | as in the package or producer table named above |
| cf_vs_full_pooled_5 | as in the package or producer table named above |
| cf_vs_full_union_5 | as in the package or producer table named above |
| calib_full_minus_frozen_5 | as in the package or producer table named above |
| mg_pooled_vs_frozen_5 | as in the package or producer table named above |
| gain_odd_full_5 | as in the package or producer table named above |
| gain_odd_self_5 | as in the package or producer table named above |
| gain_odd_disjoint_5 | as in the package or producer table named above |
| FL_pooled_9 | as in the package or producer table named above |
| FL_union_9 | as in the package or producer table named above |
| FL_pert1_9 | as in the package or producer table named above |
| FL_pert2_9 | as in the package or producer table named above |
| FL_pooled_naive_9 | as in the package or producer table named above |
| FULL_pooled_9 | as in the package or producer table named above |
| FULL_union_9 | as in the package or producer table named above |
| FULL_pert1_9 | as in the package or producer table named above |
| FULL_pert2_9 | as in the package or producer table named above |
| FULL_pooled_naive_9 | as in the package or producer table named above |
| MG_pooled_9 | as in the package or producer table named above |
| MG_union_9 | as in the package or producer table named above |
| MG_pert1_9 | as in the package or producer table named above |
| MG_pert2_9 | as in the package or producer table named above |
| MG_pooled_naive_9 | as in the package or producer table named above |
| frozen_fl_9 | as in the package or producer table named above |
| frozen_mg_9 | as in the package or producer table named above |
| gain_cf_pooled_9 | as in the package or producer table named above |
| gain_cf_union_9 | as in the package or producer table named above |
| gain_cf_pert1_9 | as in the package or producer table named above |
| gain_cf_pert2_9 | as in the package or producer table named above |
| gain_full_pooled_9 | as in the package or producer table named above |
| gain_full_union_9 | as in the package or producer table named above |
| gain_frozen_9 | as in the package or producer table named above |
| cf_vs_full_pooled_9 | as in the package or producer table named above |
| cf_vs_full_union_9 | as in the package or producer table named above |
| calib_full_minus_frozen_9 | as in the package or producer table named above |
| mg_pooled_vs_frozen_9 | as in the package or producer table named above |
| gain_odd_full_9 | as in the package or producer table named above |
| gain_odd_self_9 | as in the package or producer table named above |
| gain_odd_disjoint_9 | as in the package or producer table named above |
| cluster_order | rank of the original participant identifier within the cohort (ordinal only) |

## crossfit_pilot_folds.tsv

sha256: 5193552ebffa42dc827f33284275cb361e6a07b43bdda54a4416d8eac3816581  
rows: 60  
produced by: capscan_crossfit_10fold  
used for v28d lines: 773  
Source Data package: not in the package

MS2 spectra per fold in the six pilots.

| column | meaning |
|---|---|
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| fold | fold 0-9 |
| fold_ms2 | MS2 spectra in the fold |
| fl_fasta_sequences | sequences in that fold's selected FASTA |

## crossfit_q_validation.tsv

sha256: f6ca018041098d3dc17e5786db472b37e079a9a542fcc6931bf50abf88c1b684  
rows: 36  
produced by: capscan_crossfit_10fold  
used for v28d lines: 775, 777  
Source Data package: not in the package

q re-implementation and full rescoring against Sage on 36 unsplit searches (6 pilots x 6 searches).

| column | meaning |
|---|---|
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| search | search checked |
| psms | PSMs |
| accepted_sage | accepted by Sage |
| accepted_reimpl | accepted by the re-implemented q |
| q_reimplementation_equal | identical accepted set |
| q_threshold_disagreements | PSMs on different sides of 0.01 |
| full_rescoring_equal | identical accepted set after refitting the discriminant |
| full_rescoring_max_abs_score_diff | largest score difference |

## crossfit_summary.tsv

sha256: b9a2be7d3c2019d81527e41bbbdc729564f239d5deb187b0592de67f767bf3f6
rows: 112  
produced by: capscan_crossfit_10fold  
used for v28d lines: 777  
Source Data package: Supplementary_Figure_06/metagenome_crossfit/summary.tsv (byte copy; compared, agrees)

Producer's stratum medians and intervals.

| column | meaning |
|---|---|
| stratum | as in the package or producer table named above |
| min_length | lower peptide length bound (residues) |
| comparison | as in the package or producer table named above |
| n_acquisitions | as in the package or producer table named above |
| n_participants | as in the package or producer table named above |
| median | as in the package or producer table named above |
| ci95_low | as in the package or producer table named above |
| ci95_high | as in the package or producer table named above |
| median_exact | as in the package or producer table named above |
| ci95_low_exact | as in the package or producer table named above |
| ci95_high_exact | as in the package or producer table named above |
| n_pos | as in the package or producer table named above |
| n_neg | as in the package or producer table named above |
| n_zero | as in the package or producer table named above |
| n_nan | as in the package or producer table named above |

## crossfit_summary_by_region.tsv

sha256: 8e39f6a3bb1f5d3ad345f485bb7646051eb713c7f5d44325a5b40c901fe244cf  
rows: 140  
produced by: capscan_crossfit_10fold  
used for v28d lines: 777  
Source Data package: Supplementary_Figure_06/metagenome_crossfit/summary_capscan_by_region.tsv (byte copy; compared, agrees)

Producer's region medians and intervals.

| column | meaning |
|---|---|
| stratum | as in the package or producer table named above |
| min_length | lower peptide length bound (residues) |
| comparison | as in the package or producer table named above |
| n_acquisitions | as in the package or producer table named above |
| n_participants | as in the package or producer table named above |
| median | as in the package or producer table named above |
| ci95_low | as in the package or producer table named above |
| ci95_high | as in the package or producer table named above |
| median_exact | as in the package or producer table named above |
| ci95_low_exact | as in the package or producer table named above |
| ci95_high_exact | as in the package or producer table named above |
| n_pos | as in the package or producer table named above |
| n_neg | as in the package or producer table named above |
| n_zero | as in the package or producer table named above |
| n_nan | as in the package or producer table named above |

## donor_recount.tsv

sha256: 5510b874ec3fa9559597edbf8a9ef430e62a5b9d72e612d941f4a72a3ead7e0e  
rows: 1,890  
produced by: healthy_stool_donor_entrapment  
used for v28d lines: 727-729  
Source Data package: per-seed records behind Supplementary_Figure_05/donor_entrapment_seed_means.tsv

FDRBench explicit-peptide entrapment, one row per target, seed, length and q threshold.

| column | meaning |
|---|---|
| target | acquisition searched |
| donor | acquisition the database was selected from |
| tier | T1_within_group (donor-selected) or T0_self |
| seed | random seed of the control |
| group | sampling group: healthy-stool design or CapScan region or stool; blank for PREDICT (visit codes removed) |
| min_length | lower peptide length bound (residues) |
| q_threshold | peptide q threshold |
| n_target | T |
| n_entrapment | E |
| partner_undiscovered | E0 |
| entrapment_ranks_higher | B |
| exact_score_ties | entrapments with exactly their target's score |
| paired_numerator | E + E0 + 2B |
| lower_fraction | E / (T + E) |
| combined_estimate | 2E / (T + E) |
| paired_estimate | (E + E0 + 2B) / (T + E) |
| paired_ties_upper | paired estimate counting ties as B |

## fusion_campi_acquisition_components.tsv

sha256: fe3269de4c6bef0ea1196e89e5ce9f8eacca21900651113608f0edb50b191114  
rows: 450  
produced by: campi_fusion_375m  
used for v28d lines: 715  
Source Data package: Supplementary_Figure_05/fusion_acquisition_components.tsv, rows with cohort CAMPI_MOCK (values unchanged)

CAMPI fusion control, one row per acquisition, seed, q threshold and length window.

| column | meaning |
|---|---|
| q_threshold | peptide q threshold |
| min_residues | lower peptide length bound (residues) |
| max_residues | upper peptide length bound (residues) |
| n_canonical | N, accepted canonical sequences |
| n_target_only | accepted in targets only |
| n_entrapment_only | accepted in controls only |
| n_shared | accepted in both |
| n_unpaired_accepted | U, accepted outside the audited one-to-one paired subset |
| n_paired_entrapment | E |
| n_entrapment_alone | E0 |
| n_entrapment_q_strict_better | entrapments with strictly better q than their accepted target |
| n_entrapment_q_tied | tied q (B = strict + tied) |
| n_entrapment_score_strict_better | by score (sensitivity) |
| n_entrapment_score_tied | tied score |
| n_q_score_order_discordant | pairs whose q and score order disagree |
| paired_numerator | E + E0 + 2B |
| score_order_numerator_sensitivity | numerator with score ordering |
| paired_diagnostic | (paired_numerator + U) / N |
| score_order_diagnostic_sensitivity | same with score ordering |
| paired_component_fraction | share of the diagnostic |
| unpaired_component_fraction | share |
| q_tie_component_fraction | share |
| n_multiple_peptidoforms | sequences with several peptidoforms |
| n_unknown_shared_target_control_sequence | U reason |
| n_unknown_partner_shared | U reason |
| n_unknown_nonbijective_mapping | U reason |
| n_unknown_outside_audited_digest | U reason |
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| seed | random seed of the control |

## fusion_decomposition_current375.tsv

sha256: 2a134617f73131d75914219516723e6bae7da3253b0de9397ec12beb3aae8825  
rows: 6  
produced by: campi_fusion_375m  
used for v28d lines: 715  
Source Data package: Supplementary_Figure_05/fusion_decomposition_current375.tsv (byte copy; compared, agrees)

Producer's pooled decomposition at q <= 0.01.

| column | meaning |
|---|---|
| seed | random seed of the control |
| min | lower length |
| max | upper length |
| n | N |
| E | E |
| E0 | E0 |
| B | B |
| q_ties | tied part of B |
| paired_numerator | E + E0 + 2B |
| U | U |
| paired_pct | percent |
| unresolved_pct | percent |
| total_pct | percent |
| unresolved_share_pct | 100 x U / (E + E0 + 2B + U) |

## hstool_pooling_per_acquisition.tsv

sha256: 474de2e78bbfc37fcf132930b688069df06909699b6f7cb9a0c12bf0c42222f8  
rows: 360  
produced by: healthy_stool_pooling_controls  
used for v28d lines: 743  
Source Data package: Figure3/per_acquisition.tsv (byte copy)

Full-route control acceptance in 12 technical replicates searched against own, all-12 and other-11 databases.

| column | meaning |
|---|---|
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| q_threshold | peptide q threshold |
| min_residues | lower peptide length bound (residues) |
| n_canonical | accepted canonical peptides |
| n_target_only | assigned to ordinary targets only |
| n_shared | shared |
| n_shuffle_only | assigned only to shuffled controls |
| shuffle_only_fraction | as in the package or producer table named above |
| n_plant_only | plant-only |
| plant_only_fraction | as in the package or producer table named above |
| n_thermophilic_only | thermophile-only |
| thermophilic_only_fraction | as in the package or producer table named above |
| shuffle_only_not_in_original_reservoir | shuffle-only and absent from the original reservoir (used for line 743) |
| shuffle_only_not_in_original_reservoir_fraction | as in the package or producer table named above |
| plant_only_not_in_original_reservoir | as in the package or producer table named above |
| plant_only_not_in_original_reservoir_fraction | as in the package or producer table named above |
| thermophilic_only_not_in_original_reservoir | as in the package or producer table named above |
| thermophilic_only_not_in_original_reservoir_fraction | as in the package or producer table named above |
| arm | search or control arm |

## human_reference.tsv

sha256: 055fcdc1b40f638b5118aac801fa3579e029597bd6d1ab293f963329cda25fc4  
rows: 1  
produced by: cohort_database_sizes  
used for v28d lines: 673  
Source Data package: not in the package

The reviewed canonical human Swiss-Prot FASTA used as the size reference.

| column | meaning |
|---|---|
| label | description recorded by the producer |
| sha256 | FASTA sha256 |
| entries | FASTA entries |
| residues | residues |
| unique_exact_sequences | distinct sequences |
| isoform_accessions | isoform-suffixed accessions |
| nonhuman_headers | non-human headers |
| nonreviewed_headers | TrEMBL headers |

## joint_counts_per_acquisition.tsv

sha256: 8076dd95e23525655760ab566cc8847d906843feb1f4cf55dd47b4ca9d2c5829  
rows: 116  
produced by: healthy_stool_joint_search_q  
used for v28d lines: 765  
Source Data package: not in the package

Accepted peptides per acquisition under Sage q and the D/T estimator applied per acquisition (file_) or once over all 58 (glob_).

| column | meaning |
|---|---|
| arm | search or control arm |
| filename | mzML name |
| psms | PSM rows |
| sage | accepted with Sage peptide q |
| file_plain | per-acquisition D/T |
| file_plus1 | per-acquisition (D+1)/T |
| negctl_file_plus1 | label-permuted control with (D+1)/T |
| glob_plain | D/T once over all acquisitions |
| glob_plus1 | (D+1)/T once over all |

## joint_fileq_coverage.tsv

sha256: 4e3dc4be9774a8809f42fddff16fe1f3bfc891206d6df53c91694b7922623361  
rows: 58  
produced by: healthy_stool_joint_search_q  
used for v28d lines: 765  
Source Data package: not in the package

Per-acquisition confidence on the joint search, with the label-permutation control.

| column | meaning |
|---|---|
| sample | acquisition (raw-file stem) |
| arm | search or control arm |
| mode | fileq |
| identified_canonical_peptides | accepted |
| positive_source_cells | positive LFQ cells |
| same_file_id_source_cells | cells with same-acquisition identification |
| positive_without_same_file_id | transfer-only cells |
| psm_rows | PSMs |
| accepted_rank1_target_psms | accepted PSMs |
| dropped_length | outside 5-50 |
| negctl_permuted_label_accepted | peptides accepted after one within-acquisition label permutation (D/T) |
| sage_q_accepted | accepted with Sage q |

## label_validity_presence.tsv

sha256: 7f65028e3e94872c5a0321456e58c8ed0e3292be38673942476eec7c8ec8df91  
rows: 13,441  
produced by: campi_fusion_375m  
used for v28d lines: 721  
Source Data package: Supplementary_Methods/label_validity/PRESENCE.tsv (byte copy)

Presence of each unique accepted entrapment-only sequence (ENTRAP), its composition-preserving shuffle (SHUFCTL) and sampled positive controls (REALPOS) in the 375,333,654-protein reservoir.

| column | meaning |
|---|---|
| arm | search or control arm |
| peptide | sequence (I/L-equivalent) |
| present | 1 if an exact substring of a reservoir protein |

## label_validity_result.tsv

sha256: 452a7706ff0d0c92cf6363ec78168e5fa04bf0005de54056af3c6ee8ec3760e2  
rows: 1  
produced by: campi_fusion_375m  
used for v28d lines: 721  
Source Data package: Supplementary_Methods/label_validity/RESULT.json (flattened)

Reservoir size and positive-control recovery.

| column | meaning |
|---|---|
| reservoir_sequences | proteins scanned |
| positive_total | positive controls |
| positive_recovered | positives found |
| q | q threshold |
| status | producer status |

## microbial_group_pair_summary.tsv

sha256: 776895eadb4f32753d6521237ba2322ced6fa7831ed83d33218c274cac3caf28  
rows: 36  
produced by: healthy_stool_joint_search_q  
used for v28d lines: 769  
Source Data package: Supplementary_Figure_07/per_acquisition_confidence/microbial_group_pair_summary.tsv (byte copy; compared, agrees)

Medians across pairs.

| column | meaning |
|---|---|
| definition | as in the package or producer table named above |
| comparison | as in the package or producer table named above |
| method | as in the package or producer table named above |
| tier | as in the package or producer table named above |
| pairs | as in the package or producer table named above |
| common_rows | as in the package or producer table named above |
| log2_pearson | as in the package or producer table named above |
| median_pair_cv_pct | as in the package or producer table named above |

## microbial_group_pairs.tsv

sha256: e8d58907fb5564149344d76b09e6dbe8b3dd3062e83084ed14da083143471329  
rows: 3,078  
produced by: healthy_stool_joint_search_q  
used for v28d lines: 769  
Source Data package: Supplementary_Figure_07/per_acquisition_confidence/microbial_group_pairs.tsv (byte copy)

Per-pair agreement of microbial protein groups.

| column | meaning |
|---|---|
| definition | matched or own_rows |
| comparison | pair of strategies compared |
| method | strategy scored |
| tier | TechnicalReplicates or WorkflowReplicates (independent preparations) |
| sample_a | acquisition |
| sample_b | acquisition |
| common_rows | rows used |
| log2_pearson | Pearson r of log2 intensities |
| median_pair_cv_pct | median over rows of population SD / mean x 100 |

## minlen7_benchmark_per_acquisition.tsv

sha256: 07ca109f6eeca2115b403da608fe1472c41b36c87921f87eebfdef42cb16c4eb  
rows: 1,773  
produced by: minlen7_search_sensitivity  
used for v28d lines: 777  
Source Data package: Supplementary_Figure_04/LENGTH_SENSITIVITY_source_data.xlsx sheet a_b_benchmark_per_acquisition (exported)

Matched benchmark per acquisition at baseline 5-50, nested 9-50 and fresh minimum-length-7 searches.

| column | meaning |
|---|---|
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| arm | search or control arm |
| comparator | MG |
| sample | acquisition (raw-file stem) |
| population | baseline5, nested9 or fresh7 |
| arm_peptides | FastaLake peptides |
| comparator_peptides | metagenome peptides |
| shared | both |
| gain_pct | 100 (arm/comparator - 1) |
| comparator_identity_recovery_pct | 100 shared/comparator |

## representative_presence.tsv

sha256: a67e5fa4085d91d3d3384ae47e9e1e5bcc87f86c09e00d39d5fc87737654d5b7  
rows: 58  
produced by: peptide_centric_downstream  
used for v28d lines: 751  
Source Data package: Figure3/representative_presence.tsv (byte copy)

Quantified features assigned to a representative absent from the acquisition's own FASTA.

| column | meaning |
|---|---|
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| retained_lfq_features | quantified features |
| outside_local_fasta_features | assigned to a representative absent from own FASTA |
| outside_local_members_features | no member in own FASTA |
| outside_local_fasta_fraction | fraction |
| positive_intensity | summed intensity |
| outside_local_fasta_intensity | intensity of the outside features |

## scramble_end_to_end.tsv

sha256: 844f92521718387c35fea524b5cc8502cac8847e0064dfea04c27637293ce865  
rows: 198  
produced by: scramble_end_to_end_workflow  
used for v28d lines: 667  
Source Data package: Review_v16/Terminal_scramble/scramble_end_to_end.tsv (byte copy)

One authentic and ten shuffled full workflows on the representative healthy-stool acquisition, accepted canonical peptides by gate, q and length.

| column | meaning |
|---|---|
| index | task index of the producing run |
| arm | search or control arm |
| seed | random seed of the control |
| selected_proteins | proteins retrieved by lookup (candidates) |
| parsimony_proteins | proteins in the selected FASTA after selection parsimony |
| gate | acceptance gate: q_only, or q plus PSM PEP <= 0.01 (pep01) or <= 0.001 (pep001) |
| q | peptide q threshold |
| min_length | lower peptide length bound (residues) |
| max_length | upper peptide length bound (residues) |
| n_peptides | accepted canonical target peptides |

## scramble_lookup_comparison.tsv

sha256: fd1fb867c1ef972b08b403ee1fc4e23fd7b590ef0a99904303a3fe04dd1bc052  
rows: 210  
produced by: scramble_end_to_end_workflow  
used for v28d lines: 667  
Source Data package: Review_v16/Terminal_scramble/scramble_lookup_comparison.tsv (byte copy)

Length-specific lookup outcome of terminal-preserving and all-position scrambled queries, per seed.

| column | meaning |
|---|---|
| arm | search or control arm |
| seed | random seed of the control |
| length | query length (residues) |
| distinct_queries_terminal | distinct terminal-preserving queries (denominator) |
| matched_queries | terminal-preserving queries matching at least one reservoir protein |
| match_fraction | matched_queries / distinct_queries_terminal |
| distinct_queries_all | distinct all-position queries |
| matching_queries | all-position queries matching |
| fraction_matched | matching_queries / distinct_queries_all |
| difference_pp | difference of the two fractions, percentage points |

## scramble_lookup_summary.tsv

sha256: ec445c5ba325ff2ccd6ccc67901a9b2619cd0ff1a867f12c18daae4c0e4c40f5  
rows: 42  
produced by: scramble_end_to_end_workflow  
used for v28d lines: 667  
Source Data package: Review_v16/Terminal_scramble/scramble_lookup_summary.tsv (byte copy; used only to compare, values agree)

Minimum, median and maximum match fractions across the five seeds, per length.

| column | meaning |
|---|---|
| length | query length |
| terminal_min | fraction |
| terminal_median | fraction |
| terminal_max | fraction |
| all_min | fraction |
| all_median | fraction |
| all_max | fraction |

## selection_grid_acquisitions.tsv

sha256: 8b472997023e29d4e24b5a9f7a06d9aceda6f45fc4c45391825c075d85aa3425  
rows: 1,102  
produced by: healthy_stool_selection_sweep  
used for v28d lines: 679  
Source Data package: Supplementary_Figure_02/selection_grid_acquisitions.tsv (byte copy)

Selection-sensitivity grid, one row per acquisition and setting (58 x 19).

| column | meaning |
|---|---|
| arm | search or control arm |
| kind | length, score or reference (production top-30%, minimum length 9) |
| min_length | lower peptide length bound (residues) |
| top_percent | fraction of top-scoring predictions used as lookup queries |
| sample | acquisition (raw-file stem) |
| s3_seeds | lookup queries |
| s3_match_rate | fraction of queries matching |
| s3_proteins | proteins retrieved |
| razor_selected | proteins in the selected FASTA |
| razor_unique_peptides | peptides unique to one selected protein |
| razor_mapped_peptides | peptides mapped to selected proteins |
| psms_q01 | PSMs at q <= 0.01 |
| decoy_frac_q01 | decoy fraction at q <= 0.01 |
| peptides_q01 | accepted peptides (modifications ignored, I/L not merged) |
| proteins_q01 | proteins at q <= 0.01 |
| scored_candidates_median | median candidates scored per spectrum |
| status | run status |

## simulation_all_outcomes.tsv

sha256: e6427d549e37a18d60d26447c962c274f08c2b41a58594a66edbadc7029067a1  
rows: 288  
produced by: restricted_reference_simulation  
used for v28d lines: 733  
Source Data package: the 72 ordinary rows equal Supplementary_Figure_05/simulation_mixture_outcomes.tsv; the 216 others are not in the package

All 288 simulation outcomes (24 mixtures x 12 arms), in the original row order used by the bootstrap.

| column | meaning |
|---|---|
| arm | search or control arm |
| discoveries | target canonical peptides at q <= 0.01 |
| false_discoveries | discoveries not among the generating sequences |
| realised_FDP | false_discoveries / max(discoveries, 1) |
| decoy_winners | decoys passing |
| empty_discovery_set | no discoveries |
| unit_id | mixture |
| regime | baseline or stress |

## simulation_design.tsv

sha256: 868c66c680094cafd8afde1bc229a6c21ac14bb0f4988eed7713624d786a4e21  
rows: 8  
produced by: restricted_reference_simulation  
used for v28d lines: 733  
Source Data package: not in the package (design constants)

Design constants read from the simulation code lines and the reference FASTA.

| column | meaning |
|---|---|
| parameter | constant |
| regime | baseline, stress or all |
| value | value |
| origin_file | code file |
| origin_line | line in that file |
| origin_sha256 | sha256 of that file |

## simulation_uncertainty.tsv

sha256: f267eee37eaa0cf29e122b82a27914abdb309a8edc31a8d4fd9c397d20b0e0ec  
rows: 8  
produced by: restricted_reference_simulation  
used for v28d lines: 733  
Source Data package: Supplementary_Figure_05/simulation_uncertainty.tsv (byte copy; compared, agrees)

Producer's means and bootstrap intervals.

| column | meaning |
|---|---|
| regime | baseline or stress |
| endpoint | arm or paired difference |
| n_units | mixtures |
| mean | fraction |
| bootstrap_low | 2.5% quantile |
| bootstrap_high | 97.5% quantile |
| t_low | t interval (not printed) |
| t_high | t interval |
| total_discoveries | sum |
| total_false | sum |
| pooled_FDP | total_false / total_discoveries |

## split_half_acquisitions_710.tsv

sha256: ea5c13156581a3bcd6b6ea49a9267e7a236578d4e27d61e9348725ded4d24abb  
rows: 1,420  
produced by: cohort_spectrum_split_selection  
used for v28d lines: 755-761, 773  
Source Data package: Review_v16/Native_parity/FULL_SOURCE_DATA.tsv holds the same counts in long format

Spectrum-split counts per acquisition and length view (710 acquisitions x 2). Raw PREDICT participant identifiers replaced by the Source Data package pseudonyms PREDICT_P###; specimen column removed.

| column | meaning |
|---|---|
| index | task index of the producing run |
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| min_length | lower peptide length bound (residues) |
| full | accepted peptides, all-prediction database |
| self | same-half (odd) database |
| disjoint | even-scan database |
| D_over_S | disjoint/self |
| D_over_F | disjoint/full |
| S_over_F | self/full |
| lost_from_self | in self not disjoint |
| gained_vs_self | in disjoint not self |
| shared_with_self | in both |
| fraction_self_identities_retained | shared/self |
| disjoint_fraction_in_full | fraction of disjoint in full |
| even_predictions | predictions from even scans |
| odd_predictions | from odd scans |
| prediction_ratio | even/odd |
| even_selected_proteins | proteins selected from even |
| odd_selected_proteins | from odd |
| selected_protein_ratio | even/odd |
| conversion_recovered | one of five reconverted PREDICT inputs |
| participant | PREDICT_P### pseudonym (Source Data package scheme), CAPSCAN_<number>, blank if not recorded |

## split_half_acquisitions_726.tsv

sha256: 15efb29be1e256570fe758f3c54b58c9ecc868e4b7e3f1221d28fd9c25ce02ef  
rows: 1,452  
produced by: cohort_spectrum_split_selection  
used for v28d lines: 755-761 (v28f alternative)  
Source Data package: not in the package (CapScan re-run on 347 acquisitions)

As the 710 table with CapScan rows from the 347-acquisition re-run.

| column | meaning |
|---|---|
| index | task index of the producing run |
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| min_length | lower peptide length bound (residues) |
| full | as in the package or producer table named above |
| self | as in the package or producer table named above |
| disjoint | as in the package or producer table named above |
| D_over_S | as in the package or producer table named above |
| D_over_F | as in the package or producer table named above |
| S_over_F | as in the package or producer table named above |
| lost_from_self | as in the package or producer table named above |
| gained_vs_self | as in the package or producer table named above |
| shared_with_self | as in the package or producer table named above |
| fraction_self_identities_retained | as in the package or producer table named above |
| disjoint_fraction_in_full | as in the package or producer table named above |
| even_predictions | as in the package or producer table named above |
| odd_predictions | as in the package or producer table named above |
| prediction_ratio | as in the package or producer table named above |
| even_selected_proteins | as in the package or producer table named above |
| odd_selected_proteins | as in the package or producer table named above |
| selected_protein_ratio | as in the package or producer table named above |
| conversion_recovered | as in the package or producer table named above |
| participant | as in the package or producer table named above |

## split_half_diagnostic_sensitivity.tsv

sha256: 0df389c1b992aadc96c288bcbba403503bce8ae69f4ea960ed6cef473c16e3f7  
rows: 24  
produced by: cohort_spectrum_split_selection  
used for v28d lines: 761  
Source Data package: Review_v16/Native_parity/DIAGNOSTIC_SENSITIVITY.tsv (byte copy)

Cohort medians excluding fallback-affected acquisitions.

| column | meaning |
|---|---|
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| min_length | lower peptide length bound (residues) |
| subset | all, exclude_any_model_fallback, prediction_counts_within_2pct |
| n | acquisitions |
| median_D_over_S | median |
| median_D_over_F | median |

## split_half_model_fallbacks.tsv

sha256: 3c1706dee587086d7dc8d2b5321a1f13a908b88156e2bb67ee44f4280281e055  
rows: 13  
produced by: cohort_spectrum_split_selection  
used for v28d lines: 761  
Source Data package: SEARCH_MODEL_REVIEW.json of the producer (flattened)

Sage arm searches that used heuristic scoring after discriminant-fit failure.

| column | meaning |
|---|---|
| index | task index of the producing run |
| cohort | cohort: CAMPI or CAMPI_MOCK, HSTOOL (healthy stool), PREDICT, CAPSCAN / CapScan |
| sample | acquisition (raw-file stem) |
| arm | search or control arm |

## split_half_zero_yield.tsv

sha256: 4951b4f637dc11291d9e93dd67dc7009d6c45c09d20c58b7d9acd98402397f5e  
rows: 3  
produced by: cohort_spectrum_split_selection  
used for v28d lines: 761  
Source Data package: SEARCH_MODEL_REVIEW.json of the producer (flattened)

The CapScan acquisition with a zero-yield disjoint arm.

| column | meaning |
|---|---|
| index | task index of the producing run |
| sample | acquisition (raw-file stem) |
| arm | search or control arm |
| rank_one_targets | rank-one target PSMs |
| rank_one_decoys | rank-one decoys |
| minimum_target_peptide_q | smallest target peptide q |
| accepted_peptides | accepted at q <= 0.01 |

## technical_agreement_summary.tsv

sha256: c9d59a24319a1106ab0e1e7895b50133040d78a28faca01918c1895bf7632751  
rows: 28  
produced by: healthy_stool_common_fasta  
used for v28d lines: 769  
Source Data package: Supplementary_Figure_07/technical_agreement_summary.tsv (byte copy)

Leading annotated species agreement.

| column | meaning |
|---|---|
| method | strategy |
| tier | replicate tier |
| threshold | both_resolved or both_strict |
| original_agreement | fraction, acquisition-specific |
| test_agreement | fraction, tested strategy |

## transfer_pool_comparisons.tsv

sha256: a7618438c577a20b6037bf801b4daa89c181dcafb724066cfada229b06035947  
rows: 420  
produced by: healthy_stool_database_transfer  
used for v28d lines: 739-743  
Source Data package: per-acquisition table behind Figure 3 transfer panels (not a package sheet)

Transfer and pooling comparisons per target acquisition, arm and length.

| column | meaning |
|---|---|
| sample | acquisition (raw-file stem) |
| kind | transfer or pool |
| group | sampling group: healthy-stool design or CapScan region or stool; blank for PREDICT (visit codes removed) |
| arm | search or control arm |
| min_length | lower peptide length bound (residues) |
| own | own-database peptides A |
| peptides | test-arm peptides B |
| overlap | |A n B| |
| lost | in A not B |
| gained | in B not A |
| count_retention_percent | 100 |B|/|A| |
| own_set_recovery_percent | 100 |A n B|/|A| |
| own_sequence_present | own peptides present in the test FASTA |
| own_digest_compatible | own peptides digest-compatible there |
| own_digest_compatible_unidentified | compatible but not identified |
| database_accessions | accessions in the test FASTA |
| donor | acquisition the database came from |

## union_contributions.tsv

sha256: 5eaa7eb32db23fa8306754567589d3f9550b8a1a5feb5650ae8f3604039dfeb9  
rows: 1,042  
produced by: cohort_selected_fasta_unions  
used for v28d lines: 751  
Source Data package: not in the package; byte copy of source_data/methods/union_contributions.tsv (unit methods)

Per selected FASTA, the exact sequences not in any earlier FASTA of the roster; a roster's union size is the column sum.

| column | meaning |
|---|---|
| roster | HSTOOL58, PREDICT306, CAPSCAN331 (v28d) or CAPSCAN347 (v28f) |
| order | position in acquisition-name order |
| acquisition | acquisition |
| records | records in the FASTA |
| new_unique_sequences | sequences new to the union (I/L distinct) |
| fasta_sha256 | FASTA sha256 |


## adjacent_swap_lookup.tsv

sha256: 26c380324bb27f9b02ea516b11460ea262447d1acb75b1a4dcecc28c2d142174  
rows: 2  
produced by: bag_mass_lookup_search (adjacent_swap/, Slurm job 5954787)  
used for v30 paragraph: 344 (Supplementary Methods 6, adjacent-residue-swap sensitivity)  
Source Data package: not yet in the package (byte copy of the producer's RESULT_TABLE.tsv)

Adjacent-swap lookup of the representative healthy-stool acquisition against the 233,497,702-entry reservoir. Run A searches the original queries; run B the originals plus every single adjacent-swap variant. REAL = authentic retained queries; SIM = composition-preserving scramble, seed 20260913.

| column | meaning |
|---|---|
| arm | REAL (authentic queries) or SIM (scrambled control) |
| originals | distinct I/L-collapsed retained queries |
| swap_only_variants | swap variants not identical to an original query |
| (a) proteins_by_originals | reservoir proteins matched by the originals (run A) |
| (b) additional_by_swaps_only | further proteins matched only by swap variants |
| b/a | (b) divided by (a) |
| originals_without_exact_hit | original queries with no exact reservoir match |
| rescued_by_swap | of those, queries with at least one matching swap variant |
| (c) swap_only_in_razor_cover | swap-only candidate proteins kept by selection parsimony in run B (REAL only; - for SIM) |
| razor_cover_size_B | proteins in the run-B parsimony cover (REAL only) |
| razor_cover_size_A | proteins in the run-A parsimony cover (REAL only) |

## crossfit_pep_summary.tsv

sha256: 6ee68ef36f14689a9797718b9b3ff0628748d39656585215461fef255959a277  
rows: 48  
produced by: capscan_crossfit_10fold (cf_pep.py, cf_analyze.py; Slurm jobs 6057990, 6058010)  
used for v30 paragraphs: 125, 326  
Source Data package: not yet in the package (byte copy of CF_PEP_SUMMARY.tsv)

Cross-fitted FastaLake gain over the matched metagenome search with and without an added PEP gate, per cohort stratum. Comparisons: `gain_cf_<option>` = 100 (FL/MG - 1); `cf_vs_full_<option>` = 100 (FL/FULL - 1), FULL being the all-prediction arm. Options: baseline (pooled LDA rescoring and picked-peptide q <= 0.01), opt_a (baseline plus the per-fold Sage posterior_error <= -2, i.e. PEP <= 0.01, on the same PSM; unset PEP excluded), opt_b (baseline plus PEP <= 0.01 re-estimated from the pooled discriminant scores with the Sage PEP model). Medians with participant-cluster bootstrap 95% intervals (10,000 draws, seed 20261004).

| column | meaning |
|---|---|
| stratum | PREDICT, CAPSCAN_all, CAPSCAN_stool, CAPSCAN_intestinal (intestinal fluid) |
| min_length | lower peptide length bound (5 or 9 residues; upper bound 50) |
| comparison | `gain_cf_<option>` or `cf_vs_full_<option>` (see above) |
| n_acquisitions | acquisitions in the stratum |
| n_participants | participants (patients) in the stratum |
| median | median percentage gain |
| ci95_low, ci95_high | bootstrap 95% interval |
| n_pos, n_neg, n_zero, n_nan | acquisitions with a positive, negative, zero or undefined gain |

## crossfit_pep_summary_by_region.tsv

sha256: 73978e823366365636188d85f98c9d8b1a943d91a04440d7e9931a791e7065e9  
rows: 60  
produced by: capscan_crossfit_10fold (cf_analyze.py, Slurm job 6058010)  
used for v30 paragraph: 326  
Source Data package: not yet in the package (byte copy of CF_PEP_SUMMARY_BY_REGION.tsv)

As crossfit_pep_summary.tsv, with CapScan split by sampling region (CAPSCAN_duodenum, CAPSCAN_jejunum, CAPSCAN_ileum, CAPSCAN_colon, CAPSCAN_stool). Same columns.

## crossfit_pep_per_acquisition.tsv

sha256: 457d13cc117e555197d02dd4bbb29cad808f617772324298c6b5c80766714e35  
rows: 591  
produced by: capscan_crossfit_10fold (cf_analyze.py, job 6058010), pseudonymised by export_pep_source_data.py  
used for v30 paragraphs: 125, 326  
Source Data package: not yet in the package. Differs from the producer's CF_PEP_PER_ACQUISITION.tsv only in `subject_id` (raw identifier replaced by the pseudonym of crossfit_per_acquisition.tsv) and `group` (blank for PREDICT)

| column | meaning |
|---|---|
| i | row of the 578-acquisition matched roster |
| cohort | PREDICT or CAPSCAN |
| sample | acquisition (raw-file stem) |
| group | CapScan sampling region; blank for PREDICT |
| subject_id | pseudonymised participant (PREDICT_P###, CapScan participant id) |
| baseline_equal_all | True if the recomputed baseline sets equal the frozen pooled sets in all three arms |
| fallback | arms whose pooled rescoring fell back to Sage's heuristic score (empty if none) |
| `<ARM>_<option>_<L>` | accepted canonical peptides of length L-50 (L = 5 or 9) in arm FL, FULL or MG under option baseline, opt_a or opt_b |
| `gain_cf_<option>_<L>` | 100 (FL/MG - 1) |
| `cf_vs_full_<option>_<L>` | 100 (FL/FULL - 1) |
