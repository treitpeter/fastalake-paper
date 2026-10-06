# Source data: Supplementary Methods 1-5

Origins, package comparison and transformations: `text_statistics/supplementary_methods_1_5/PROVENANCE.md`. All tables were written by `text_statistics/supplementary_methods_1_5/prepare_source_data.py`. Tables read from other units (not repeated here): `Fig1/stage1_stats.json`, `Fig2/alphanovo_accuracy_per_acquisition.tsv`, `Fig2/reservoir_protein_presence_by_proteome.tsv`, `SFig02/bags_bounded48_per_acquisition.tsv`, `SFig02/bags_bounded48_orderings.tsv`.

Boolean columns hold `True` or `False`. Sample names are instrument file names. No clinical covariate is present; PREDICT patients appear only as pseudonyms `U01`-`U77`.

## parser_audit_acquisitions_726.tsv
sha256: 28c39e5eb91f4472730e2777fcc12b96b8a3641a482f13dc193698b162443248
rows: 726 (one row per AlphaNovo prediction file; indices 0-709 are the v28d audit, 710-725 the 16 participant-15 CapScan files)
produced by: prediction_parser_audit

| column | unit | meaning |
|---|---|---|
| index | id | task index |
| cohort | id | CAMPI, HSTOOL, PREDICT or CAPSCAN |
| sample | id | acquisition |
| prediction_sha256 | id | sha256 of the prediction file |
| prediction_bytes | bytes | file size |
| prior_hash_matches | boolean | file hash equals an earlier recorded hash; empty if none was recorded |
| csv_errors, unequal_field_count | count of rows | CSV parse errors and rows with a wrong field count |
| eligible_sequence_parser_differences | count of rows | lookup-eligible rows whose parsed sequence differs between parser versions |
| missing_score, nonfinite_score, unparseable_score | count of rows | score problems |
| nonpositive_score, positive_finite_score | count of rows | rows by score sign |
| nonzero_beam_rank | count of rows | rows that are not the top-ranked beam |
| rows | count of rows | prediction rows in the file |
| sequence_parser_differences | count of rows | rows whose parsed sequence differs between parser versions (any row) |
| top30_sequence_parser_differences | count of rows | as above among rows retained by the top-30% rule |
| lookup_rows_exact | boolean | archived and current parser give identical lookup rows and scores |
| selection_exact | boolean | archived and current parser give identical selection evidence |
| top30_unique | count of sequences | distinct top-30% lookup queries |
| selection_unique | count of sequences | distinct selection-evidence peptides |
| selected_records | count of proteins | records in the selected FASTA |
| lookup_only_queries | count of sequences | top-30% lookup queries absent from selection evidence |
| seconds | s | run time of the check |

## search_input_completeness.tsv
sha256: 6422821e85bd7c8be6d56bd0da8991d5b2e12c2d522c8a09c88ec9fd88304966
rows: 711 (one row per Sage search: 306 PREDICT and 347 CapScan public-sequence searches, 58 healthy-stool searches)
produced by: sage_input_completeness

| column | unit | meaning |
|---|---|---|
| lineage | id | `primary_public_current` (public-sequence searches) or `healthy_stool_58` |
| cohort | id | PREDICT, CAPSCAN or HSTOOL |
| sample | id | acquisition |
| repaired | boolean | the mzML is one of the five reconverted PREDICT files |
| flag | boolean | flagged by the original check (missing closing tag or coverage below 0.9) |
| closing_tag | boolean | the mzML ends with its closing tag |
| nul_in_tail | boolean | NUL bytes in the last block of the mzML |
| declared_spectra | count of spectra | spectrumList count declared in the mzML |
| max_scan | scan number | highest scan scored by Sage |
| max_scan_fraction | dimensionless | max_scan / declared_spectra, rounded to 4 decimals by the producer (the notebook recomputes it) |
| rows | count of PSMs | rows in results.sage.tsv |
| target_peptides_q01_modified | count of peptides | target peptides at peptide q <= 0.01, modified forms distinct |
| bytes | bytes | mzML size |
| software | id | conversion software recorded in the mzML |
| sage_version | id | version reported by Sage |
| n_mzml | count | mzML files in the search |
| error | text | error message of the check; empty if none |

## hstool97_stage1_stats.tsv
sha256: 93012c9a2b1e686e1c056bd2256a98d21d8390dbfbca6707646ed96c666e5c26
rows: 13 (one row per statistic of the original 97-acquisition study-level lookup)
produced by: healthy_stool_97_study_lookup

| column | unit | meaning |
|---|---|---|
| key | id | statistic name as written by the lookup; `proteins_matched` = study-supported proteins, `proteins_searched` = reservoir proteins, `total_peptides` = pooled distinct queries |
| value | as named | value; `database` is the reservoir file name |

## hstool97_pooled_prediction_files.tsv
sha256: a491b3a9acd1c4219823d327a277532999fb5bc1cd7e2e04c373383c96025668
rows: 97 (one row per prediction file pooled by the original production lookup)
produced by: healthy_stool_97_study_lookup

| column | unit | meaning |
|---|---|---|
| prediction_file | id | `<acquisition>_predictions.csv` |
| bytes | bytes | file size |

## origin_rows_HSTOOL.tsv, origin_rows_PREDICT.tsv, origin_rows_CAPSCAN331.tsv, origin_rows_CAPSCAN347.tsv
sha256: 7ffd1119dec76ee115178a3b7c66d2e8ba1731bfc3b784c64e0c6b94873a7004 (HSTOOL); 7dab235331618ea4970b3366f2950d6f4d1a783d087bef2ec9b8a42beb54bfa5 (PREDICT); 3547fee389ef9a63f8c857b55f8345d3891210726f04910c1f7487fba576155d (CAPSCAN331); d5c3fdeb3230b499997e348afe40786a816c74a202b2b344e6e92574da6f133b (CAPSCAN347)
rows: 1,953; 6,086; 4,463; 4,499 (one row per annotated-origin class)
produced by: origin_annotation_share

| column | unit | meaning |
|---|---|---|
| row | id | row index in the producer's taxonomy table |
| rank | id | rank of the deepest shared lineage |
| state | id | `supported_lca`, `annotated_parent_lca` (some compatible proteins had no origin annotation), `root_only`, `unannotated_taxonomy` |
| lineage | id | GTDB-style lineage, `;`-separated |
| human_species | boolean | lineage ends in `s__Homo sapiens` |
| signal_pooled | LFQ intensity | intensity of features in this class, summed over all acquisitions |
| acquisitions | count of acquisitions | acquisitions pooled |

## mgnify_seed_ortholog_check_by_genome.tsv
sha256: 7f3e39ec9626d6fa75fbbd31b678663a50afe6b17f14aead924831624c869770
rows: 4,744 (one row per UHGG species-representative genome)
produced by: mgnify_seed_ortholog_check

| column | unit | meaning |
|---|---|---|
| genome | id | MGnify genome accession |
| rows | count of annotations | MGnify eggNOG annotation rows checked |
| no_alignment | count of annotations | DIAMOND found no alignment to the seed ortholog |
| fail_ratio_below_half | count of annotations | DIAMOND bitscore below half the reported bitscore (includes no_alignment) |
| fail_query_shorter_than_half_seed | count of annotations | failed rows whose protein is shorter than half the seed ortholog |
| fail_seed_equals_adjacent | count of annotations | failed rows whose seed is also the seed of the adjacent gene |
| pass_query_shorter_than_half_seed, pass_seed_equals_adjacent | count of annotations | the same among passing rows |

## annotation_mgnify_coverage_by_arm.tsv
sha256: ea585a0069f05eeec5a783204e9085a814d7fc3f09c5998ad38d69e52416b003
rows: 9 (one row per study arm; CapScan once for the v28d 331 roster and once for the 347 roster)
produced by: annotation_repair_coverage

| column | unit | meaning |
|---|---|---|
| arm | id | study arm |
| roster | id | `v28d` or `347` |
| acquisitions | count of acquisitions | acquisitions in the arm |
| compatible_parents | count of proteins | compatible proteins of accepted LFQ features |
| native_gene_sequence_verified | count of proteins | proteins with an MGnify annotation joined by gene identifier and verified sequence |
| native_full_sequence_annotation_transfer | count of proteins | proteins annotated through an identical full-length sequence |
| native_not_available | count of proteins | proteins without an MGnify annotation |

## campi_supplemental_annotation_queries.tsv
sha256: 96fc30bd22463fe881e7fb1ec76d64fec08b1218d12e0ff12f1b4f234317dbf5
rows: 17,935 (one row per CAMPI reference protein sent to the supplementary eggNOG-mapper run)
produced by: campi_supplemental_eggnog

| column | unit | meaning |
|---|---|---|
| accession | id | FastaLake sequence accession |
| sequence_sha256_IL | id | sha256 of the I/L-normalised sequence |
| sequence_verified | boolean (1/0) | sequence matched its source record |
| terminal_stop_removed | boolean (1/0) | a terminal stop character was removed |
| has_eggnog_record | boolean | eggNOG-mapper returned an annotation |

## campi_target_fasta_records.tsv
sha256: 4b089729d14ec5bf7147977272e5f6d5f05a523c4a4ebb5931fbab2be484d5eb
rows: 49,224 (one row per record: 29,673 reference, 19,551 MG)
produced by: supplied CAMPI files (no upstream analysis)

| column | unit | meaning |
|---|---|---|
| database | id | `reference` or `MG` |
| record | index | 0-based position in the FASTA |
| header_token | id | first token of the header |
| header_kind | id | `uniprot` (UniProt accession), `assembly_contig` (k-mer contig gene) or `other` (contaminant) |
| header_has_trailing_whitespace | boolean | header line ends in whitespace |
| length | residues | sequence length |
| sequence_sha256 | id | sha256 of the sequence |

## allocation_coverage_HSTOOL.tsv, allocation_coverage_PREDICT.tsv, allocation_coverage_CAPSCAN347.tsv
sha256: 2441b85054ba37b6786a459c633ae59d6dada7bc18f4315b9a8e9df0e2c4cec9 (HSTOOL); 70085379b026377ebd709a947ee45bb864ac8c5019d30bb0c116c5a1c97e28ee (PREDICT); 2c04660ce6631fe61e3274f649ef8768b2c3c761befaf93206177308c734730c (CAPSCAN347)
rows: 1,392; 7,344; 8,328 (one row per acquisition, level, profile and allocation rule)
produced by: functional_allocation_sensitivity

| column | unit | meaning |
|---|---|---|
| arm | id | study arm |
| level | id | KO, OG or Pathway |
| profile | id | `all` or `other_unresolved` (human-assigned signal withheld) |
| variant | id | allocation rule: `equal`, `unique`, `complete`, `unique_complete` |
| sample | id | acquisition |
| total_LFQ | LFQ intensity | total intensity of the acquisition |
| allocated | LFQ intensity | intensity allocated to functional terms |
| unassigned | LFQ intensity | intensity without an eligible term |
| human_withheld | LFQ intensity | human-assigned intensity withheld (0 in profile `all`) |

## allocation_capscan_own_stool_distances_347.tsv
sha256: fbebb380c85f4b6e7b417380f4749b693b5180b83235aabce93058fe730b0ae7
rows: 2,925 (13 settings x 15 x 15 participant pairs)
produced by: functional_allocation_sensitivity

| column | unit | meaning |
|---|---|---|
| variant | id | allocation rule, or `NEGATIVE_CONTROL_stool_rotated` |
| prevalence | fraction | OG prevalence filter across the 75 participant-region cells |
| features | count of OGs | OGs kept by the filter |
| intestinal_participant | id | participant (1-15) whose four intestinal-region profiles are compared |
| stool_participant | id | participant whose stool profile is compared |
| mean_hellinger_distance | dimensionless | mean Euclidean distance of square-root relative profiles over the four regions |

## allocation_predict_distances.tsv
sha256: 6631628d94357fc8203a209a1fe38862a3cff018c6a3d828e01110d00c1b233e
rows: 79,002 (27 distance matrices x 2,926 patient pairs)
produced by: functional_allocation_sensitivity

| column | unit | meaning |
|---|---|---|
| level | id | KO, OG or Pathway |
| layer | id | metaG, metaT, microbial_metaP, host_metaP or metabolites |
| variant | id | allocation rule for microbial_metaP (or the permuted negative control); `all` for the other layers |
| unit_i, unit_j | id | patient pseudonyms U01-U77 (i < j, in scipy pdist order) |
| distance | dimensionless | Euclidean distance between transformed profiles |

## allocation_hstool_ko_pairs.tsv
sha256: 4d750a489699f6cfdba75c2fb97a9410bd2cd3258ed76de7122e724f8935139a
rows: 2,844 (316 within-design acquisition pairs x 9 settings)
produced by: functional_allocation_sensitivity

| column | unit | meaning |
|---|---|---|
| profile | id | `all` or `other_unresolved` |
| variant | id | allocation rule, or `NEGATIVE_CONTROL_feature_shuffle` |
| design | id | TechnicalReplicates, WorkflowReplicates, OnePiece, Intra, Interindividual |
| i, j | id | the two acquisitions |
| r_per_pair | dimensionless | Pearson r of log2 KO intensity over KOs positive in both |
| cv_per_pair_pct | % | median coefficient of variation over those KOs |
| n_features_per_pair | count of KOs | KOs positive in both |

## allocation_fraction_medians_producer.tsv, allocation_capscan_own_stool_producer.tsv, allocation_predict_crossomics_producer.tsv, allocation_hstool_ko_precision_producer.tsv
sha256: e18fa519328dfa4981100b404e8119160c6da45a425a1649d82722fc8a3c6311; f220676ce7c4342325ad3b4c57adca41e879da87b6fa3a8b2cea12fef7ed69c7; a26ed8de60f2923da590c49f481371bd3d90be01df933aaf76ec4f9219d21b80; aba3ca25544b5003acf25ed2b10bb80b4c9831c86304e7260223d6eaa394062c
rows: 36; 13; 60; 45
produced by: functional_allocation_sensitivity

The producers' own result tables, copied byte for byte. The notebook compares its recomputed results with them (`stats.py:check_against_producers`); they are not used to compute any printed number. Columns as written by the producer: allocated-fraction medians (`arm, level, variant, allocated_fraction`); own-stool results for the 331 and 347 rosters (`variant, prevalence`, then `features, own_closer, mean_own, mean_other, mean_advantage, p` with suffix `_331` or `_347`, `participants`, `acquisitions_without_allocated_signal_*`); PREDICT rows (`variant, cohort, scope, level, layer_a, layer_b, prevalence, n_units, n_donors, rho, p`); healthy-stool medians (`design, r_per_pair, cv_per_pair_pct, n_features_per_pair, profile, variant`).

## header_audit.tsv
sha256: 7a66dc68c0fe38b5885d9412e82dc7a0c71388d9e58aec5deaa536ab3ea751ce
rows: 8 (one row per cohort and accession namespace; PREDICT 296 specimens, CapScan 250 specimens of the 282 matched acquisitions)
produced by: metagenome_header_audit (the two `full/<cohort>/header_audit.tsv` outputs of worker.py, concatenated with a cohort column; values unchanged)

Sequence identity is the SHA-256 of the exact amino-acid string (no I/L collapse); header identity is the first whitespace-delimited header token. Namespaces: `FL_public_as_searched` (public-sequence FastaLake FASTAs as searched), `MG_as_searched` (matched metagenome plus human Swiss-Prot FASTAs as searched, hashed accessions), `MG_original_prefixed_headers` (specimen assembly proteins under their original sample-prefixed contig names), `MG_assembly_local_headers` (the same with the sample prefix stripped, from `k141_`). Rates are empty where the denominator is zero.

| column | unit | meaning |
|---|---|---|
| cohort | id | PREDICT or CAPSCAN |
| namespace | id | accession namespace, see above |
| distinct_sequences | count of sequences | distinct sequences in the namespace across all specimens |
| sequences_shared_across_specimens | count of sequences | distinct sequences present in more than one specimen |
| shared_sequences_with_multiple_headers | count of sequences | shared sequences carrying more than one header |
| alias_rate_among_shared_sequences | fraction | shared_sequences_with_multiple_headers / sequences_shared_across_specimens |
| alias_rate_among_all_sequences | fraction | shared_sequences_with_multiple_headers / distinct_sequences |
| distinct_headers | count of headers | distinct header tokens |
| headers_shared_across_specimens | count of headers | header tokens present in more than one specimen |
| shared_headers_with_multiple_sequences | count of headers | shared header tokens naming more than one sequence |
| collision_rate_among_shared_headers | fraction | shared_headers_with_multiple_sequences / headers_shared_across_specimens |
| collision_rate_among_all_headers | fraction | shared_headers_with_multiple_sequences / distinct_headers |
