# Source data: Fig4

Tables read by `figures/Fig4/Fig4.ipynb`. They were cut from the analysis outputs by `figures/Fig4/prepare_source_data.py`; origins and checks against the published Source Data package are in `figures/Fig4/PROVENANCE.md`, input hashes in `EXTRACTION_RECEIPT.json`. PREDICT acquisitions are `PREDICT_A###` (numbered by sorted acquisition name over the 306 PREDICT acquisitions) and PREDICT patients `PREDICT_P###` (the pseudonyms of the Source Data package); no clinical covariate is in any table. CapScan participants keep their published numbers 1-15. Sampling sites of CapScan are the source labels duodenum, jejunum, ileum, colon (the manuscript's regions 1-4) and stool. Tables rewritten by the extraction store floats with 17 significant digits; `panels.py:load_table` reads them so that the values equal those the original producers used.

## crossfit_per_acquisition.tsv
sha256: 5bcdbe883a9ddb41582de92d11c7ad3bc739fa272b71241fcc86244cd2b9b094
rows: 578 (one row per acquisition with a matched metagenome: 296 PREDICT, 282 CapScan)
produced by: capscan_crossfit_10fold

| column | unit | meaning |
|---|---|---|
| crossfit_index | id | row index of the original cross-fit table (row order is the order the bootstrap and the plot jitter use) |
| cohort | id | PREDICT or CAPSCAN |
| acquisition_id | id | PREDICT_A### or the CapScan raw-file name |
| participant_id | id | PREDICT_P###, PREDICT_NOID_1-4 (acquisitions without a patient number; each its own bootstrap cluster) or CapScan participant 1-15 |
| participant_sort_rank | dimensionless | rank of the original participant identifier within the cohort; the bootstrap orders clusters by it, as the original analysis did |
| sampling_site | id | stool (all PREDICT) or CapScan duodenum, jejunum, ileum, colon, stool |
| lda_fallback | id | arms whose pooled rescoring fell back to the heuristic score (empty if none) |
| FL_pooled_5, MG_pooled_5 | count of peptides | unique peptides of 5-50 residues accepted (peptide q ≤ 0.01) after pooled rescoring of the 10 cross-fitted folds, FastaLake and metagenome search |
| frozen_fl_5, frozen_mg_5 | count of peptides | same-spectrum (production) FL and MG peptide counts, 5-50 residues |
| FL_pooled_9, MG_pooled_9, frozen_fl_9, frozen_mg_9 | count of peptides | as above, 9-50 residues |
| recovery_cf_5_pct, recovery_cf_9_pct | % | share of the cross-fitted MG peptides also accepted by cross-fitted FL, 5-50 and 9-50 residues |

## same_spectrum_benchmark.tsv
sha256: cb8248852f55dee7a8ac5fca054e88fd83cd85e2a84edfaed8a6543c2ec2e16b
rows: 1156 (one row per acquisition and length window: 296 PREDICT and 282 CapScan acquisitions x 2)
produced by: clinical_public_only_comparison

| column | unit | meaning |
|---|---|---|
| cohort | id | PREDICT or CAPSCAN |
| acquisition_id | id | as above |
| participant_id | id | PREDICT_P### (empty for the four acquisitions without a patient number) or CapScan participant |
| sampling_site | id | as above |
| min_length | residues | 5 (5-50 residues) or 9 (9-50 residues) |
| n_fl | count of peptides | peptides accepted with the public-sequence FastaLake database, same spectra |
| n_mg | count of peptides | peptides accepted with the specimen-matched metagenome database |
| overlap | count of peptides | peptides accepted by both |
| count_ratio | dimensionless | n_fl / n_mg |
| mg_recovery_fraction | fraction | overlap / n_mg |

## predict_acquisition_roster.tsv
sha256: 2622c52461d03bfa05055eae517c78a779c766e7460f64513ec28125901a78da
rows: 306 (one row per PREDICT acquisition)
produced by: clinical_public_only_comparison

| column | unit | meaning |
|---|---|---|
| acquisition_id | id | PREDICT_A### |
| participant_id | id | PREDICT_P### (empty if no patient number) |
| matched_metagenome | boolean | True for the 296 acquisitions with a specimen-matched metagenome database |

## spectrum_audit_per_acquisition.tsv
sha256: 590dac79f3ea0ac369c899db85998bbcd118b005e91d387370fafed179fe404e
rows: 296 (one row per matched PREDICT acquisition)
produced by: predict_fl_only_spectrum_audit

| column | unit | meaning |
|---|---|---|
| acquisition_id | id | PREDICT_A### |
| total | count of spectra | spectra accepted (peptide q ≤ 0.01) with FL but not with MG |
| A_absent | count of spectra | FL sequence absent from the MG database |
| D_displaced | count of spectra | FL sequence in the MG digest; a different MG peptide was the top target and failed the threshold |
| B_not_cleavable | count of spectra | FL sequence present in the MG database but outside the audited tryptic digest |
| C_top_fails_q | count of spectra | the same peptide was the top MG target but failed the threshold |
| E_no_target | count of spectra | no MG target candidate reported (up to five reported candidates per spectrum) |
| n_fl_spectra_q01, n_mg_spectra_q01, n_shared_spectra_q01 | count of spectra | spectra accepted with FL, with MG, and with both |

## clinical_paired_statistics.tsv
sha256: 9518b0117267eb3198b1aafe785dab47b53b020bd4deb0d9d2fb24f67ccc3bc7
rows: 2 (one row per endpoint; aggregate over 42 patients)
produced by: predict_clinical_refit

| column | unit | meaning |
|---|---|---|
| metric | id | Human (human-assigned share of LFQ intensity) or pg_parsimony (quantified protein groups) |
| n_patients | count of patients | patients sampled off and on antibiotics |
| mean_off, mean_on | % or count of protein groups | mean over patients of the per-state values |
| mean_on_minus_off | percentage points or count of protein groups | mean paired difference |
| CI95_low, CI95_high | as estimate | percentile bootstrap 95% CI (50,000 draws, patients resampled) |
| units | id | unit of the estimate |
| p_signflip_MC | dimensionless | two-sided Monte Carlo sign-flip P (200,000 flips) |
| q_two_clinical_endpoints | dimensionless | Benjamini-Hochberg q over the two endpoints |

## clinical_composition_means.tsv
sha256: c9c324f5be270ef15114c8c7ded33bba537edcd14e376cf1fb14512703d805a8
rows: 2 (Off, On; mean over 42 patients)
produced by: predict_clinical_refit (means written by the D12 figure builder)

| column | unit | meaning |
|---|---|---|
| condition | id | Off or On antibiotics (aggregate; renamed from antibiotic_state) |
| Human, Microbial_species, Microbial_above_species, Other_annotated_origins, Unresolved_origin | % of LFQ intensity | mean share of each origin class |

## heatmap_features.tsv
sha256: edabc3ec0e045136b12e12c96a541611311cf1c4aa2558c5779145ec3f58dc5e
rows: 20 (one row per feature of the fixed overview selection; 19 with n_positive_85 > 0 are drawn)
produced by: predict_eds_screen (selection: predict_eds_feature_selection)

| column | unit | meaning |
|---|---|---|
| block | id | Host proteome, Microbial KEGG pathway or Metabolomics |
| feature_id | id | UniProt entry, KEGG map or metabolomics feature id |
| label | id | name printed in Fig. 4h |
| n_positive_85 | count of patients | patients with a positive value of the feature |

## heatmap_z.tsv
sha256: af4e6a5413c3baaf88142c3e76b7188366c2499490f7fb777db02c3266a55415
rows: 1700 (one row per feature and patient; 20 x 85)
produced by: predict_eds_screen

| column | unit | meaning |
|---|---|---|
| feature_id | id | as in heatmap_features.tsv |
| display_position | id | column of Fig. 4h, 1-85, patients ordered by EDS |
| z | dimensionless | within-feature z score (empty where the feature is missing; drawn grey) |

## heatmap_tracks.tsv
sha256: 6aad0547cc064ff53582efbb8e68276a0bda252936c150c109a2f6c2a6a8e764
rows: 85 (one row per patient, in display order)
produced by: predict_eds_screen

| column | unit | meaning |
|---|---|---|
| display_position | id | 1-85 |
| eds_score | dimensionless (0-1) | metagenomic Enterotype Dysbiosis Score |
| human_fraction | fraction | human-assigned share of LFQ intensity |

## crossomics_native_correlations.tsv
sha256: dd1ddf3ff79eb8845de0162bcac3f9e5a722826226fba2ba9803728ebffaba00
rows: 32 (one row per cohort, analysis and layer pair)
produced by: predict_crossomics_refit

| column | unit | meaning |
|---|---|---|
| cohort | id | PREDICT or CAPSCAN |
| analysis | id | primary_prevalence10 (features in at least 10% of profiles) or sensitivity_all_features (all observed features; Fig. 4i) |
| layer_a, layer_b | id | DNA, RNA, Non-host protein, Host protein, Metabolites |
| rho | dimensionless | Spearman correlation of the two distance matrices |
| p | dimensionless | permutation P (1,999 permutations) |
| permutations | count | permutations |
| n_profiles, n_participants | count of patients / participants | profiles and participants |
| q | dimensionless | Benjamini-Hochberg q within cohort and analysis |

## crossomics_functional_correlations.tsv
sha256: 4c05b200a276d48e3730276f328749bde336191d373bd7a58c609b5aec1276b6
rows: 120 (one row per cohort, scope, level, layer pair and prevalence filter)
produced by: predict_crossomics_refit

| column | unit | meaning |
|---|---|---|
| cohort | id | PREDICT or CAPSCAN |
| scope | id | displayed77 (PREDICT, Fig. 4i) or CapScan scopes |
| level | id | OG (eggNOG orthologous group), KO (KEGG ortholog), Pathway (KEGG pathway) |
| layer_a, layer_b | id | metaG (DNA), metaT (RNA), microbial_metaP (protein excluding human-assigned), host_metaP (human protein), metabolites |
| prevalence | fraction | minimum prevalence of features (0 = all observed; Fig. 4i) |
| n_units, n_donors | count of patients | profiles and patients |
| rho, p, q | dimensionless | as above; q within cohort, scope and prevalence |

## crossomics_distances_predict77.tsv
sha256: 9267915b97de32a65a40f54388cca27225ccbf9dda4136172df4f7a5e39bedf1
rows: 2926 (one row per pair of the 77 PREDICT patients)
produced by: predict_crossomics_refit

| column | unit | meaning |
|---|---|---|
| patient_a, patient_b | id | patient position 1-77 in the cross-omics roster (a < b) |
| native__DNA, native__RNA, native__Non-host protein, native__Host protein, native__Metabolites | distance (dimensionless) | profile distance, native features, all observed features |
| OG__metaG ... Pathway__metabolites (15 columns) | distance (dimensionless) | profile distance at eggNOG OG, KEGG KO and KEGG pathway level, all observed features |

## capscan_dna_yield_by_acquisition.tsv
sha256: 351bdf31dc30f108e11796f15614cb05299ba8d37a8991d6ca4f10f7fa1a9c4d
rows: 282 (one row per matched CapScan acquisition)
produced by: capscan_metagenome_yield_join

| column | unit | meaning |
|---|---|---|
| sample | id | CapScan raw-file name |
| region | id | duodenum, jejunum, ileum, colon or stool |
| subject | id | participant 1-15 |
| biosample | id | NCBI BioSample of the metagenome |
| assembly_proteins | count of proteins | predicted proteins in the specimen-matched MG database |
| dna_ng | ng | recorded DNA amount (empty if not recorded) |
| volume_ul | µl | sample volume |
| FL_peptides, MG_peptides | count of peptides | same-spectrum peptides, 5-50 residues |
| gain_percent | % | 100 x (FL_peptides / MG_peptides - 1) |
| MG_identity_recovery | fraction | share of MG peptides also accepted by FL |

## EXTRACTION_RECEIPT.json

Input files (LAB_ROOT-relative) with sha256, and the sha256 of every output and of the two controlled files, written by `prepare_source_data.py`.
