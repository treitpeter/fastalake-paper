# Source data: Supplementary Data 1

Inputs of `tables/SData1/SData1.ipynb`. The `v28d_<sheet>.tsv` tables are the data sheets of the published
Supplementary Data 1 workbook (Source Data package copy, sha256 5ee76e28adb7...), one table per sheet, cell values
unchanged. PREDICT patient and tube codes are the package pseudonyms (PREDICT_Pnnn, PREDICT_Tnnnn).
Withheld (controlled, see `controlled/CONTROLLED_INPUTS.md`, file SData1_PREDICT_withheld_cells.tsv): the columns
august_visit_id, august_group and later_visit_number (PREDICT metadata, mapping changes, repeated tubes),
recorded_visit_number (PREDICT stack eligibility), and the PREDICT_ANTIBIOTIC42 rows of PREDICT membership (306 rows)
and PREDICT patients (93 rows). The Cohorts and Sheet directory sheets are computed by the notebook and are not here.
The v28f CapScan roster also reads `source_data/Fig5/design_347.tsv` (unit Fig5).

## v28d_analyses.tsv
sha256: c7098876aff21b2837de7e619dadd05703a5c23cd5ec84272f5cca34f228156c
rows: 30 (one row per analysis)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| analysis_id | id | analysis or design subset; one row of the Analyses sheet |
| cohort | id | dataset label |
| purpose | text | question the analysis answers |
| analysis_unit | text | statistical or aggregation unit of the analysis |
| n_MS_acquisitions | count of acquisitions | acquisitions in the analysis roster |
| n_recorded_proteomic_specimens | count of specimens | recorded proteomic specimens, or a statement why none are recorded |
| n_participants | count of participants or patients | participants or patients, or a statement on donor identity |
| n_DNA_runs | count of DNA sequencing runs | DNA runs, or 'Not applicable' |
| n_DNA_specimens | count of specimens | DNA specimens, or a pointer |
| n_participant_region_cells | count of participant-region cells | participant-region cells, or 'Not applicable' |
| roster_sheet | id | workbook sheet that lists the members |
| roster_filter | text | filter to apply to roster_sheet |
| exclusion_sheet | id | workbook sheet that lists nonmembers |
| exclusion_filter | text | filter to apply to exclusion_sheet |
| inclusion_rule | text | rule that defines membership |
| aggregation | text | how acquisitions are aggregated |
| limitations | text | what the roster does not establish |
| provenance | id | producing source labels or table |

## v28d_campi_acquisitions.tsv
sha256: 12a68ec2fdd8ff26af71833e5648ddc98849853d29995c4b7d9ac63b959a3985
rows: 15 (one row per acquisition)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| cohort | id | dataset label |
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| deposit_raw_filename | id | raw-file name in PRIDE PXD023217 |
| publication_workflow_id | id | CAMPI workflow label S01-S12 |
| publication_file_set | id | file names in the CAMPI publication |
| fraction | count | fraction number within S11 (blank otherwise) |
| fractions_in_workflow | count of fractions | fractions in the workflow |
| material | text | sample material |
| source_specimen_id | id | specimen identifier in the source metadata (blank if none) |
| source_participant_id | id | participant identifier in the source metadata (blank if none) |
| participant_status | text | why participant identity does not apply |
| instrument_from_acquisition | id | instrument model read from the acquisition |
| submitting_team | id | laboratory that produced the acquisition |
| source_extraction_design | text | extraction design in the source |
| pride_accession | id | PRIDE accession |
| identity_note | text | note on specimen identity |
| source_ids | id | short labels of the producing source records, ';'-separated |

## v28d_campi_search_arms.tsv
sha256: c3e23a483d2a7016839cdf1030ccede882882327a4199415f49a7015e2f30a98
rows: 120 (one row per acquisition and search arm)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| analysis_id | id | analysis or design subset; one row of the Analyses sheet |
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| arm | id | database arm (transfer/pool label, or CAMPI search arm) |
| engine | id | search engine |
| source_config | id | search result or configuration file name |
| source_id | id | short label of the producing source record |

## v28d_capscan_acquisitions.tsv
sha256: 045c279f18f681783b1fffc545007376be7b7736dfa0e62dbc0655430e259d13
rows: 331 (one row per acquisition)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| specimen_id_biology | id | specimen in the biology analyses ('main:<aliquot main id>') |
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| region | id | sampling group: duodenum, jejunum, ileum, colon (CapScan regions 1-4) or stool |
| aliquot | id | aliquot type (pellet for intestinal fluid, main for stool) |
| aliquot_main_id | id | device aliquot (Main) number in the source workbook |
| specimen_id_crossomics | id | specimen in the cross-omics analyses |
| link_basis | text | how the acquisition was linked to specimen, participant and region |

## v28d_capscan_dna_eligibility.tsv
sha256: 9fec76d3b867647e4921a28ac479b1974d7e079cca4d4007ee42330211aa8094
rows: 210 (one row per participant, region and representation)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| region | id | sampling group: duodenum, jejunum, ileum, colon (CapScan regions 1-4) or stool |
| representation | id | functional representation (OG, KO, Pathway) |
| eligible_DNA_runs | count of DNA runs | eligible DNA runs in the cell |
| eligible_DNA_specimens | count of specimens | eligible DNA specimens in the cell |
| DNA_region_available | boolean | DNA available for the cell |
| in_complete12 | boolean | participant in the 12-participant complete DNA set |

## v28d_capscan_dna_materials.tsv
sha256: 1bde51285dceb71e35d9ce1a8bac144de8e3ce1204d7704224e16fec0cae67c9
rows: 6 (one row per DNA source region or material)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| source_region_or_material | id | DNA source region or material |
| DNA_runs | count of DNA runs | DNA sequencing runs |
| DNA_specimens | count of specimens | DNA specimens |
| participants | count of participants | participants |

## v28d_capscan_dna_runs.tsv
sha256: 3a8978311c2fc38e03ce7cfeb6998c3c6f868df2bcb44c20e2d4ef850f47593f
rows: 301 (one row per DNA sequencing run)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| run_accession | id | public sequencing run accession (DNA rows) |
| specimen_id_DNA | id | DNA specimen |
| participant_id_DNA | id | participant of the DNA specimen |
| region_DNA | id | sampling group or material of the DNA specimen |

## v28d_capscan_exclusions.tsv
sha256: 25387080b9c0aa36007493fd3bec339f4d6b536034453337a1f4b362e8ad0f0a
rows: 140 (one row per excluded entity per analysis)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| analysis_id | id | analysis or design subset; one row of the Analyses sheet |
| entity_type | id | acquisition, specimen, participant, participant_region, DNA_run or DNA_specimen |
| entity_id | id | identifier of the entity |
| reason | text | recorded reason for inclusion or exclusion |
| source_file | id | file name of the producing table |
| source_sha256 | id | sha256 of the producing table |
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| region | id | sampling group: duodenum, jejunum, ileum, colon (CapScan regions 1-4) or stool |
| specimen_id | id | specimen (PREDICT tube pseudonym; CapScan 'main:<aliquot>') |
| missing_region | id | sampling group missing for the participant |
| disposition_type | id | class of nonmembership |

## v28d_capscan_h_acquisitions.tsv
sha256: 636f54a49daeb9a2722d3f58e8d94c8b4f6c220faae3bf6bdc26555e530557e9
rows: 331 (one row per acquisition)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| histidine_specimen_id | id | specimen in the histidine analysis |
| specimen_id_crossomics | id | specimen in the cross-omics analyses |
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| region | id | sampling group: duodenum, jejunum, ileum, colon (CapScan regions 1-4) or stool |
| aliquot | id | aliquot type (pellet for intestinal fluid, main for stool) |
| in_quantitative330_roster | boolean | in the 330 nonempty quantitative roster |
| histidine_scope | text | scope of the histidine roster |
| source_file | id | file name of the producing table |

## v28d_capscan_h_finite_n.tsv
sha256: a4bbf72ee5aece2273afbe7fcaccd956c060215801e900c87f2142d6acdd8273
rows: 7 (one row per histidine endpoint family)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| family | id | statistical family |
| endpoint | id | histidine endpoint (route|terminus|source class|residue) |
| saved_tests | count of tests | tests saved in the family |
| minimum_saved_n | count of participants | smallest saved finite n |
| maximum_saved_n | count of participants | largest saved finite n |
| saved_estimable_tests | count of tests | tests with an estimate |
| scope | text | what the statement covers and does not cover |
| source_file | id | file name of the producing table |

## v28d_capscan_h_participants.tsv
sha256: 2ce20879be7ebe8233fa10746bac9898015357d41f23b2d4184b5f19b0257913
rows: 467 (one row per histidine endpoint and region)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| endpoint_kind | id | kind of histidine endpoint |
| route | id | selection route |
| side | id | peptide terminus |
| origin | id | source class |
| region | id | sampling group: duodenum, jejunum, ileum, colon (CapScan regions 1-4) or stool |
| finite_participants | count of participants | participants with a finite value |
| participant_ids | id | those participants, ';'-separated |
| source_rows | count of rows | source rows behind the set |
| source_file | id | file name of the producing table |

## v28d_capscan_id_crosswalk.tsv
sha256: 82eac43283744af32bf88d39e51db301535fa3dcfbc19b8cc11af9d67fa028b9
rows: 331 (one row per acquisition)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| specimen_id_biology | id | specimen in the biology analyses ('main:<aliquot main id>') |
| specimen_id_crossomics | id | specimen in the cross-omics analyses |
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| region | id | sampling group: duodenum, jejunum, ileum, colon (CapScan regions 1-4) or stool |
| link_basis | text | how the acquisition was linked to specimen, participant and region |

## v28d_capscan_membership.tsv
sha256: 5e66c8368394013f1c40e5b79079ca3c4cf75381ef5cc89e3804466675a55daf
rows: 7599 (one row per included entity per analysis)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| analysis_id | id | analysis or design subset; one row of the Analyses sheet |
| entity_type | id | acquisition, specimen, participant, participant_region, DNA_run or DNA_specimen |
| entity_id | id | identifier of the entity |
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| specimen_id_biology | id | specimen in the biology analyses ('main:<aliquot main id>') |
| specimen_id_crossomics | id | specimen in the cross-omics analyses |
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| region | id | sampling group: duodenum, jejunum, ileum, colon (CapScan regions 1-4) or stool |
| run_accession | id | public sequencing run accession (DNA rows) |
| source_file | id | file name of the producing table |
| source_sha256 | id | sha256 of the producing table |
| inclusion_rule | text | rule that defines membership |

## v28d_capscan_mg_pairs.tsv
sha256: d6bf44441bc1c40732b2a4847fb24183cf84745de47ca05ad4a8d47071399f51
rows: 295 (one row per matched acquisition)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| specimen_id_A16 | id | specimen key used by the A16 benchmark (participant|capsule label) |
| specimen_id_biology | id | specimen in the biology analyses ('main:<aliquot main id>') |
| specimen_id_crossomics | id | specimen in the cross-omics analyses |
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| region | id | sampling group: duodenum, jejunum, ileum, colon (CapScan regions 1-4) or stool |
| n_test | count | value of the A16 task field n_test (test arm, public FastaLake) |
| n_reference | count | value of the A16 task field n_reference (specimen-matched metagenome arm) |

## v28d_capscan_ms_nonmembership.tsv
sha256: 71578d7ac53adbbd8fceb52002812383dbf850c5267874a721d4dacec5bbbfcd
rows: 188 (one row per acquisition not in an analysis)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| analysis_id | id | analysis or design subset; one row of the Analyses sheet |
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| specimen_id_biology | id | specimen in the biology analyses ('main:<aliquot main id>') |
| specimen_id_crossomics | id | specimen in the cross-omics analyses |
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| region | id | sampling group: duodenum, jejunum, ileum, colon (CapScan regions 1-4) or stool |
| disposition_type | id | class of nonmembership |
| reason | text | recorded reason for inclusion or exclusion |

## v28d_capscan_regional_pairs.tsv
sha256: fdfcbcd5cf09d283c2187129c04f2512313da69bde29cb262c5100ae12c04668
rows: 56 (one row per participant and capsule region)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| region | id | sampling group: duodenum, jejunum, ileum, colon (CapScan regions 1-4) or stool |
| paired_stool_region | id | stool group paired with the capsule region |
| source_file | id | file name of the producing table |

## v28d_hs_acquisitions.tsv
sha256: e8231bd74d416a14b8994bf0d97fb651794138b954fa75ebfe03461b8ab3b395
rows: 58 (one row per acquisition)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| cohort | id | dataset label |
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| mzml_file | id | converted spectrum file name |
| source_specimen_id | id | specimen identifier in the source metadata (blank if none) |
| source_participant_id | id | participant identifier in the source metadata (blank if none) |
| analysis_participant_label | id | donor label used in the analysis (HS_shared_donor or Ind1-Ind10) |
| participant_label_basis | text | basis of the donor label |
| design | id | healthy-stool design: TechnicalReplicates, WorkflowReplicates, OnePiece, Intra, Interindividual |
| biological_question | text | question the design answers |
| sample_relationship | text | documented relationship to other designs |
| identity_limit | text | what is not known about the sample identity |
| metadata_source_id | id | source metadata record |
| identity_source_ids | id | records that establish the identity statements |

## v28d_hs_campi_membership.tsv
sha256: 032e99fdfb1311507b5bb944b591355308fe03f3cbd14671c77fc8b1b366eb77
rows: 591 (one row per acquisition per analysis)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| analysis_id | id | analysis or design subset; one row of the Analyses sheet |
| cohort | id | dataset label |
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| included | yes/no | acquisition is in the analysis roster |
| reason | text | recorded reason for inclusion or exclusion |
| source_participant_id | id | participant identifier in the source metadata (blank if none) |
| analysis_participant_label | id | donor label used in the analysis (HS_shared_donor or Ind1-Ind10) |
| design | id | healthy-stool design: TechnicalReplicates, WorkflowReplicates, OnePiece, Intra, Interindividual |
| source_ids | id | short labels of the producing source records, ';'-separated |

## v28d_hs_excluded.tsv
sha256: 3d722509b96f579543f0a9bc03e165d6e835b85e4636c604b9395372cfd407e0
rows: 39 (one row per excluded acquisition)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| cohort | id | dataset label |
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| included_in_primary | boolean | in the 58-acquisition primary roster |
| reason | text | recorded reason for inclusion or exclusion |
| mzml_file | id | converted spectrum file name |
| source_id | id | short label of the producing source record |

## v28d_hs_pools.tsv
sha256: 10748ef2191b98b5e720f80c9f6c6dfc119301b297c9df97494eabc890bb0522
rows: 288 (one row per target acquisition, arm and database-source acquisition)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| target_acquisition | id | acquisition searched |
| arm | id | database arm (transfer/pool label, or CAMPI search arm) |
| database_source_acquisition | id | acquisition whose de novo peptides selected the database |
| source_id | id | short label of the producing source record |

## v28d_hs_relationships.tsv
sha256: d4d5ccd72d5b0cc88ce96830126fd82bb2912b8140e183aeb96e7246302e57ff
rows: 3 (one row per documented relationship)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| subject | id | participant (1-15; blank if not linked) |
| relationship | text | documented relationship |
| object | id | design or label related to |
| source | id | record that documents the relationship |
| scope | text | what the statement covers and does not cover |

## v28d_hs_transfer.tsv
sha256: fe087c3fbcd285e37b88c85e3fc5127a2b21ba98aebdfe9d1a75685e8be2fa02
rows: 174 (one row per target acquisition and arm)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| target_acquisition | id | acquisition searched |
| arm | id | database arm (transfer/pool label, or CAMPI search arm) |
| database_source_acquisition | id | acquisition whose de novo peptides selected the database |
| target_design | id | design of the target acquisition |
| source_design | id | design of the database-source acquisition |
| target_participant_label | id | donor label of the target |
| database_source_participant_label | id | donor label of the database source |
| biological_independence | text | whether target and source are from the same donor |
| input_fasta | id | database file name |
| source_id | id | short label of the producing source record |

## v28d_predict_85_to_77.tsv
sha256: 4bf7b6279d6eb5242e54c78a14ef762aa2d21f65d37b483ae64aa906afd5985e
rows: 8 (one row per patient)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| n_later_samples | count of samples | samples of the patient in the later metadata |
| n_EDS_eligible_samples | count of samples | samples eligible for the shared-profile (EDS) analysis |
| n_native_complete_samples | count of samples | samples with the complete assay stack |
| reason | text | recorded reason for inclusion or exclusion |
| absent_assays_across_EDS_eligible_samples | text | missing assay profiles, as missing/eligible samples |

## v28d_predict_h_patients.tsv
sha256: e4cc6af0db07401cbcc7df588e8d742b78357059300cc260b55eb8233a76290c
rows: 1116 (one row per patient per histidine endpoint)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| analysis_id | id | analysis or design subset; one row of the Analyses sheet |
| endpoint | id | histidine endpoint (route|terminus|source class|residue) |
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| status | id | included or excluded (PREDICT), or the census status (discrepancy table) |
| reason | text | recorded reason for inclusion or exclusion |

## v28d_predict_mapping_changes.tsv
sha256: e20452df063f030fa51ae4de6fa9c790df586a939ff1ca582800de7d7f7f5386
rows: 4 (one row per acquisition)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| recorded_sample_number | id | sample number recorded in the acquisition metadata |
| august_specimen_tube_id | id | specimen tube pseudonym (PREDICT_Tnnnn) in the August metadata |
| august_participant_id | id | patient pseudonym (PREDICT_Pnnn) in the August metadata |
| august_join_status | id | join result in the August metadata |
| august_aliquot | id | aliquot type in the August metadata |
| later_analytic_sample_id | id | analytic sample id in the later (268-sample) metadata |
| later_specimen_tube_id | id | specimen tube pseudonym in the later metadata |
| later_participant_id | id | patient pseudonym in the later metadata |
| later_join_basis | text | how the later metadata were joined |
| later_mapping_status | id | join result in the later metadata |
| identity_status | text | agreement of the two metadata versions |

## v28d_predict_membership.tsv
sha256: dcb0b12ba1f9fabfeae9528301a734d577ecaa291a8436b53eb1f2e8c1f958f5
rows: 2142 (one row per acquisition per analysis (antibiotic-pair rows withheld))
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| analysis_id | id | analysis or design subset; one row of the Analyses sheet |
| status | id | included or excluded (PREDICT), or the census status (discrepancy table) |
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| analytic_sample_id | id | analytic sample id |
| recorded_sample_number | id | sample number recorded in the acquisition metadata |
| specimen_id | id | specimen (PREDICT tube pseudonym; CapScan 'main:<aliquot>') |
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| identity_metadata_version | id | metadata version that supplies the identity (August or later) |
| reason | text | recorded reason for inclusion or exclusion |
| blank_identity_reason | text | why the identity is blank |
| source_ids | id | short labels of the producing source records, ';'-separated |
| reason_type | id | class of reason |

## v28d_predict_metadata.tsv
sha256: ed617c34c305a103119b5a78f0bec1520696dd2529621ca51ee09a489068a8a3
rows: 306 (one row per acquisition)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| recorded_sample_number | id | sample number recorded in the acquisition metadata |
| august_specimen_tube_id | id | specimen tube pseudonym (PREDICT_Tnnnn) in the August metadata |
| august_participant_id | id | patient pseudonym (PREDICT_Pnnn) in the August metadata |
| august_join_status | id | join result in the August metadata |
| august_aliquot | id | aliquot type in the August metadata |
| later_analytic_sample_id | id | analytic sample id in the later (268-sample) metadata |
| later_specimen_tube_id | id | specimen tube pseudonym in the later metadata |
| later_participant_id | id | patient pseudonym in the later metadata |
| later_join_basis | text | how the later metadata were joined |
| later_mapping_status | id | join result in the later metadata |
| identity_status | text | agreement of the two metadata versions |

## v28d_predict_mg_exclusions.tsv
sha256: 64cf1502359c898a48146c578b85720c2f707fed1f5f6791ce7bdb2d87a73f5a
rows: 10 (one row per acquisition)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| recorded_sample_number | id | sample number recorded in the acquisition metadata |
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| MG_set_in_frozen_manifest | boolean | matched metagenome peptide set present in the frozen comparison manifest |
| MG_config_in_frozen_manifest | boolean | matched metagenome configuration present in the frozen manifest |
| available_arms | id | search arms available for the acquisition |
| reason | text | recorded reason for inclusion or exclusion |

## v28d_predict_patients.tsv
sha256: b2e4d1d1b2ad6c27485143e90aacd4aa3338972600a5b3d08a1fcad4f9d5945f
rows: 372 (one row per patient per analysis (antibiotic-pair rows withheld))
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| analysis_id | id | analysis or design subset; one row of the Analyses sheet |
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| status | id | included or excluded (PREDICT), or the census status (discrepancy table) |
| reason | text | recorded reason for inclusion or exclusion |

## v28d_predict_repeated_tubes.tsv
sha256: ca496194a0181230476808aaefb8f0bee28bc7880e007064d20dd8096a3afa29
rows: 4 (one row per acquisition)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| recorded_sample_number | id | sample number recorded in the acquisition metadata |
| august_specimen_tube_id | id | specimen tube pseudonym (PREDICT_Tnnnn) in the August metadata |
| august_participant_id | id | patient pseudonym (PREDICT_Pnnn) in the August metadata |
| august_join_status | id | join result in the August metadata |
| august_aliquot | id | aliquot type in the August metadata |
| later_analytic_sample_id | id | analytic sample id in the later (268-sample) metadata |
| later_specimen_tube_id | id | specimen tube pseudonym in the later metadata |
| later_participant_id | id | patient pseudonym in the later metadata |
| later_join_basis | text | how the later metadata were joined |
| later_mapping_status | id | join result in the later metadata |
| identity_status | text | agreement of the two metadata versions |

## v28d_predict_stack_eligibility.tsv
sha256: 8bcec9d7e42903251849ea6e6a5feb724963e8791213b09705829b12498416d0
rows: 268 (one row per analytic sample)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| analytic_sample_id | id | analytic sample id |
| participant_id | id | participant (CapScan: study number 1-15; PREDICT: pseudonym PREDICT_Pnnn) |
| specimen_id | id | specimen (PREDICT tube pseudonym; CapScan 'main:<aliquot>') |
| metadata_prediction262 | boolean | sample in the 262-sample prediction metadata |
| missing_required_prediction_fields | text | names of required metadata fields that are missing (names only, no values) |
| has_DNA_OG_profile | boolean | sample has a DNA OG profile |
| has_RNA_OG_profile | boolean | sample has a RNA OG profile |
| has_DNA_KO_profile | boolean | sample has a DNA KO profile |
| has_RNA_KO_profile | boolean | sample has a RNA KO profile |
| has_DNA_Pathway_profile | boolean | sample has a DNA Pathway profile |
| has_RNA_Pathway_profile | boolean | sample has a RNA Pathway profile |
| has_MicrobialProtein_OG_profile | boolean | sample has a MicrobialProtein OG profile |
| has_MicrobialProtein_KO_profile | boolean | sample has a MicrobialProtein KO profile |
| has_MicrobialProtein_Pathway_profile | boolean | sample has a MicrobialProtein Pathway profile |
| has_HostProtein_profile | boolean | sample has a HostProtein profile |
| has_Metabolites_profile | boolean | sample has a Metabolites profile |
| EDS_shared_profile_eligible | boolean | has every profile needed for the shared-profile analysis |
| EDS_primary_selected | boolean | selected as the patient's primary sample |
| EDS_not_selected_reason | text | why not selected |
| native_complete_profile_eligible | boolean | has the complete native profile stack |
| native_primary_selected | boolean | selected as the patient's complete-stack sample |
| native_not_selected_reason | text | why not selected |

## v28d_read_me.tsv
sha256: 15be0c76da66ab4859bbf16decce3da0907e0bbc0e4c9fcd1a807ee3db197119
rows: 11 (one row per workbook topic)
produced by: sample_rosters_workbook; pseudonymised by source_data_package
origin: sheet of Source Data package Supplementary_Data_1/Sample_rosters.xlsx; values agree with the package (checked by the notebook)

| column | unit | meaning |
|---|---|---|
| topic | text | workbook topic |
| explanation | text | how to read the workbook for that topic |

## CAPSCAN_raw_header_audit.tsv
sha256: 8f39a925b8cb774a0fd1925ea372214a984e8e357bc8e00a6438c5c5c85076c5
rows: 368 (one row per CapScan RAW file (368 audited))
produced by: capscan_raw_header_audit
origin: RAW header audit all_RAW_header_names.tsv; six columns kept, local file paths and checksums dropped

| column | unit | meaning |
|---|---|---|
| sample | id | raw-file name of the CapScan acquisition |
| current_331 | boolean | acquisition in the 331-acquisition CapScan roster |
| subject | id | participant (1-15; blank if not linked) |
| header_original_names | id | original file name stored in the instrument record of the RAW file |
| header_specimen_number_changed | boolean | specimen number in the instrument record differs from the file name |
| header_status | id | result of reading the instrument record |

## CAPSCAN_eligibility_audit_347.tsv
sha256: 3c58ea360554a7631bdef00c21fad5b394ebd0b4fbfab72b6ff537e585ceb36b
rows: 347 (one row per acquisition (CapScan 347 roster))
produced by: capscan_crossomics
origin: A82 output eligibility_audit.tsv (347-roster cross-omics run), unchanged

| column | unit | meaning |
|---|---|---|
| acquisition | id | raw-file name of the LC-MS/MS acquisition |
| previous_positive | boolean | positive non-host OG profile in the 331-roster run |
| current_positive | boolean | positive total quantified profile in the 347-roster run (A82) |

## CAPSCAN_participant15_join.tsv
sha256: 688dea12752c41fb95015bb6350156e5d1d9522511cd201a09546ea3718e559c
rows: 16 (one row per acquisition (participant-15 intestinal fluid))
produced by: capscan_participant_roster
origin: A80 output participant15_added.tsv (participant-15 join fix), unchanged

| column | unit | meaning |
|---|---|---|
| sample_id | id | raw-file name of the LC-MS/MS acquisition |
| region | id | sampling group: duodenum, jejunum, ileum, colon (CapScan regions 1-4) or stool |
| region_ordinal | dimensionless | region order 1-5 (duodenum to stool) |
| subject_id | id | participant (1-15) |
| how | id | join route used |
| overview_rows | count of rows | matching rows in the sample overview |
| workbook_Main | id | device aliquot (Main) in the source workbook |
| workbook_Super | id | supernatant aliquot in the source workbook |
| overview_Pellet | id | Pellet field of the sample overview row |
| workbook_Subject | id | participant in the source workbook |
| workbook_Set | id | capsule set in the source workbook |
| workbook_Type | id | capsule type (Type1-4) in the source workbook |
| overview_Type | id | capsule type in the sample overview |
| subject_and_type_agree | boolean | workbook and overview agree on participant and capsule type |

## CapScan_specimen_number_discrepancies.tsv
sha256: 2f4a14503dd71c1d5ce7440949f2406eefdc36efb1ce0910841de22ba83087e6
rows: 16 (one row per acquisition with discordant specimen numbers)
produced by: source_data_package
origin: Source Data package Supplementary_Data_1/CapScan_specimen_number_discrepancies.tsv, unchanged; recomputed by the notebook from CAPSCAN_raw_header_audit.tsv

| column | unit | meaning |
|---|---|---|
| acquisition | id | raw-file name of the LC-MS/MS acquisition |
| in_analysis_roster | yes/no | acquisition is in the CapScan analysis roster (v28d: 331) |
| roster_participant | id | participant in the roster |
| specimen_number_in_file_name | id | specimen number in the raw-file name |
| specimen_number_in_instrument_record | id | specimen number in the instrument record inside the raw file |
| status | id | included or excluded (PREDICT), or the census status (discrepancy table) |

## v30 inputs

Copied once by `tables/SData1/prepare_v30_inputs.py` (Slurm job 6075320) for the v30 workbook, in which every CapScan analysis is on the 347-acquisition roster.

## v30_capscan_crossomics_specimens.tsv
sha256: d6b45a61b385f19ae6380b8bb956bd9c7c33519fe4553de53d8c9456d61c6553
rows: 303 (one row per specimen)
produced by: capscan_crossomics (347 cross-omics run, Slurm 6063849)
origin: capscan_specimen_roster.tsv of the run, columns specimen, participant and region only (device and covariate columns dropped)

| column | unit | meaning |
|---|---|---|
| specimen_id_crossomics | id | CapScan Main specimen id |
| participant_id | id | CapScan participant 1-15 |
| region | id | duodenum, jejunum, ileum, colon or stool |

## v30_capscan_crossomics_participants.tsv
sha256: 5ace9f60f1f56061577bf2e4a1e83d73f50b7f5e86da906a2d4406851f6c8239
rows: 15 (one row per participant)
produced by: capscan_crossomics (Slurm 6063849)
origin: COMPLETE.json of the run (protein_metabolite_participants, dna_participants)

| column | unit | meaning |
|---|---|---|
| participant_id | id | CapScan participant 1-15 |
| protein_metabolite_complete | bool | valid protein and metabolite specimens in all five sampling groups |
| dna_complete | bool | also DNA profiles in all five sampling groups |

## v30_capscan_metabolite_specimens.tsv
sha256: 81c18e05ab7ce277bfd8b8ebb088b7636c8ed1c27647faae80ce8761832877b6
rows: 318 (one row per specimen)
produced by: capscan_crossomics (extended metabolite table of the 347 run)
origin: index of CAPSCAN_metabolite_intensity.tsv; identifiers only, no values

| column | unit | meaning |
|---|---|---|
| specimen_id_crossomics | id | specimen with a row in the metabolite table |

## v30_capscan_h_design.tsv
sha256: 1d24a77c8d50694cadda71434e9bc1b6a430e42506b1f31b2d023411b2019740
rows: 347 (one row per acquisition)
produced by: terminal_histidine_inference (fixed-347 A28 run, Slurm 6064031)
origin: CAPSCAN tasks of inputs/manifest_revision4.json

| column | unit | meaning |
|---|---|---|
| acquisition_id | id | raw-file name of the LC-MS/MS acquisition |
| participant_id | id | CapScan participant |
| region | id | sampling region |
| histidine_specimen_id | id | specimen id used by the histidine run (main:<Main id>) |

## v30_capscan_h_finite_n.tsv
sha256: b9d374e66792e2b74ba165e784fe12980835f5cb142f932e2f8b972f8229029e
rows: 7 (one row per testing family and endpoint)
produced by: terminal_histidine_inference (Slurm 6064031)
origin: complete_statistical_families.tsv of the fixed-347 run, summarised by prepare_v30_inputs.py:summarise_histidine (code of build_histidine_scope.py)

| column | unit | meaning |
|---|---|---|
| family | id | statistical family |
| endpoint | id | endpoint |
| saved_tests | count of tests | tests saved in the family |
| minimum_saved_n | count of participants | smallest saved n |
| maximum_saved_n | count of participants | largest saved n |
| saved_estimable_tests | count of tests | tests marked estimable |
| scope | text | what the summary establishes |
| source_file | id | producing table |

## v30_capscan_h_participants.tsv
sha256: 42de614359e490e70818471159831ffe15409ba97ea2bd6ef17ed3512413d32f
rows: 467 (one row per endpoint cell)
produced by: terminal_histidine_inference (Slurm 6064031)
origin: participant_region_enrichment.tsv (residue H) and participant_reciprocal_controls.tsv of the fixed-347 run, summarised by prepare_v30_inputs.py:summarise_histidine

| column | unit | meaning |
|---|---|---|
| endpoint_kind | id | H_source_difference or complementary_product_difference |
| route | id | identification route |
| side | id | boundary side (blank for reciprocal controls) |
| origin | id | origin class |
| region | id | sampling region |
| finite_participants | count of participants | participants with a finite difference |
| participant_ids | id | those participants, semicolon-separated |
| source_rows | count of rows | rows of the producing table |
| source_file | id | producing table |

## v30_inputs_receipt.tsv
sha256: e0352be2814111f8da3a084412e253d15c55f80f019a5c21425d7f1d37a402b9
rows: 8 (one row per file read)
produced by: tables/SData1/prepare_v30_inputs.py

| column | unit | meaning |
|---|---|---|
| file | path | file read, as LAB_ROOT/<relative path> |
| sha256 | hash | sha256 of the file |
