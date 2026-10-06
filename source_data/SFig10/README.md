# Source data: Supplementary Figure 10

Tables read by `figures/SFig10/SFig10.ipynb`. They were derived from the original
producer inputs by `figures/SFig10/prepare_source_data.py` (origins and hashes in
`figures/SFig10/PROVENANCE.md`). Patients are keyed by their display position (1-85, EDS
order) and the pseudonym used in the Source Data package. Antibiotic state, clinical
group, raw patient identifiers and sample identifiers are not here; the notebook reads
antibiotic state and clinical group from the controlled input
`patient_order_and_covariates.tsv` (see `controlled/CONTROLLED_INPUTS.md`).

Read floats with `pandas.read_csv(..., sep="\t", float_precision="round_trip")`; the
values then equal the producer's arrays bit for bit.

## patients.tsv
sha256: 288ecd87333c04f381ef23ca51eb4c85581df2fa7bff34865c8c25277e24e9bf
rows: 85 (one row per patient, earliest matched proteome-metabolome observation, in EDS order)
produced by: predict_eds_focused_pairs

| column | unit | meaning |
|---|---|---|
| patient_position | id | display position in Supplementary Fig. 10 and Fig. 4h (1 = lowest EDS) |
| patient_pseudonym | id | pseudonymised patient identifier (PREDICT_Pnnn), as in the Source Data package |
| dysbiosis_score | dimensionless | metagenomic Enterotype Dysbiosis Score (0-1) |
| Human | fraction | human-assigned share of total LFQ intensity of the acquisition |
| pg_parsimony | count of protein groups | quantified protein groups of the acquisition (protein-group depth) |
| ko_total_signal | LFQ intensity | sum of the KO-allocated signal over all KOs (denominator of the annotated-relative normalization) |

## focused_values.tsv
sha256: a4607f929b739084a081f23f1afd2bc16dbcf90a7a4cc4dffd5dcce11cb6e53c
rows: 510 (6 features x 85 patients)
produced by: predict_eds_focused_pairs

| column | unit | meaning |
|---|---|---|
| feature | id | metabolite feature id (HILIC_/C18_) or KEGG orthology id (K) |
| patient_position | id | as in patients.tsv |
| raw | metabolite: intensity (a.u.); KO: LFQ intensity | measured value; 0 or empty = not detected |
| z | dimensionless | (value - mean) / SD over positive finite values of the feature (population SD); empty where not detected |

## focused_features.tsv
sha256: 7945c40a934016cbfab1952c42b02cdf2854fb923b3815b166b308a81e352aae
rows: 6 (one row per displayed feature of the figure (v30 single panel, formerly panel a), in plotting order)
produced by: predict_eds_focused_pairs

| column | unit | meaning |
|---|---|---|
| feature | id | feature id |
| label | id | row label in the figure |
| layer | id | Metabolite or KO |
| n_positive | count of patients | patients with a positive finite value |
| EDS_rho | dimensionless | Spearman rho with EDS over positive values |
| EDS_full_family_q | dimensionless | Benjamini-Hochberg q in the complete native unadjusted EDS family of the layer |

## legacy_85_values.tsv
sha256: 3f26766a6c50dae457c5cf84e260960bbfa89002314ba9bc68cfe58efc14dcea
rows: 1700 (20 original overview features x 85 patients)
produced by: predict_eds_screen (feature selection: predict_eds_feature_selection)

| column | unit | meaning |
|---|---|---|
| feature | id | UniProt entry, KEGG map id or metabolite feature id |
| patient_position | id | as in patients.tsv |
| raw | human protein and pathway: LFQ intensity; metabolite: intensity (a.u.) | measured value |
| z | dimensionless | within-feature z score over positive finite values; empty where not detected |

## legacy_features.tsv
sha256: b883e1659743f5437c3d1d3981abd07a36d1ae41ee5d87b6321596fe90968e51
rows: 20 (one row per feature of the original EDS overview selection)
produced by: predict_eds_screen (feature selection: predict_eds_feature_selection)

| column | unit | meaning |
|---|---|---|
| block | id | Host proteome, Microbial KEGG pathway or Metabolomics |
| feature_id | id | feature id |
| label | id | row label in the former panel b (v28f; removed in v30, identical to Fig. 4h) |
| n_positive_85 | count of patients | patients (of 85) with a positive value; 0 for HV366, which is not drawn |
| rho_EDS_85 | dimensionless | Spearman rho with EDS over positive values in the 85 patients |

## legacy_259_values.tsv
sha256: 0e3dc2a157f9fdd6c01cbd992f65b454e15c0e28c86e988d9b337bc6b613893b
rows: 5180 (20 features x 259 matched observations, including repeated visits)
produced by: predict_eds_screen

| column | unit | meaning |
|---|---|---|
| feature | id | as in legacy_features.tsv |
| observation_index | id | 1-259, order of the original supplementary sample view; observation identifiers withheld |
| raw | as in legacy_85_values.tsv | measured value |

## augmented_pair_screen.tsv
sha256: d7e59ab759d62aa5b2068ac5fbe0b774ece288ad3894bd633cfd219cac2a3967
rows: 11496 (one row per linked metabolite-KO or metabolite-pathway pair, normalization and model)
produced by: predict_eds_screen

| column | unit | meaning |
|---|---|---|
| layer | id | KO or Pathway |
| normalization | id | native or annotated_relative |
| model | id | unadjusted, host_depth_antibiotic_clinical, or host_depth_antibiotic_clinical_EDS |
| metabolite | id | metabolite feature id |
| name | id | metabolite annotation |
| compound | id | KEGG compound |
| protein_feature | id | KO or KEGG map id |
| n | count of patients | patients with positive finite values (and complete covariates) |
| rho | dimensionless | Spearman or partial rank correlation |
| p | dimensionless | two-sided P, t approximation on residual df; empty when n < 43 |
| q | dimensionless | Benjamini-Hochberg q before the post hoc literature links were added |
| KO | id | KO of a post hoc literature link (empty otherwise) |
| source | id | reference of a post hoc link |
| identity_scope | id | identity qualification of a post hoc link |
| direction | id | stated reaction direction of a post hoc link |
| q_augmented | dimensionless | Benjamini-Hochberg q in the augmented family (layer, normalization, model); the q used in the manuscript |

## all_EDS_associations.tsv
sha256: dcd15919dde69cce09d3fdd3e96c5517e35f1786b824bb5575ca7dc0131008dd
rows: 49798 (one row per feature, normalization and model)
produced by: predict_eds_screen

| column | unit | meaning |
|---|---|---|
| layer | id | Metabolite, KO or Pathway |
| normalization | id | native or annotated_relative |
| model | id | unadjusted or host_depth_antibiotic_clinical |
| feature | id | feature id |
| n | count of patients | patients with a positive finite value |
| rho | dimensionless | Spearman or partial rank correlation with EDS |
| p | dimensionless | two-sided P; empty when n < 43 |
| q | dimensionless | Benjamini-Hochberg q within (layer, normalization, model) |

## selected_pair_summary.tsv
sha256: 42e89ff32dc0a4d3977feec112f3e1e827cb5940d82df4062d45df7f9cb8e8c1
rows: 3 (one row per displayed pair)
produced by: predict_eds_focused_pairs (kept for comparison; the notebook recomputes every value)

| column | unit | meaning |
|---|---|---|
| metabolite, KO | id | pair |
| metabolite_label, enzyme_label | id | names |
| reaction, reaction_url | id | KEGG reaction linking the pair |
| n | count of patients | pair-complete patients |
| rho | dimensionless | Spearman rho |
| full_family_q | dimensionless | q_augmented, native unadjusted KO family |
| joint_adjusted_rho | dimensionless | partial rank correlation, joint adjustment |
| joint_adjusted_full_family_q | dimensionless | q_augmented, native joint-adjusted KO family |

## separate_adjustments_all_models.tsv
sha256: fe72b674538a021d48f14fc66c50d8aaa9e454960ca636969276b99a2133917b
rows: 192 (3 pairs x 2 subsets x 2 normalizations x 16 covariate combinations)
produced by: predict_eds_focused_pairs (kept for comparison; the notebook recomputes every value)

| column | unit | meaning |
|---|---|---|
| metabolite, KO | id | pair |
| scope | id | pair_complete or common_84 |
| normalization | id | native or annotated_relative |
| model | id | covariates joined by '+', or none |
| n | count of patients | subset size |
| rho | dimensionless | partial rank correlation |
| delta_from_unadjusted | dimensionless | rho minus the unadjusted rho of the same subset and normalization |

## producer_plotted_ab.tsv
sha256: 80ffe71be9cfe29bc8853db569ced7b43ff002aaaa43e14339997ba3cca19512
rows: 2125 (heatmap cells: 6 x 85 in panel a, 19 x 85 in the former panel b; v30 draws only the 510 panel-a cells)
produced by: the original Supplementary Fig. 10 producer (header and panel a and b rows of the Source Data sheet Supplementary_Figure10_plotted_values.tsv, copied byte for byte; the track rows, which carry per-patient antibiotic state, are withheld)

| column | unit | meaning |
|---|---|---|
| panel | id | a or b |
| row | id | heatmap row, 1 = top |
| label | id | row label as drawn |
| patient_position | id | as in patients.tsv |
| z | dimensionless | within-feature z score; empty where not detected |
| plotted | dimensionless | value drawn, z clipped at ±2; empty where not detected |
