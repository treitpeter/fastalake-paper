# Source data: Supplementary Table 5

Inputs of `tables/STab5/STab5.ipynb`. PREDICT cohort, 77 patients (one earliest complete
specimen each). The first four files are byte copies of the analysis outputs; the profile
files are the patient-level inputs of the correlation step, from which the notebook
recomputes every test. How each file was made from the lab files: `tables/STab5/prepare_source_data.py`
and `tables/STab5/PROVENANCE.md`. No clinical covariates and no raw patient or sample
identifiers are included; patients carry the Source Data package pseudonyms.

Profile values: metaG and metaT OG values are TPM; metaG and metaT KO and pathway values
are the per-sample abundances of the producer tables (unit not recorded in the producer
files); microbial_metaP values are LFQ intensity allocated to the function (equal share
among a feature's distinct eligible terms; human-assigned features excluded, microbial,
other and unresolved origins retained: profile `other_unresolved`), averaged over the
patient's acquisitions. The correlation step normalizes each layer to relative signal over
the function identifiers shared by the layer pair, so only relative values matter.

## crossomics_same_function_correlations.tsv
sha256: bbab94eb3dcbc52e24b841b070316cc5c5614080bcec0138a112659f9c90c18f
rows: 30581 (one row per eligible function and layer pair, revised allocation (sheet All_30581_tests))
produced by: predict_function_correlations

| column | unit | meaning |
|---|---|---|
| level | id | functional representation: OG, KO or Pathway |
| layer_a | id | first layer: metaG (DNA) or metaT (RNA) |
| layer_b | id | second layer: metaT (RNA) or microbial_metaP (non-host protein) |
| feature | id | function identifier (eggNOG OG token, KEGG ortholog K#####, or KEGG map#####) |
| n_patients | count of patients | patients with positive values in both layers (35 to 77) |
| rho | dimensionless | Spearman correlation across those patients |
| p | dimensionless | two-sided P (scipy.stats.spearmanr, t approximation, n - 2 df) |
| q | dimensionless | Benjamini-Hochberg q within the level x layer-pair family |

## previous_allocation_same_function_correlations.tsv
sha256: 3cf4e8e0b22812d156d201ec926c95e53023a7467653363b53ce3dfc657ef8a8
rows: 17914 (one row per eligible function and layer pair under the previous allocation (17,914 tests; input of Eligibility_common))
produced by: predict_function_correlations

| column | unit | meaning |
|---|---|---|
| level | id | functional representation: OG, KO or Pathway |
| layer_a | id | first layer: metaG (DNA) or metaT (RNA) |
| layer_b | id | second layer: metaT (RNA) or microbial_metaP (non-host protein) |
| feature | id | function identifier (eggNOG OG token, KEGG ortholog K#####, or KEGG map#####) |
| n_patients | count of patients | patients with positive values in both layers (35 to 77) |
| rho | dimensionless | Spearman correlation across those patients |
| p | dimensionless | two-sided P (scipy.stats.spearmanr, t approximation, n - 2 df) |
| q | dimensionless | Benjamini-Hochberg q within the level x layer-pair family |

## same_function_summary.tsv
sha256: a7f59dbe3703c0f3d9be55d6d451aca737a813bbbfe04e78977ef713bd106a82
rows: 9 (one row per level x layer-pair family (sheet Nine_family_summary))
produced by: predict_function_correlations

| column | unit | meaning |
|---|---|---|
| level | id | functional representation |
| layer_a | id | first layer |
| layer_b | id | second layer |
| features | count of functions | tests in the family |
| median_rho | dimensionless | median Spearman rho across the family |
| significant | count of functions | tests with q < 0.05 |

## eligibility.tsv
sha256: ccba7b357e04b450e61d586a8bf8934a719a96b0162a2324c9c9cd941d8d2b15
rows: 184386 (one row per candidate function and layer pair (function present in both layers' vocabularies; sheet All_eligibility))
produced by: predict_function_correlations

| column | unit | meaning |
|---|---|---|
| level | id | functional representation |
| layer_a | id | first layer |
| layer_b | id | second layer |
| feature | id | function identifier |
| n_patients | count of patients | patients with positive values in both layers |
| eligible | boolean | True if n_patients >= 35 and the correlation was finite |

## profiles_OG_metaG.tsv.gz
sha256: 18df033bbe840ef2f12b42a1fc13e6a05659be59f9846ead1a42f5e4f0ada0b6
rows: 1140547 (one row per patient and function with a non-zero OG value, metagenome (DNA); 77 patients)
produced by: predict_metagenome_metatranscriptome_functions

| column | unit | meaning |
|---|---|---|
| patient | id | pseudonymised PREDICT patient (PREDICT_Pnnn, as in the Source Data package) |
| feature | id | function identifier |
| value | see below | abundance of the function in that patient; only non-zero values are listed, absent pairs are zero |

## profiles_OG_metaT.tsv.gz
sha256: 5858906c163bf6b8b09d2383cbc9957f914b9aadb946ce6847806e8cafc5a2dd
rows: 597417 (one row per patient and function with a non-zero OG value, metatranscriptome (RNA); 77 patients)
produced by: predict_metagenome_metatranscriptome_functions

| column | unit | meaning |
|---|---|---|
| patient | id | pseudonymised PREDICT patient (PREDICT_Pnnn, as in the Source Data package) |
| feature | id | function identifier |
| value | see below | abundance of the function in that patient; only non-zero values are listed, absent pairs are zero |

## profiles_OG_microbial_metaP.tsv.gz
sha256: 922bf65773de69fa001ede80933ae4bd912a3c14704c0efb52cfe6be5be7d65d
rows: 1260681 (one row per patient and function with a non-zero OG value, non-host protein, mean over the patient's acquisitions; 77 patients)
produced by: predict_verified_functional_profiles

| column | unit | meaning |
|---|---|---|
| patient | id | pseudonymised PREDICT patient (PREDICT_Pnnn, as in the Source Data package) |
| feature | id | function identifier |
| value | see below | abundance of the function in that patient; only non-zero values are listed, absent pairs are zero |

## profiles_KO_metaG.tsv.gz
sha256: 5d579de10bdb89e514e982f716ebcecbf7726f876511992eedecbd5bfd0bcd0c
rows: 381129 (one row per patient and function with a non-zero KO value, metagenome (DNA); 77 patients)
produced by: predict_metagenome_metatranscriptome_functions

| column | unit | meaning |
|---|---|---|
| patient | id | pseudonymised PREDICT patient (PREDICT_Pnnn, as in the Source Data package) |
| feature | id | function identifier |
| value | see below | abundance of the function in that patient; only non-zero values are listed, absent pairs are zero |

## profiles_KO_metaT.tsv.gz
sha256: c7cf1716a11b128600ccef083e04000f382324883930a26b3bf434930b7d2249
rows: 303877 (one row per patient and function with a non-zero KO value, metatranscriptome (RNA); 77 patients)
produced by: predict_metagenome_metatranscriptome_functions

| column | unit | meaning |
|---|---|---|
| patient | id | pseudonymised PREDICT patient (PREDICT_Pnnn, as in the Source Data package) |
| feature | id | function identifier |
| value | see below | abundance of the function in that patient; only non-zero values are listed, absent pairs are zero |

## profiles_KO_microbial_metaP.tsv.gz
sha256: 54314052201f0f9984eda5254145d01e28eb4c8d78bfb55d1efbd1794cd7ac92
rows: 475912 (one row per patient and function with a non-zero KO value, non-host protein, mean over the patient's acquisitions; 77 patients)
produced by: predict_verified_functional_profiles

| column | unit | meaning |
|---|---|---|
| patient | id | pseudonymised PREDICT patient (PREDICT_Pnnn, as in the Source Data package) |
| feature | id | function identifier |
| value | see below | abundance of the function in that patient; only non-zero values are listed, absent pairs are zero |

## profiles_Pathway_metaG.tsv.gz
sha256: aa3ee72655bb382b099c5919f6b4c7d9f46a3156b4e682426c10254f63de8cd7
rows: 21482 (one row per patient and function with a non-zero Pathway value, metagenome (DNA); 77 patients)
produced by: predict_metagenome_metatranscriptome_functions

| column | unit | meaning |
|---|---|---|
| patient | id | pseudonymised PREDICT patient (PREDICT_Pnnn, as in the Source Data package) |
| feature | id | function identifier |
| value | see below | abundance of the function in that patient; only non-zero values are listed, absent pairs are zero |

## profiles_Pathway_metaT.tsv.gz
sha256: 48a9aeb54f565239c22fa96ffea79760b6d244c7a33b5a186cb5faa620768bdd
rows: 19801 (one row per patient and function with a non-zero Pathway value, metatranscriptome (RNA); 77 patients)
produced by: predict_metagenome_metatranscriptome_functions

| column | unit | meaning |
|---|---|---|
| patient | id | pseudonymised PREDICT patient (PREDICT_Pnnn, as in the Source Data package) |
| feature | id | function identifier |
| value | see below | abundance of the function in that patient; only non-zero values are listed, absent pairs are zero |

## profiles_Pathway_microbial_metaP.tsv.gz
sha256: 9284464703a06d68b81c92a707476d580801829c8f8fb17c53bd51c6e3e602fc
rows: 26111 (one row per patient and function with a non-zero Pathway value, non-host protein, mean over the patient's acquisitions; 77 patients)
produced by: predict_verified_functional_profiles

| column | unit | meaning |
|---|---|---|
| patient | id | pseudonymised PREDICT patient (PREDICT_Pnnn, as in the Source Data package) |
| feature | id | function identifier |
| value | see below | abundance of the function in that patient; only non-zero values are listed, absent pairs are zero |

## profile_features.tsv
sha256: b6c2f66c7942a64e2c0428f03be62b596c5b6bd8cb40b1e5c792383536e1586d
rows: 196140 (one row per feature column of the profile matrices (level, layer), including features that are zero in all 77 patients)
produced by: predict_function_correlations

| column | unit | meaning |
|---|---|---|
| level | id | functional representation |
| layer | id | metaG, metaT or microbial_metaP |
| feature | id | function identifier |

## patients.tsv
sha256: b266650fe8e08a82b941f20f16ae97fa0f4888287950e16fd2bead66bb66f5ea
rows: 77 (one row per patient (77))
produced by: predict_function_correlations

| column | unit | meaning |
|---|---|---|
| patient | id | pseudonymised PREDICT patient |
