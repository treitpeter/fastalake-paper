# RECEIPTS

## Original Slurm jobs

| role | job |
|---|---|
| baseline | 6044302 |
| refit | 6048399 |

## Code

| original (LAB_ROOT-relative) | file here | sha256 original | sha256 here |
|---|---|---|---|
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/crossomics.py` | `crossomics.py` | `5d0322e7f1fbe3c4f73001a2840d6fae1572035fe69d97dbb8f9072342e86404` | `5d0322e7f1fbe3c4f73001a2840d6fae1572035fe69d97dbb8f9072342e86404` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/provenance/check_prevalence.py` | `check_prevalence.py` | `b05f65dac2aaf53f0f4d238df674e90d17d7632363e3c8d74bc5fcbabc352c10` | `257c4e6968f0292e08cfa7fd0917f5908501e4264b2fd1640bca7cd32a035f22` |

## Inputs (sha256 recorded by the original run)

| input (LAB_ROOT-relative) | sha256 |
|---|---|
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/CAPSCAN_noclust/host_representative_intensity.tsv` | `4f0ee2b452a1b957bcf515c47d9602d2917aa6a7a911cecd725949fcba8539a7` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/CAPSCAN_noclust/microbial_KO.tsv` | `07b03f2a00328458fae14893409274a77eb30a1e69b736acd4f33f58e9df1eaf` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/CAPSCAN_noclust/microbial_OG.tsv` | `14d0eb1f7362e7650c41501385235b71b9c0ed21bc9816e49e2e38fb678dee29` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/CAPSCAN_noclust/microbial_Pathway.tsv` | `972e638c9bba8f4a2e6ddacc625e1011c950ed4e0d5d08150f81dec45c0ebdd6` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/PREDICT_noclust/host_representative_intensity.tsv` | `bf78627a6b6657263a4d52210f39bd966a90de096616ddf9a8f2f94e543fd1cb` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/PREDICT_noclust/microbial_KO.tsv` | `474088cfdbb5688b45e10fefa6007fab55a58b0d90b83dbd79fdcd0db108d672` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/PREDICT_noclust/microbial_OG.tsv` | `61df97764f988071c3bc2c1ffd5797f51f1e1460b066bf60a19603731b303b22` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/PREDICT_noclust/microbial_Pathway.tsv` | `51bfe74644fcb47ae5029f8a544748adca5789588b69401b84b346dc416fb61b` |
| `analyses/capscan_metag_analysis_2026-09-10/aggregate/metaG_KO.tsv` | `17b81e780ba40fb7d90993650f2cd5fa3194a96f945426fc3e59296d2414feab` |
| `analyses/capscan_metag_analysis_2026-09-10/aggregate/metaG_OG.tsv` | `82bc95cf8f010a8aca6caa93a5b0963fcffa8e5ffd5cb404bea439f46748279c` |
| `analyses/capscan_metag_analysis_2026-09-10/aggregate/metaG_Pathway.tsv` | `8d13d06f5399078c848ca50cac1fec1a1ca31148666a6e6022e87ad0c742763e` |
| `analyses/capscan_metag_analysis_2026-09-10/manifests/runs.tsv` | `0edba6053720e2ac1ccb81f15283aed1bf927416ecf6b772cc90951f63a96306` |
| `analyses/enterotypes_2026-09-07/CAPSCAN_metabolite_intensity.tsv` | `1ed653a2d738d6cfb8757a73e472d54f87514577caad29ebcfe038be45a6de13` |
| `analyses/enterotypes_2026-09-07/CAPSCAN_specimen_design.tsv` | `75bfffdb1007257acdf4f2c6b632356f24cf06074590231f6967786876d7ed15` |
| `analyses/enterotypes_2026-09-07/PREDICT_fullstack_primary_samples.tsv` | `5e967b0aa182f2dfd6b2d902d2d1b393dd96caf43ef619575ee5c5b972aa09a7` |
| `analyses/enterotypes_2026-09-07/PREDICT_metaG_KO.tsv` | `077b3f48920f8041e579acaf80d0ff285de4dc4188e8161593ebaf00ebcef034` |
| `analyses/enterotypes_2026-09-07/PREDICT_metaG_Pathway.tsv` | `cc692926acf327b0ccc9f9e71489d4a13a2138432e4d4cd6c1b6d4cab89d3d52` |
| `analyses/enterotypes_2026-09-07/PREDICT_metaT_KO.tsv` | `4cb6c0000a46daefeb7073944377d5779efc085b6b5d617851247e9fe7015a51` |
| `analyses/enterotypes_2026-09-07/PREDICT_metaT_Pathway.tsv` | `d8d267bc41a0610247e4c8fbb46e5c2fdf4dc131a31072cdfc306ccfc5fedb85` |
| `analyses/enterotypes_2026-09-07/PREDICT_metabolite_intensity.tsv` | `0663b71545e31e677368c1fcbfe7d2c503e987f4e5a4eff677da6d8848418cc9` |
| `analyses/enterotypes_2026-09-07/metaG_OG_TPM.tsv` | `c39df737891e118b92c3251aca87be433fae886fe16ceca2ad4e8eb55ebe1f73` |
| `analyses/enterotypes_2026-09-07/metaT_OG_TPM.tsv` | `351a6cbeead4b21b7b5aa40e4c27918170f8a6093b2fa76aab8dc3083f792266` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/full/functional_correlations.tsv` | `ea6dc8d014c6a12c396c4b6cb6beac641351e42b98654ae4b4214fcd0d1cd5bc` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/full/native_SOURCE_DATA_correlations.tsv` | `a02fb8ed34d2f9f68771b75aca4db137edd2a1965a138b57b528dcf6e901a693` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/functional_distance_matrices.npz` | `dd9eb6e595e474ae05d57614e6fbfe1567ba7382211019ffc12acff4d615b632` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/native/INPUT_MANIFEST.json` | `437c76019adcf128d7f2bdf5d7caaf97222b76c15bc1d0c1df8bc925d2df3fb8` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/native/native_pyramids.py` | `64ba668b266b9be42bb0a3f273878a2b93daba53022cc3c5a70ab04d4d56c8da` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/native_CAPSCAN_distance_matrices.npz` | `8cbb306cccd34b78639975fb549749c1f909b74aff31f6b8bcc8e51671ad5b5c` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/native_PREDICT_distance_matrices.npz` | `d579ad4338c51486f3a5d0fb3147da003c35cb9d56c19d70465057fefe43cf3e` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/provenance/check_prevalence.py` | `b05f65dac2aaf53f0f4d238df674e90d17d7632363e3c8d74bc5fcbabc352c10` |

## Outputs (sha256 recorded by the original run)

| output (LAB_ROOT-relative) | sha256 |
|---|---|
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_baseline/capscan_specimen_roster.tsv` | `8b430520e7583964e3cb8d45e13373d8c82b014b725a050b487f9b6ac0478b07` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_baseline/signal_retention.tsv` | `7cb9c93c5376e6b0e228b74e73232f6eddcad61ef14fbf2a97e902073857a2af` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_baseline/native_method_comparison.tsv` | `ee93f6c5da509b011960ee4af0c8eae015235de8500ca191b1f2724709e3401c` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_baseline/functional_distances.npz` | `dee3fef696256dbbe1fab68808c310f7e0c86189b58f543162fcc1ef9c0b200f` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_baseline/predict_patient_roster.tsv` | `5e967b0aa182f2dfd6b2d902d2d1b393dd96caf43ef619575ee5c5b972aa09a7` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_baseline/native_quality.json` | `cdc44b930056861f20822b8029f6040d92a567cb1f8c4924a48c45802a80bd15` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_baseline/functional_correlations.tsv` | `ea6dc8d014c6a12c396c4b6cb6beac641351e42b98654ae4b4214fcd0d1cd5bc` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_baseline/functional_method_comparison.tsv` | `52ad8785c6c385592a4a4f6cc8704f2d73e1891e8c6e64a1e724db21a9f02665` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_baseline/feature_coverage.tsv` | `9497200d7a063885837856b5dd959bc041770862608c35a88d58150d97c6a945` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_baseline/native_distances.npz` | `74a1eb121ff85fa8d8c51289e21fd1ffa6d985a26ab523bef61cd438654fbfb3` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_baseline/eligibility_audit.tsv` | `c6f7a4c206c2c43d69cdc5b65d231360c8e01700e3bb23c1ae4faa707c8d858e` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_baseline/native_correlations.tsv` | `5454e4abfb020ad98f43ca9c320df6d22074d708c0d542f815b81f23b474b2ae` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_refit/capscan_specimen_roster.tsv` | `8b430520e7583964e3cb8d45e13373d8c82b014b725a050b487f9b6ac0478b07` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_refit/signal_retention.tsv` | `fe6fade1a1accf479268488ca1ac45eeecdb3478de9fe082e0e321d363760878` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_refit/native_method_comparison.tsv` | `301143b3df1d1befa29fd5deef2ef8a3e656ce0f956aaab24e4a8d5a1326faf0` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_refit/functional_distances.npz` | `a8e0d5a9f7f2aa52f9dfbfad87a6ac720934a5458bb6f341bbabc49bd53648d0` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_refit/predict_patient_roster.tsv` | `5e967b0aa182f2dfd6b2d902d2d1b393dd96caf43ef619575ee5c5b972aa09a7` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_refit/native_quality.json` | `7d050cf90e7e7ca237556e8d5a155a0706e9a15a0ec44b276bb099aa05b2bcab` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_refit/functional_correlations.tsv` | `4c05b200a276d48e3730276f328749bde336191d373bd7a58c609b5aec1276b6` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_refit/functional_method_comparison.tsv` | `a9e90426ebedb70f92db3dcd1f8d006f963a287d9454f9d8bd02747962b82847` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_refit/feature_coverage.tsv` | `51be35073903f16c664058aef0bf0c595f4a6e3fcfec72e7b366e4afcf83884b` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_refit/native_distances.npz` | `8928b6d572273e1763cf3422fcd5fa28fd381fbbfce0e92ab5bd9d06efa65e51` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_refit/eligibility_audit.tsv` | `c6f7a4c206c2c43d69cdc5b65d231360c8e01700e3bb23c1ae4faa707c8d858e` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_refit/native_correlations.tsv` | `dd1ddf3ff79eb8845de0162bcac3f9e5a722826226fba2ba9803728ebffaba00` |

## Path aliases

Lab path segments are shortened in code and tables (`jess_integration_2026-09-29` to `integration_2026-09-29`, `v10_mm_philip_comment_audit_2026-09-22` to `v10_comment_audit_2026-09-22`, `coauthor_comment_check_v25_2026-10-05` to `comment_check_v25_2026-10-05`). `FASTALAKE_DATA_DIR` must hold the frozen inputs under the shortened names; the original names are the lab directory names.
