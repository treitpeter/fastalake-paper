# Receipts for peptide_length_window_recounts

Paths are relative to LAB_ROOT. Recount outputs are under LEN950 = review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/len950.

## Jobs

| stage | Slurm job id |
|---|---|
| clinical own vs common FASTA recount (10-task array) | 6072755 |
| healthy-stool common-FASTA recount (srun) | 6072769 |
| other recounts and renders (srun, ids not recorded) | not recorded |

Upstream jobs referenced by the recounts (not part of this pipeline): 6048637 (corrected joint search, joint_fileq extraction), 6048638 (compare_fileq.py), 5928104 (A21 comparison), 6057982 (CapScan region test).

## Inputs

| path | bytes | sha256 |
|---|---|---|
| MM_review_analysis/end_of_day_2026-09-15/common_cache/comparison/inputs/manifest.json | 910088 | 06ea28bc9560a0f6bb0ce98b2f32057bcbfc66aee7f64588012fb9f2b681f2d5 |
| MM_review_analysis/A21_7_hstool_quant_comparison/inputs/manifest_revision2.json | 867050 | c29aaba72d84efc9462b64d87e9c75819b88c8ae9e4ce10b19945998f0d93da3 |
| MM_review_analysis/A21_7_hstool_quant_comparison/results/joint_corrected/accepted_identities.tsv.gz | 15371689 | cf436eb977a11c12721074e5649714a118e5dbbafb79d2b5ec798abf46054603 |
| MM_review_analysis/A21_7_hstool_quant_comparison/results/original/accepted_identities.tsv.gz | 11925775 | 98a1d7d0cbc23f1cc0edd996beae28bf7ffdab2e961f52e485b2c78faf84bc2b |
| MM_review_analysis/A21_7_hstool_quant_comparison/results/union_individual/accepted_identities.tsv.gz | 13073387 | 685d066002a447eaf22555d82d781210336364148f554de916aa230db1bb2891 |
| MM_review_analysis/A21_7_hstool_quant_comparison/results/joint_corrected/same_file_id.features.npz | 20669357 | 57f37dfcfd55cd39b73f4edd0e89236cd552d5dd88e9d422937c236991fa559f |
| MM_review_analysis/A21_7_hstool_quant_comparison/results/original/same_file_id.features.npz | 20571500 | c39cfb43cc767dc9a677b7969557c50f3c1814cbd7b93f879aad1a353ba36281 |
| MM_review_analysis/A21_7_hstool_quant_comparison/results/union_individual/same_file_id.features.npz | 21666407 | e4c10a3c2bd06c2da083f6241b0a5a6a235ad72d9d77f040595e4b6317a3d750 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/out/joint_fileq/accepted_identities.tsv.gz | 13065426 | 9819835bd164048bf8d982139ee08adb2f2b502a0e1b40d0b93b1e70a65acc29 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search/out/joint_fileq/same_file_id.features.npz | 17598391 | dd689a4aac4b306a2ee8e255f8de3b64b725f285790a403f11bc61cbb5cea0bf |
| MM_review_analysis/recovery_2026-09-15/public_comparison_recovery/A16/inputs/manifest.json | 3056551 | 7dba4e517eb0fbef5eaaa0b2bd7e1859bd00fa4829718e08e9ef8c3ec9c28092 |
| reviews/supplementary_completion_2026-09-13/assignment_exports/PREDICT_296_MLAB/SEARCH_INPUTS.json | 346741 | 79ff5198f873b9c4049ff0f9539929db17ac8129406bc939b6b9b8057fe61200 |

## Outputs

| path | bytes | sha256 |
|---|---|---|
| LEN950/Fig2/NEW_VALUES.tsv | 18396 | 27035ec3c9a6095384e7cbff73321fcf28fdb7c685aceee95e447d14470bd084 |
| LEN950/Fig2/upset_intersections_FL_9_50.tsv | 11715 | 1f7e4d0e4b01a59eb11a9059641a26274134e48adf83f85db032092c3004035b |
| LEN950/Fig2/work/control/BUILD_LOG.json | 5573 | 6045b115f3db28a105755beb9e575f168a65422a9999b0c5c8d93dace8e81b4a |
| LEN950/Fig2/work/len950/BUILD_LOG.json | 5341 | efc5933b8179ced56e3bf101d227e51bee92547ff65ed231a4a4013176c7ace2 |
| LEN950/Fig3/Figure3_len950_receipt.json | 3080 | 1ffdf358fa226596b4bb174f1cf572675c190644bddd09db867eeca3660aea2e |
| LEN950/Fig3/NEW_VALUES.tsv | 7574 | 95efc2dc00ba99e66a8c0b454fdf7c50810fe31d086378dae11c1f4eacd769d8 |
| LEN950/Fig3/work/control5/Figure3_control5_receipt.json | 3081 | 57296064666397c0e606e057fa9809a977774ec8a437d5d56a9e234a69a9911e |
| LEN950/Fig4/build/capscan_regions.json | 919 | 456785812f5ed7630f82c76ae99f1c24b1f60edc06f1c4762790b21985f41313 |
| LEN950/Fig4/build/same_vs_crossfit.json | 243 | c706adcdbbf78aea8ff52418e284bc025974b024c04f3078bdfcd396fbb85e35 |
| LEN950/Fig4/NEW_VALUES.tsv | 1771 | 376ff86a791e11b28eab474c861e6576df3ec019e0d3657878b15ce9ef677211 |
| LEN950/SFig03/NEW_VALUES.tsv | 2327 | 21c126cff0619aec6dc8b1e470125009ff877407ea5f65f72bf48e427125ff3c |
| LEN950/SFig03/work/control/BUILD_LOG.json | 1705 | c119c611c76a1eb393ec5389d662f3df4e2a9bc4735707aef61742c0058dc6b5 |
| LEN950/SFig03/work/control/plotted_values.tsv | 11133 | 935a51d7b3c4ded186e2147482f301b4c0a4845bb3582539bee24aa8d83790d7 |
| LEN950/SFig03/work/len950/BUILD_LOG.json | 2005 | 6c5f2b0f0459eeaf73d407d6cd9deb94eff67a32450e6fd688bfb23ce7b17869 |
| LEN950/SFig03/work/len950/plotted_values.tsv | 11120 | 137c3fe8ff187884c6be6dfa3ff4cac7deb9dd18496da7c9838fa1545104fc56 |
| LEN950/SFig06/build/render_5_710/panel_de_plotted_values.tsv | 116101 | 45eadc7038017f3ca22147e3d17b63ef8df0c63f56260f873f0c4cb7a784b827 |
| LEN950/SFig06/build/render_5/panel_de_plotted_values.tsv | 116101 | 45eadc7038017f3ca22147e3d17b63ef8df0c63f56260f873f0c4cb7a784b827 |
| LEN950/SFig06/build/render_9_726/panel_de_plotted_values.tsv | 116031 | ac8fcc41cbd6ae4d68eb24d575fddcd5da5361809e2aef4c25fae366e3400628 |
| LEN950/SFig06/build/render_9/panel_de_plotted_values.tsv | 116031 | ac8fcc41cbd6ae4d68eb24d575fddcd5da5361809e2aef4c25fae366e3400628 |
| LEN950/SFig06/build/stats_5_710.json | 2603 | f5d30b6b379c6ded2858f52186970c47771d799a54035911e0c188f1cdb7144c |
| LEN950/SFig06/build/stats_5.json | 2759 | cc37e77c22e74229a106ddc5dda99c636e435bb9c861d85d1cfd4e35706498c0 |
| LEN950/SFig06/build/stats_9_726.json | 2649 | df020cbe515029132755e433298ab6e6c6aff293423174b4565109ce60da163d |
| LEN950/SFig06/build/stats_9.json | 2804 | 25b1bdc400b5bd276374a4894eb7461f11e7f9b3e9296320ff10bf93bbdb9ce4 |
| LEN950/SFig06/NEW_VALUES.tsv | 4237 | e5e3e0a6335b4bb08d566b58395172826b3244c5958645ba29c4f68b0d11b856 |
| LEN950/SFig06/SupplementaryFigure6de_same_spectrum_9_50_plotted_values.tsv | 116031 | ac8fcc41cbd6ae4d68eb24d575fddcd5da5361809e2aef4c25fae366e3400628 |
| LEN950/SFig07/build/plotted_values.tsv | 1725 | 73a5ca8615f8107c00c44005a53e78f6d7ba633b7f510ce34ed6b80ffe2ade09 |
| LEN950/SFig07/build/SFig07_receipt.json | 1478 | ecb5f1305f2268ccbd71d8bb4e7f598c6a83c9501abf0e90d9726ee250e6423b |
| LEN950/SFig07/build/source_data_len950/arm_definitions.tsv | 735 | 18e520c0595b21296a832d96bf09330cff2fba3c5c8ba56b63dffe0248e1c93c |
| LEN950/SFig07/build/source_data_len950/common_search_pair_CV.tsv | 272879 | d0fd92b1481b938e3de5967496b80c47b234a976e16524ba4d240020b2c7f1e2 |
| LEN950/SFig07/build/source_data_len950/common_search_participants.tsv | 31785 | 4d53cd6c590d2f354972ba834c087bee23def3e719c3eefe823734fefe2b672c |
| LEN950/SFig07/build/source_data_len950/completeness_all.tsv | 2104 | 0d73c3788ea060b63c16a068c0b4ba5a74def2f78dc0ef094e8e0e0b3d451f60 |
| LEN950/SFig07/build/source_data_len950/identity_comparison.tsv | 67259 | c5a863acb448d698d0a89b440e5c0a42064e55a5475b0de3d793f6f836287090 |
| LEN950/SFig07/build/source_data_len950/identity_per_acquisition.tsv | 28351 | 95e1ea4bb104f6cd41f978db5ef790cbb0b71bdfde42d0420ead16be03413b67 |
| LEN950/SFig07/build/source_data_len950/identity_summary.tsv | 473 | f4b9f08e71920b6e08aad21436cc3ed487ebf2ae9049048aef1f168602e0aa7a |
| LEN950/SFig07/build/source_data_len950/microbial_group_pairs.tsv | 642027 | e8d58907fb5564149344d76b09e6dbe8b3dd3062e83084ed14da083143471329 |
| LEN950/SFig07/build/source_data_len950/microbial_group_pair_summary.tsv | 4129 | 776895eadb4f32753d6521237ba2322ced6fa7831ed83d33218c274cac3caf28 |
| LEN950/SFig07/build/source_data_len950/technical_agreement_summary.tsv | 2549 | c9d59a24319a1106ab0e1e7895b50133040d78a28faca01918c1895bf7632751 |
| LEN950/SFig07/NEW_VALUES.tsv | 9459 | 8b329f7421f76ca132198108f3efd4dee6317f888610184e9b1dfc24660bea24 |
| LEN950/SFig07/work/clin/part_000.tsv | 5208 | f41afbc3415b206c44bc2d790e57d671cc747e13b514084bd6a9e0796063774a |
| LEN950/SFig07/work/clin/part_001.tsv | 5220 | eadae600ede93dfb5135f3e990c072c7d2d43fcd9ec38cd652ffb286d9b95198 |
| LEN950/SFig07/work/clin/part_002.tsv | 5121 | fa3aa54a793061d451545d4fe49cbefc1d2fffa28bae0b720dcfd4eacca765e8 |
| LEN950/SFig07/work/clin/part_003.tsv | 5200 | f36319742b963fa465f4a1f8f4e78678b4802e089857f2bab98bae5e9b223a21 |
| LEN950/SFig07/work/clin/part_004.tsv | 5216 | 89e6c3126da881f94fe41e568a071da9038af8f6e9e258800e52911e96d9ecaf |
| LEN950/SFig07/work/clin/part_005.tsv | 5204 | 9d67c03292abb430b7b57e1e925429734a6602b0a8bc9a5b611731434f2c9d4d |
| LEN950/SFig07/work/clin/part_006.tsv | 5259 | 628af4d21550e7bac0758fedb18afb702d719f367e2223a57aea3c6efdedfc91 |
| LEN950/SFig07/work/clin/part_007.tsv | 5103 | 831189aad8fd5595c516626a31d58ce0bc6d34160d8d0c156b4405130e25f0fa |
| LEN950/SFig07/work/clin/part_008.tsv | 5072 | f971d0783618c4eb84cec4b8d2cdd0a96f603be57e69a909f912febea0327610 |
| LEN950/SFig07/work/clin/part_009.tsv | 5148 | 64230f3f9562cc1678a506d21318d1041436d72bd0858c0f3463e8f6db8567c3 |
| LEN950/SFig07/work/hstool/completeness_len.tsv | 897 | cb7f9080323d48c26254c6398d68db495e575825f81d136df62a6ad59a797a48 |
| LEN950/SFig07/work/hstool/identity_per_acquisition_len.tsv | 57646 | e2cf4d30b9139ff44f55a908e8402a9008bfad66cc45f551401ecb79be1e78f9 |
| LEN950/SFig07/work/hstool/identity_summary_len.tsv | 599 | 47f2bc3f60d1e0a587eddbdf2e043aebede2dd774adfa8c91eef0c9d13720f87 |
| LEN950/SFig08/build/stats.json | 2888 | 2b633b0115562d7fafc74de21d51bb3c1b64f61c0e7ddd34c7feba4b43c16677 |
| LEN950/SFig08/curation_counts_5_9.tsv | 311737 | faf0e097d3845f94f02caa8ac46a0f01107b2abbd8d6778d27904944f51cb49f |
| LEN950/SFig08/NEW_VALUES.tsv | 2440 | 2b669a010e94de3bdc371188c01a6e7f61709dd3aa177b311e76bffc4d56341f |
| LEN950/SFig09/build/per_specimen_5.tsv | 12909 | 800c0e370c61f8557539c1ec149d981965c49c893797910301108ca49cc144e9 |
| LEN950/SFig09/build/per_specimen_9.tsv | 12931 | b1c12f1d432bc1baf100c7db06abb4114a72a24562561fa7f7a8ab9917188709 |
| LEN950/SFig09/build/stats.json | 1764 | 8a1f942c402995af9fb8c02a900c54c9d2a99105ff27f90753bde57d13d0b18d |
| LEN950/SFig09/build/stats_ylim.json | 325 | 6b67f7ef30d7c6481fe5510816df097c5e0d605e3ab986a3586726042ba17839 |
| LEN950/SFig09/build/tertiles_5.tsv | 182 | 9d9b3e2bf774245e74406e5111e69ad02654eebe4eea746641c48b4c71144796 |
| LEN950/SFig09/build/tertiles_9.tsv | 183 | 761b3e620045ad96c4cdcaefb49ca4eb3bc6b8eeefb48c3e4e4b2719b8f2638f |
| LEN950/SFig09/capscan_acquisitions_9_50.tsv | 28285 | 5bb9ec898849e853e93ef0d01da9b7eaebab64db2b1a075b89f8d9d44fe5ad39 |
| LEN950/SFig09/NEW_VALUES.tsv | 1353 | 06fc1cd2af0fde4a9b1bcb44d0a3b77e6972d3d7621103d092286536b9759766 |
