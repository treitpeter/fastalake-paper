# RECEIPTS

## Original Slurm jobs

| role | job |
|---|---|
| baseline | 6044323 |
| refit first attempt (failed) | 6044327 |
| refit (production, PREDICT rows and 331 CapScan rows) | 6044329 |
| make_configs / control331 / negative control / run347 / collect | 6064852 / 6064853 / 6064857 / 6064864 / 6064865 |

## Code

| original (LAB_ROOT-relative) | file here | sha256 original | sha256 here |
|---|---|---|---|
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/distance_sensitivity.py` | `distance_sensitivity.py` | `f23a5468d394757c4ba73eee691bde3d19f0896a1716dfe6b37f1d7205992d65` | `413464a65aa2d3c2a7d3f49aa30c297736d70bfc966cfd01d423fcf2b53d40b8` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/common.py` | `common.py` | `e7e8564e77b125cebec0bfa4f38e7f6c57b00f07845eb6e13522831caa2fb232` | `b9cafb3902139dcfa45417bbedfe3835283ab4e960a7737e45a2f34c6e7519c4` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/reruns_347/A4_controlled_distance/code/experiment.py` | `experiment.py` | `382d8f70e3352f493e8b71b328d82bd9fa91086d390b7d945ecc26c78eb62d62` | `51a544ebe56f5ab5ad7a99aeea034423e9df99166efd047ad6cce6c5c87e2f9b` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/reruns_347/A4_controlled_distance/code/distance_capscan347.py` | `distance_capscan347.py` | `d81c1293a9da1f045bf66b970716d88383a018e41b4dc2dd084fe363552a9fc3` | `d81c1293a9da1f045bf66b970716d88383a018e41b4dc2dd084fe363552a9fc3` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/reruns_347/A4_controlled_distance/code/make_configs.py` | `make_configs.py` | `6d25a260888bff0e2bb82d7a565fb8c38c6afe934e15ad2a4b6e592dc6b716fc` | `9192c3aa5828da3ab004a73ee166bb1e39bc03fd72c57241b669b4380168e1a8` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/reruns_347/A4_controlled_distance/code/collect.py` | `collect.py` | `36f8cff072f9c03c6ccd8a1410689d7482262181d8e5a1df1ad91c0e701b50bf` | `3e22e5888f1b1b37fbe629e3b2aa20dacb648a96ad2a78824cb52a1da9fb82c0` |

## Inputs (sha256 recorded by the original run)

| input (LAB_ROOT-relative) | sha256 |
|---|---|
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/distances/CAPSCAN/COMPLETE.json` | `fa46832362992abd5be7bf8d6dd90b624ca1082398ba0f9eaf7824d99d1ed3d9` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/distances/CAPSCAN/comparisons.tsv` | `4df675ac2253e21fa4c139fba89090d44289ac47efbd7016f0e6b76e92ca91fd` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/distances/CAPSCAN/distance_matrices.npz` | `fba77197c69d9ea3f3e6e3f79d2fb3d393d2cf2c9b348db8e82c7167ecb52195` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/distances/PREDICT/COMPLETE.json` | `f21489a4e0dccf2b6a24745734c5173ca1dcb4fe9cde724dff63492ee29fc919` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/distances/PREDICT/comparisons.tsv` | `af3f5d49d4cbbc00eaccab82a1625917b8c7d151b1a0b31434c8dc6853d29d5e` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/distances/PREDICT/distance_matrices.npz` | `6a50407e9969278cb17096aa3be887aed8a5e271352e357ea92cadf246217271` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/experiment.py` | `382d8f70e3352f493e8b71b328d82bd9fa91086d390b7d945ecc26c78eb62d62` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/prepared/CAPSCAN/COMPLETE.json` | `d68e746e67a4b6ecbe2e74f7ec08fec750d60ca69de3682c8f54f39a1bcc91d0` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/prepared/CAPSCAN/Host_protein.npz` | `ebda2c13e2c6e92e375e7501a7fc0261f3d9d579781a7d18ac85816b426b05cf` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/prepared/CAPSCAN/Metabolites_log_profiles.npy` | `267037d4712a4d86d7bb2cd2d757bc771031ec64278c6eb44c3c0d436c719148` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/prepared/CAPSCAN/Non-host_protein.npz` | `40576921ebdaab9055fa3fbba14b85e4b3e4e01cd636af48e349762d3d913304` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/prepared/CAPSCAN/profiles.tsv` | `2b1e4d70240d09b52cadf6d4fc9d718c2c3ad01d44a35d9344f41dd737f20b79` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/prepared/PREDICT/COMPLETE.json` | `1764b27649e51ecdfc18ce9b3b0900b5e167f057c8118361b8d88a7d30007605` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/prepared/PREDICT/Host_protein.npz` | `bd9440d21c6ea516e59f99b026043b9370b502ebc088c5a9d091b88a8e714814` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/prepared/PREDICT/Metabolites_log_profiles.npy` | `d28bcd777be65ee549597e8c99ecc5f005395d6176c8c301540f88c7698310b2` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/prepared/PREDICT/Non-host_protein.npz` | `90cf6279ded7313001e7577def282453aba8a806c10deed04833a9263634caf4` |
| `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/figure5_clarity_review_2026-09-23/native_pyramids_2026-09-23/method_comparison_mofa_2026-09-23/prepared/PREDICT/profiles.tsv` | `0a675328294fcab9f2f59b9402558f7e5648b058451b7421bb3bb216480d61cf` |
| `analyses/enterotypes_2026-09-07/CAPSCAN_specimen_design.tsv` | `75bfffdb1007257acdf4f2c6b632356f24cf06074590231f6967786876d7ed15` |
| `analyses/enterotypes_2026-09-07/PREDICT_fullstack_primary_samples.tsv` | `5e967b0aa182f2dfd6b2d902d2d1b393dd96caf43ef619575ee5c5b972aa09a7` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/cross_refit/capscan_specimen_roster.tsv` | `8b430520e7583964e3cb8d45e13373d8c82b014b725a050b487f9b6ac0478b07` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/native/native_pyramids.py` | `64ba668b266b9be42bb0a3f273878a2b93daba53022cc3c5a70ab04d4d56c8da` |

## Outputs (sha256 recorded by the original run)

| output (LAB_ROOT-relative) | sha256 |
|---|---|
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/distance_refit/PREDICT_quality.json` | `d0c361dfd2bdd81135fd60e5fb167e48cb47df869bd92cac6ea2402a298a98b7` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/distance_refit/controlled_176_comparisons.tsv` | `380f1fb46ce2921166baa70ffcc929fe43c4410c8c4ae15e7cf5464be2a916d2` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/distance_refit/PREDICT_distances.npz` | `f596d528271810b8ad604b48c5a1f1b2f218870ba3b76d83d1c8661c2ed4c7ee` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/distance_refit/PREDICT_method_comparison.tsv` | `0f21df6acc1854aea3b630b42623121b92016d12afcda9d6ba4565c01a7d449e` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/distance_refit/CAPSCAN_method_comparison.tsv` | `39b03ab50e9ff340a31f164c9dedbc57e6f336538f446cbeb5b8cbe23bb6e746` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/distance_refit/CAPSCAN_distances.npz` | `58c350c32f140762e535f12d2cb03e8ac6a27ecba986ccff2d4d372a8f13d74b` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/distance_refit/comparisons.tsv` | `8d5e5671532b046ea2b178f3b0a5fdbca301d461db91ab79b7e7ceac28ba2045` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/distance_refit/CAPSCAN_quality.json` | `60ca1a55e1d33c002cf6b4dc76d3dbe66f9e3be7cdd9c29909d7289c3af6609e` |

| output (LAB_ROOT-relative) | sha256 |
|---|---|
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/reruns_347/A4_controlled_distance/controlled_176_comparisons_347.tsv` | `e0757c104a4e14526ee59f1fe237a4d43b546c4577250a17e28aef6d7dad289d` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/reruns_347/A4_controlled_distance/CAPSCAN_controlled_331_vs_347.tsv` | `c6c9ce726635719efdaf48c6fce11b7fe6924aaa832891d5929d65577fbc357f` |

## Path aliases

Lab path segments are shortened in code and tables (`jess_integration_2026-09-29` to `integration_2026-09-29`, `v10_mm_philip_comment_audit_2026-09-22` to `v10_comment_audit_2026-09-22`, `coauthor_comment_check_v25_2026-10-05` to `comment_check_v25_2026-10-05`). `FASTALAKE_DATA_DIR` must hold the frozen inputs under the shortened names; the original names are the lab directory names.
