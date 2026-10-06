# crossomics_distance_sensitivity

Internal ID: A85 (previous folder name A85_controlled_distance_sensitivity)

Does the agreement between layers (Spearman correlation of sample distance matrices) depend on the distance recipe? 220 comparisons are computed over three feature variants, two metrics and, for CapScan, region adjustment; 176 are the controlled estimates. PREDICT rows come from the production run; CapScan rows come from the 347-acquisition rerun that adds participant 15. Used by SFig11 e,f and supplementary methods 8-11.

Status: code copied as it ran (path edits only), not re-run in this repository. See RECEIPTS.md for original job ids and hashes.

## Inputs

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

The 347 rerun reads configs written by `make_configs.py` with the sha256 of every input (C2 groups, specimen design, roster, DNA native cells, extended metabolite table); see `JOBS`/`STATUS.md` of the original folder.

## Outputs

`distance_refit/controlled_176_comparisons.tsv`, `PREDICT_method_comparison.tsv`, `CAPSCAN_method_comparison.tsv`, `*_distances.npz`, `PREDICT_quality.json` (production); `controlled_176_comparisons_347.tsv`, `CAPSCAN_controlled_331_vs_347.tsv` (347 rerun). Hashes in RECEIPTS.md.

## Slurm resources and run time (original logs, `sacct`)

| step | job | cores | requested mem | peak RSS | wall |
|---|---|---|---|---|---|
| baseline | 6044323 | 2 | 12G | 3.5 MB (sampled) | 26 s |
| refit, failed first attempt | 6044327 | 2 | 12G | 3.4 MB (sampled) | 12 s, FAILED |
| refit (production) | 6044329 | 2 | 12G | 3.4 MB (sampled) | 19 s |
| make_configs | 6064852 | 2 | 2G | n/a | 2 s |
| control331 | 6064853 | 2 | 8G | 2.40 GB | 34 s |
| run347 | 6064864 | 2 | 2800M | 0.14 GB (sampled) | 32 s |
| collect | 6064865 | 2 | 1G | 3.4 MB | 1 s |
| negative control, participant 5 (expected to fail) | 6064857 | n/a | n/a | n/a | n/a |

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| feature variants | all observed, prevalence 10%, Wuyts recipe | recipe sensitivity | none |
| metrics | Hellinger, rank (Wuyts recipe rank only) | recipe sensitivity | none |
| adjustment | unadjusted; region-adjusted for CapScan | repeated regions | none |
| comparisons | 220 total, 176 controlled | full grid retained | no MOFA fit, no causal ranking of methods |
| CapScan roster | 13 DNA-complete participants (347) against 12 (331) | participant 15 added | 331 reproduced exactly by control331 (atol 1e-10) |
| PREDICT | production refit unchanged | no CapScan input | none |

## Notes

Code edits relative to the original: lab absolute paths replaced by `FASTALAKE_DATA_DIR`; lab-tree roots located by a marker file replaced by `FASTALAKE_DATA_DIR`. No statistics changed. The scripts keep their original integrity gates (`verify()`, `gate()`, hash checks against `*_CODE.json` and `*_INPUTS.json` manifests of the frozen run). Those manifests contain absolute lab paths and are not in the repository; the hashes they held are listed in RECEIPTS.md. Because the code here carries the path edits described below, its own sha256 differs from the original (both are in RECEIPTS.md), so a re-run needs the manifests regenerated for this copy.`common.py` here (the one of the original `biology_refit/`) provides the shared helpers and must be on the work path. `native_pyramids.py` (A86) and `distance_sensitivity.py` are byte-identical to the copies in the 347 rerun folder; one copy is kept. The mem and run times of the 347 rerun are those of the original `STATUS.md` and `sacct`.
