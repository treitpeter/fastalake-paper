# Receipts for clinical_second_engine

Paths are relative to LAB_ROOT. Run directory: `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all`.

## Jobs

| stage | Slurm job id |
|---|---|
| manifest build | 6064459 |
| preflight manifests | 6064463 |
| preflight array (CAMPI reproduction, two negative controls) | 6064464 |
| benchmark array (offset 0, indices 198, 255, 512, 551, 888) | 6064480 |
| benchmark array (offset 1000, indices 1228, 1235, 1499, 1517, 1599, 1605) | 6064481 |
| FL PREDICT (tasks 0-295) | 6064500 |
| FL CapScan (tasks 592-1269: FL_public 331, FL_cohortMG 347) | 6064501 |
| FL healthy stool (tasks 1552-1609) | 6064502 |
| MG PREDICT (tasks 296-591) | 6064503 |
| MG CapScan (tasks 1270-1551) | 6064504 |
| collector (RESULTS.md) | 6064505 |
| participant-15 preflight (task 888 reproduced its receipt; tampered FASTA sha256 failed; dry append and end-to-end pilot T760_1799 passed) | 6065736 |
| append participant-15 tasks 1610-1625 | 6065847 |
| FL CapScan participant 15 (tasks 1610-1625) | 6065848 |
| collector FL347 (RESULTS_CAPSCAN_FL347.md) | 6065849 |
| PREDICT mzML truncation check | 6068722 |

The first benchmark submission was rejected (array index above MaxArraySize) and resubmitted with an offset.

## Outputs and inputs hashed

| path | bytes | sha256 |
|---|---|---|
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/RESULTS.md | 2090 | 379e9bbc2119273890af92a3ea07bd88b7b13e47a1a63cfc0140f55e2bfa7f4c |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/RESULTS_CAPSCAN_FL347.md | 1224 | 0a1d66c8eef8f0721fee162072537ae2c3cef6506bad7284ea0da99b43587532 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/results_per_task.tsv | 234693 | d53f8b30a54dfb47fc722b8390005467d21dd21c5e0c2f790c1e0d46de116e56 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/results_paired.tsv | 133068 | 623e005d3df4e036b4b8aa2eae655bb404e7ede6f2cd0ef378818cde03e1bf85 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/results_hstool.tsv | 6278 | d17a10b6ff4852b01800cdc92928f695875b89db8f93dd150879400a8f8a8865 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/results_capscan_fl347.tsv | 41377 | b453aabcb1e909d9483c9c1bc3ffded6e83c5ee4bf2458852eccfd5c48f81933 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/manifest.tsv | 1838058 | 17fdd3afd98ed0a438ae02d9f86b2eb100fbfbaf347dec0fc19a5b97ee520f84 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/MANIFEST_SUMMARY.md | 1699 | 23e9cbd42e9a61739c6f139ca2092d5cd54581bc914964ab2a8c3a5c6be79237 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/inputs/manifest.json | 2181821 | 70607dc98a8fa289ad60b2ff90aa335e03883069ab8fee64c721b557790d7a07 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/inputs/manifest_1610_before_p15.json | 2156546 | ecfe9931e6a14bfdc7e28a88c07408ec8dd87c691bb75cb36934441345ab15c3 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/preflight/PREFLIGHT_RESULT.txt | 662 | d234d872eb6509c104a9fc6d422aa12f6cb04659bceb9e75b20c54ffa1d0e502 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/truncation_check/predict_mzml_endings.tsv | 18796 | 171632dfa5f160d9daf8ee6ca9f21b0564ae0540fecc4dd40982090ac47e5755 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/truncation_check/copies.tsv | 7982 | 499363bc1a51048418c296eb24200ac108af7d251a9dce53429ec7c44c48c464 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/frozen/search.py | 8100 | 3238804ec0ef949b524d7d1212c847c92f5408736f7889da87450223020297bc |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/frozen/lib.py | 3671 | a46c404fc225664624942ef1e380407054bab675ee63bad96fc944f3d74f5736 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/frozen/peptide.py | 6189 | 40e344f8bb3bc900b1528333496a35ec203bd85b79d81d15c339c61446f03eb1 |
| review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/second_engine_all/frozen/closed_fragger.params | 10164 | 774a69b1d086842c4d524dca5427062960078173eb4f3c09df07f0b84358e298 |

The frozen search.py and lib.py hashes are of the run copies, which carry absolute interpreter and java paths; the repository copies in campi_second_engine and campi_near_miss_audit differ from them only in those paths. peptide.py and closed_fragger.params are byte-identical to the repository copies.

Per-task receipts: results/<task>/RECEIPT.json (1,613 PASS) in the run directory.
