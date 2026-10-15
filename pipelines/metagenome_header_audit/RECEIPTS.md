# Receipts for metagenome_header_audit

Paths are relative to LAB_ROOT.

## Jobs

| stage | Slurm job id |
|---|---|
| pilot | 6025421 |
| full PREDICT | 6025426 |
| full CAPSCAN | 6025427 |
| figures | 6025428 |
| CapScan 295: manifest build | 6085728 |
| CapScan 295: negative control (expected failure) | 6085729 |
| CapScan 295: 282 control | 6085730 |
| CapScan 295: full run | 6086073 |

worker.py sha256 331ad8e43660a0e4a25d9d27cbd80b577d296b46519aca02a22ad4ef997a7fe2 (= code_sha256 in COMPLETE.json).

## Outputs and inputs hashed

| path | bytes | sha256 |
|---|---|---|
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/manifest.json | 2658135 | 02b020362f3d65708d5a187effda21e9e83973fd20e07d4420ad6df9152f1249 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/CAPSCAN/COMPLETE.json | 592 | 15446469f952e3e244ff574025323010f00230d61a7ea46fcb6b6d83872127b6 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/CAPSCAN/header_audit.tsv | 755 | d685580bc6a26991ade8ce47b1586defd1e5cf4c956e481b2f942cefde0eda9a |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/CAPSCAN/MG_SEQUENCE_PRESERVATION.json | 31145 | be8af33dd2cb2988d39ac4cbec826571c4d5ac4576872f9bff606d504e94d09c |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/PREDICT/COMPLETE.json | 592 | a61bd6ccd309843cb0869bde68373ee237f05bf1b6416c430c44ffb6080f25ee |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/PREDICT/header_audit.tsv | 815 | c9b2560dc1a42f204701e2fdfdb46a68125e75f7b00ec5354f3d11723972a25a |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/PREDICT/MG_SEQUENCE_PRESERVATION.json | 35050 | eafacec77b8bcdaa7416d03fcc4bd8e89f9edded29b07450ae8e86b598e1a5a0 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/figures/header_inconsistency_audit.pdf | 19310 | fbfa962dc4582e26b4a92e5f573059b78e9e6e55591e9b63318bd99f560bb6f3 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/worker.py | 8399 | 331ad8e43660a0e4a25d9d27cbd80b577d296b46519aca02a22ad4ef997a7fe2 |

| LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item3_collision/lab/review/v15_integration/analyses/run295/manifest.json | 1226992 | 69e1d874396f4cefd7632fb30d54931510ec98fcc800512b7d9f848131798fc9 |
| LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item3_collision/lab/review/v15_integration/analyses/run295/full/CAPSCAN/COMPLETE.json | 592 | 54e11207adc291ace1cafe581f4f8f7fc7de8f490181e09d380c716652978712 |
| LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item3_collision/lab/review/v15_integration/analyses/run295/full/CAPSCAN/header_audit.tsv | 754 | 1341faa59f07d1fa47736a06b4503e7f469fec90438c3152f00afd5d7267e233 |
| LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item3_collision/lab/review/v15_integration/analyses/run295/full/CAPSCAN/MG_SEQUENCE_PRESERVATION.json | 32575 | e0bf88dc14780d5ab9eccb862921b66c7a73ababcd56b01babf1380fe3e2e69b |

Repository copy: source_data/supplementary_methods_1_5/header_audit.tsv, sha256 f31ac81ead3ec695c159c93cde09db9dc1f2f2bce5725e0408984325699085b9 (PREDICT rows from the 296-specimen run, CapScan rows from the 295-acquisition run; the 282 version was f31ac81ead3ec695c159c93cde09db9dc1f2f2bce5725e0408984325699085b9).
