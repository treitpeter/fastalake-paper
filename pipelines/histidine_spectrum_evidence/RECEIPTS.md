# Receipts

Original Slurm jobs, measured resources (from `sacct`, via `sacct -P -j`), and sha256 of every output file the figures and tables can draw on. Paths are relative to `LAB_ROOT`. Per-acquisition directories that carry a `COMPLETE.json` are summarised by the sha256 of the sorted list `<sha256 of COMPLETE.json>  <path>`; each `COMPLETE.json` itself records the sha256 of every output of its acquisition. Hashing was done on Slurm with `hashjob.py` (see the end of this file).

## Slurm jobs

| job id | stage | states | CPUs allocated | memory requested | wall time max (min) | peak RSS (GB) |
|---|---|---|---|---|---|---|
| 5937565 | repaired_search | COMPLETED: 2 | 16 | 480G | 36.5 | 257.8 |
| 5940795 | case_audit_pilot | COMPLETED: 1 | 2 | 16G | 0.3 | 0.0 |
| 5940796 | case_audit_full | COMPLETED: 1 | 2 | 16G | 0.5 | 0.0 |

## Output hashes

### `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2`

| sha256 | file |
|---|---|
| `9f6debac730c2a8e9b87c4c7db5161aad8a52fe71a272d807602802d1744ee42` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/COMPLETE.json` |
| `a61e2fd0484f296e7274b043a603cd0d5d9374a30b69626fc0f466e66a0cf18e` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_001/PROVENANCE.json` |
| `cf264a7e938a5e26e66af83ff09aed608410191be6dc59b91152c2184e8a2be2` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_001/bounded_sequence_alternatives.tsv` |
| `1953937baf885b00726db73c9c0905832460a87e6d14d797022058abb28d8aca` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_001/matched_fragments.tsv` |
| `ef100abd29ec12cda06a2e4f9eb7995df639e8f18c0e0f43342946a89dece229` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_001/observed_peaks.tsv.gz` |
| `0fcc036e1dfe27dd0af548d04428338033eab9ccdce1d9df4876f71ea08e2193` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_001/processed_peaks.tsv` |
| `ddbd957b040e33c82a8ab4426e1009501e5ef7a2f397cc37a8635bad8fe9d01c` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_002/PROVENANCE.json` |
| `a3b216fec79925ad1bddcad2eeae2be542826638255b1f684b9cfae688843955` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_002/bounded_sequence_alternatives.tsv` |
| `231673024146b08570c86aaddfe077e158a7c864344ccfee2ce42c4d8bf9a780` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_002/matched_fragments.tsv` |
| `f2688c7953b9e0fc68442e1eec382d6bb8a729bc6c05164a69a8c28fd1a23808` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_002/observed_peaks.tsv.gz` |
| `bdc5026cae3d19406dd1659b74151f639f4b83192469b37d5fe542ae564b0c75` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_002/processed_peaks.tsv` |
| `d153a6633ec65f03214a7485f51af4cd1c4eca0e49565c4675c276b2429afd1b` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_003/PROVENANCE.json` |
| `9a66f75c9827cd62c0293e081216e3478f13df273612b46d4b1385d6a2a94357` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_003/bounded_sequence_alternatives.tsv` |
| `8fd892669c8ea3b99db65d3ed07a04a92fb461e924493c4ad1298341b5505ef2` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_003/matched_fragments.tsv` |
| `299f7cbe0680cd074acce30465050a268386ec45f076867ef36205e3cd497acb` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_003/observed_peaks.tsv.gz` |
| `72de668f230552ff1acce0f0b699031113ab5840cfeb1c16ca5908ceadab5fa9` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_003/processed_peaks.tsv` |
| `43814665b3aa225dc05bd0b4a74b24ceaf7d082ae3fab0b7eadf36ffb78a6465` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_004/PROVENANCE.json` |
| `6923dae470c4a2a1a781b9b78f558a24452cd7fec388bfc9a8051923df68fd8b` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_004/bounded_sequence_alternatives.tsv` |
| `c9cf243201f03ad40ee6684ef1276673b4a6a97755b02897702f3a13208bb4ac` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_004/matched_fragments.tsv` |
| `dc064dbcd11e0b8e9df417f8788b835070b158e64b48d959c00119404f7acf98` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_004/observed_peaks.tsv.gz` |
| `d612e0bac96e2806bbb86a1df58b6700f0d13bdfa1117d0fbab93246642948fd` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_004/processed_peaks.tsv` |
| `11bff2af1ea4552805c3d5fd60e0da317d820b68350debb9621a7ac83ae02ef8` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_005/PROVENANCE.json` |
| `a360cc3427aed03dbbdfbc364ffd6b100480a256a845e04521fd2b7678e2cf7c` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_005/bounded_sequence_alternatives.tsv` |
| `2157d0e5bd1231357b0f1f30a595c22ae61dba1d5d9c72f995b9e03d1d9d228f` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_005/matched_fragments.tsv` |
| `432701280f015c5b318199627919bf31f94fd7b446c579c16dcc428bda6b7328` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_005/observed_peaks.tsv.gz` |
| `a9470d5a038276aba0b38f7aeb13f00019a1c3be09c127c0b52863c909c3c461` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_005/processed_peaks.tsv` |
| `06069f7c37a8ebd97072d3d6ac42b487700d81669d02deac8a4f8d0a3cb70833` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_006/PROVENANCE.json` |
| `15b01b2a32e94a05a23171a57b3df30ba4408652ae472b65a799f48db27771ba` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_006/bounded_sequence_alternatives.tsv` |
| `69ab507f3bfbc1304c08e208219d49c475da3a8ae39587d07235b947e06c5ec2` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_006/matched_fragments.tsv` |
| `7e4bf2e6dfe4e0cce9071914f9e1bcc4e59ef090b96b38b298bdc61544367487` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_006/observed_peaks.tsv.gz` |
| `cdb42d462179658168274275a5a6806899bf6fcc41c95ecdef241c4379c8baa1` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_006/processed_peaks.tsv` |
| `82532e8a3d7e1bc3a8bb8fad3381e89d2a76bad2292963342e3b53a9ba783b31` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_007/PROVENANCE.json` |
| `7521657a6bea0dca3f31af96d507e85717376c9f636dda9092df57687866ff59` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_007/bounded_sequence_alternatives.tsv` |
| `65285379d3ea41aa19256b5526433043b4fd9375a5322255927e6fadfa7917b0` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_007/matched_fragments.tsv` |
| `5ac4600376d7f4c47df7270e4bbe63c46dcc52953df394a4777ab9a343460ae9` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_007/observed_peaks.tsv.gz` |
| `0d8ad4c5036b637be731cd32aaa95c9315dfa77c3a9a18797c7e5214e04b9a88` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_007/processed_peaks.tsv` |
| `a71a68a2ceb1dc7ce5d4171a87dacd7d65f1a5f2c86228b9c88f0fd8e69c60dc` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_008/PROVENANCE.json` |
| `6cf8d20a90f5dc897a96fc67d43514b07fca363e5eacbd3fbdab5f1aa892438e` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_008/bounded_sequence_alternatives.tsv` |
| `1ec6bba6b3bafd6568bf529b1b387d57c79d41d668cd530965d89c27e5a4f946` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_008/matched_fragments.tsv` |
| `0cdb909a24d56f8c651af563967cc5ebe1586a5fc1c418c22c11062a085e65e0` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_008/observed_peaks.tsv.gz` |
| `7d09894327f687275d0621d3bfaaaaceb74710d88a99f527959eabed91409822` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/HIS_008/processed_peaks.tsv` |
| `2debd849014c3f2d384028ba45697b7e2f67ae802b919dfe7dc0ffb65a6b3c15` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/REPORT.md` |
| `bd1d022c52cb0ee67969a0654724690af360c20c2db01be5d0e25ce70af4e751` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/all_strata.tsv` |
| `7a38b8669dae897c69a33c86fdaa726759896ca08cf51c937963e617892c0c43` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/complete_candidate_counts.tsv` |
| `4405090ffad09ea88119fadb84b84420f79befce989766c7af71a892b7c2edd5` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/eligible_case_selection_psms.tsv.gz` |
| `933f9320447975a84dd61767958bd73c450b406c9aded6016f26da012eb451dc` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/same_scan_all_search_calls.tsv` |
| `4d7ca8b53f203897ba1a4b24dbc25c89d6aa439279f7471612eb2a7829e19872` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/selected_cases.tsv` |
| `a56ed9a61f5d555374c7df84bf4273131cb94d242d5e3debae98c80238f8492a` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/CAPSCAN/spectrum_audit.tsv` |
| `28b3e770ce24ea8f0487c9b0292a2f016c0be712b8ec922b036a8123229286d4` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/COMPLETE.json` |
| `911c803051057135f128107fd3a4ad0a5fd322b13b1c61113724063b03bbe691` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_001/PROVENANCE.json` |
| `380d365dd0df62b3090b2f857fe941a80eba675e84cfb154d77dd4116eef66db` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_001/bounded_sequence_alternatives.tsv` |
| `86d03b45b9e1e4036bfa851c771f9601c2f8ff2b6e78057c2883054eb094d667` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_001/matched_fragments.tsv` |
| `7bcdf426909b0500e8732807d71c5bc5e38eeeec54fae6c74b5d2abd580417bf` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_001/observed_peaks.tsv.gz` |
| `cf94030301d3409db8abdb9d9a3c0dda95ba1fbf702f18b4a89d91d9f505a2ad` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_001/processed_peaks.tsv` |
| `296721e54f52118e4dc38fd1c547ad67fb076a66756046bbbfc060565f3d9ed5` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_002/PROVENANCE.json` |
| `20c8dc7326d26ec44980a5cb638a4b87f2210e34caa9ad91603a650be24200a3` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_002/bounded_sequence_alternatives.tsv` |
| `dba6b9b5569492f64cf1d833c6dc39ff347e136aca4bc2ef4e47b9cea4fd0eb1` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_002/matched_fragments.tsv` |
| `f8fac8a692e23f9c97d036663ad20fd8adb41bbc39157d6adb6e062b3f16ea7b` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_002/observed_peaks.tsv.gz` |
| `6b4ba4c6adfc34e7ba361107f302114266e6f3394bdde64b19203a4f7fe4bde1` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_002/processed_peaks.tsv` |
| `5c775a3f74301d40b5a693a729939563db3edc715ea891ae38fbcf7dbc02ff2a` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_003/PROVENANCE.json` |
| `be7bdcf9a8170177af6135f26f4416c857caaa0d2ea3ef387163901b9ffedacf` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_003/bounded_sequence_alternatives.tsv` |
| `56502d9ca29d4204e03eb5c26f0b79e3575bca529158880539affd457049c9a9` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_003/matched_fragments.tsv` |
| `d17664261cb1dd917149897870f9e2917ad0c095ce5e8a1408ffda7cfb71e7e6` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_003/observed_peaks.tsv.gz` |
| `efd955b9147a404552d1008e561d0ad935a414a400e15c820ab6e2eae64a54ba` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_003/processed_peaks.tsv` |
| `59019fafe7c33b5b4e989f5616b858ecb5c4d13bf8cd4852f11e08330f2812ff` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_004/PROVENANCE.json` |
| `5ad69e87dfe46b5586150d02b4248b53f4fb210159b63714257fc72c949ac710` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_004/bounded_sequence_alternatives.tsv` |
| `f25f823ea3941ea0cb349a38bc13aa6821cb9c7013dfdf49cfe1b643fe064645` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_004/matched_fragments.tsv` |
| `847dcf73358cd353228c3c193e496365bb5f9b25c673f0cae693fddf39f38e8d` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_004/observed_peaks.tsv.gz` |
| `e442f8260acd5b27ec58c74fc8f694000bdc6237e50dda70ea5b32fdc0ab758d` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_004/processed_peaks.tsv` |
| `417911d7a764207f9a542cb37c0fdd04ada0cde22d44acd74ecc542de910998d` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_005/PROVENANCE.json` |
| `c68e672b19da93a98dc734c694cae1def2b279191d6e11fe8ff9c3147c9c0977` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_005/bounded_sequence_alternatives.tsv` |
| `d28c0a0800ae32d08952e04d2ddda6c0f7d36f70f09df336f475afb44ce4b131` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_005/matched_fragments.tsv` |
| `dbec43d7a0a66c9c3254a96a5cf7e7db66d3cffd413f85c76aac341b9fd665d5` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_005/observed_peaks.tsv.gz` |
| `197654182c3a19e1b6704f61a5616edad45592511b50c4e68c5e3c50634aee04` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_005/processed_peaks.tsv` |
| `6d390bd973df589e1ef449cb1f9dfea24b803c3619258a3bf60e0732f84ec5ff` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_006/PROVENANCE.json` |
| `8284ae2042762a18ecc89245ba634ed7db752d34d4a206543cb3680f824a1bfc` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_006/bounded_sequence_alternatives.tsv` |
| `61b8e74873a78d4b6b00cb828c5318600f8b12ac79317e3b58ef3bae16c61742` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_006/matched_fragments.tsv` |
| `af98781fa4420b8b8da0f23ce0df372848e8e977528c38f6951bf34667f21447` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_006/observed_peaks.tsv.gz` |
| `533b303423d0438a10486d930ed86bf471ccf53e5569edd15f114d66c78bcfa1` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_006/processed_peaks.tsv` |
| `cbac0c3ec94adfee3269e9c073b10ee02f1ecf48e0302528289f87f6918bf9bf` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_007/PROVENANCE.json` |
| `c5a49ff476c86bfe435587a38b6dfce44ddf5253484ce97bcb9bc75ada853650` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_007/bounded_sequence_alternatives.tsv` |
| `5cb09ccaa9352ba648b186958a869e0994626523495f3c57a8e31886ce9c9a9c` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_007/matched_fragments.tsv` |
| `d0e299205e520ddb64702161a0724636321d64e0dd9ed47707c8fe9383484a63` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_007/observed_peaks.tsv.gz` |
| `5529627f507532d0bf733455153a677846a77756f992ae06ad4518b2508a3bf6` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_007/processed_peaks.tsv` |
| `89be6c9a72fc2b5028e4c0a0612de22bf34aded88e22c6ac5fbbc640590c5a53` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_008/PROVENANCE.json` |
| `6e9207a501c2195e22262ccafdfe281b4f364048e0999f7e7e015cbf2eb17d23` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_008/bounded_sequence_alternatives.tsv` |
| `454c735aeb73018dd2ecac5b4aab3a68d7c8d57a6cc3b32aad143f036b64b190` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_008/matched_fragments.tsv` |
| `f4ab83666a64499c82d8765509366b3f34987a7b9116aada984dcea7cb176023` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_008/observed_peaks.tsv.gz` |
| `cb14d6f8c847d640e826dfd75908404863e8f524e63e96732f595d8ffd4d328b` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_008/processed_peaks.tsv` |
| `993809991399084d76959dfbfe6161e74b943216be85dd6a49847c05cf33f8f4` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_009/PROVENANCE.json` |
| `97a064352fb854a6eb0d0a21e5402c1cf1c10281fa39a36d7a5945a3acb5c270` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_009/bounded_sequence_alternatives.tsv` |
| `3e7bd994af3c8983873bcd55455ef8815ca47bc6bad8edd4d2230958dbe283f8` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_009/matched_fragments.tsv` |
| `331e2c953babc9f0fce4ca58a16ca8cc50fe6e40fe7804e95f89ae56ea991a19` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_009/observed_peaks.tsv.gz` |
| `f1b9c5f98cc9f3c8e5e266638b31ae65e9cc14bbcd53104b08f602440fcdf649` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_009/processed_peaks.tsv` |
| `f857065f31147439888a718ef0ab73804ba7a4f1322a4be07e0be379792b1a83` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_010/PROVENANCE.json` |
| `418fb5a0c68df2fbe48692bf3ba55ad34c3b8c432f88d004d98e5dc2e64de90a` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_010/bounded_sequence_alternatives.tsv` |
| `03bbfe569b95ffffdf64d23f8a117f47e7b817f9006de03ca8ed7d343752acf2` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_010/matched_fragments.tsv` |
| `d8675da0b572e5efe82d29868e66e09e024b1bcc6c60ba6070d472a14d9b5f7a` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_010/observed_peaks.tsv.gz` |
| `d7d3b799b70cd6ca6b19369ded753d6f31e1059978d37d7d0f5a1ac0aaef0775` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_010/processed_peaks.tsv` |
| `e4b72e1746f13cf522ca6b7b746fa32c8927f57a2f0fedaf59720550dce9261f` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_011/PROVENANCE.json` |
| `97e750e3ef875d16df845d559857d226b6aa142e15652b43c95e8b9ddd015bba` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_011/bounded_sequence_alternatives.tsv` |
| `5c19014c9d4b49adba1db34fa8dcd5e6e9747b7e2d8fef2eabc661aa6bbc5a4a` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_011/matched_fragments.tsv` |
| `8837614afa4a93f1a8b4ec7406318ca388061ea58dabbfb89ba171ea13e5539c` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_011/observed_peaks.tsv.gz` |
| `ea48941329059fe3e3c4491ca56057ce8ee5297b834c91c8a2b6a6b96958923c` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_011/processed_peaks.tsv` |
| `84b5ee1bd34e8bf098e798e8fed98f84f169b79382d07598d0a101374cd2542c` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_012/PROVENANCE.json` |
| `911a8f23e5337881b35d4396b98719521dcc2205698d862ffbf6e90c58fe1925` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_012/bounded_sequence_alternatives.tsv` |
| `17b95e38c4adc5cc61ae16735d29b29d8936e9ee9a1baaaedf61365c8f792efd` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_012/matched_fragments.tsv` |
| `6c8993b38346e285b1da834a994aeb1a83f8ca01e625ffa83c2265d972a07ca3` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_012/observed_peaks.tsv.gz` |
| `e987bf3639c68d1b1a7bb9c22790787372a6e654c794c47db6e94bbd065313ba` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/HIS_012/processed_peaks.tsv` |
| `2debd849014c3f2d384028ba45697b7e2f67ae802b919dfe7dc0ffb65a6b3c15` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/REPORT.md` |
| `ac9de48f223a7e7cd2cad4cdb2e0f83f79c14c59f051557e588f32cf7a5f1d05` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/all_strata.tsv` |
| `8c3f8e18e002b14984aee08d5807beea49c1c490856ba0d0d9b1cc1bc52cfa31` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/complete_candidate_counts.tsv` |
| `0a6e0c6b670dcbb6d59251878fba0fd50cd656dccb155575a1faeced4487eb5f` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/eligible_case_selection_psms.tsv.gz` |
| `a597266166854e4a6269f72640dceb5c036d760ed2ca083b1f5d13601668e945` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/same_scan_all_search_calls.tsv` |
| `ad4a2e0e82ac49ee98fcb0e4867a9e4dc90348ccad1f0c50ec7adf2a082be104` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/selected_cases.tsv` |
| `728409e2538978f61388bee2dfb8ab0a828a6f13eef6e82a4c735cdad37a310c` | `LAB_ROOT/MM_review_analysis/A25_4_7_spectrum_evidence/results_revision2/pilot/spectrum_audit.tsv` |

### `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2`

| sha256 | file |
|---|---|
| `db4ceae61845b67592de24e8a280834a678ff9e5def51a090533fa514164edc0` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/COMPLETE.json` |
| `3f7da0215372ee5180859e604f131dab227f1e62ea247691aa1b48d37010dd6d` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_001/PROVENANCE.json` |
| `faa42f67f6564d18363a3615c4a1e6ff10658c5ea233a87885fc7b19b877cebf` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_001/bounded_sequence_alternatives.tsv` |
| `ad8564088d3e603c0f76428c8d39d0067909e1c5bc2d18f631ed4fbb3e0c3b1a` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_001/matched_fragments.tsv` |
| `24587300e52c1a504dd5251b77966a0e86b1f441f87869825614f09accddce35` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_001/observed_peaks.tsv.gz` |
| `524424e9b8c8eea1921768c4717476b7cd171700b4504465ae2d25a1271de64c` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_001/processed_peaks.tsv` |
| `6facb3e4eae5ec490d4e906464fb11e057ca6c0d7fc2964e7ab590eaa8a0df51` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_002/PROVENANCE.json` |
| `3780d3648749517cff8f7d750260f2dacd7abc09330adeda090abc09107da7b7` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_002/bounded_sequence_alternatives.tsv` |
| `e52dacaad2807aebd5ef9526eae7cfa34757727384883e763633b7818d32cc4c` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_002/matched_fragments.tsv` |
| `38cf0d2b2a52473898e0fc01cf70ecde422e91e80b5e225e0a300e6928158350` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_002/observed_peaks.tsv.gz` |
| `062827e3579abdba1869585faf0ff08e258deb57ca04ba1cab3d54ae27b5e1c7` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_002/processed_peaks.tsv` |
| `4d8d5dd6aa08820fb5e10b79d55322d4b83b85f6fe2efdbe74dd4932d3654743` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_003/PROVENANCE.json` |
| `cf0fdf3f4be75168b80e391e4db582e21b848587c496c3b849d3e7b462c8d70d` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_003/bounded_sequence_alternatives.tsv` |
| `ce847bb980a96fbcc111b7a9513e6931dfdaa381e036abb8771694e323324e08` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_003/matched_fragments.tsv` |
| `f9ab97b74b14790829b35467627392520d8a5dc226da2c45cbc923dd63914e3d` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_003/observed_peaks.tsv.gz` |
| `2062af57461fe3768b9e7c15e58a2d0a1a427ee480ec558b96bd57869d247342` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_003/processed_peaks.tsv` |
| `2dde24f5fd1e63c3139c90ed9647b046ca61d9c9b34b69b80985dea14188c537` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_004/PROVENANCE.json` |
| `8e9ed927e00b275fb618af2bc1769dcbab787318bbab9f1f09e8911de6051c31` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_004/bounded_sequence_alternatives.tsv` |
| `d68208d35f0229fb18c177bfc6bd189173d3dbdf429ca015cd2a4ba20a2d4649` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_004/matched_fragments.tsv` |
| `4b2741259600b9f5813a629293619e035010b63d172559626a85f3beeacc63df` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_004/observed_peaks.tsv.gz` |
| `0bb92c37371cd69c80e38c678999b84a1df3e3e0b2b884d8e0553cc1ecd12585` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/HIS_004/processed_peaks.tsv` |
| `2debd849014c3f2d384028ba45697b7e2f67ae802b919dfe7dc0ffb65a6b3c15` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/REPORT.md` |
| `8d1f3263b5219deb52b06e07ecb8dace03091bd9c92fe2847795db290350e328` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/all_strata.tsv` |
| `0280c6d595144d9125385c717e3c35b096fa632fef42adf80387e8fd5497c706` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/complete_candidate_counts.tsv` |
| `0e9870476d7bbd556bd269b0506945b93c920adb80e8b6d047364ff5a0df0a1d` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/eligible_case_selection_psms.tsv.gz` |
| `c41eaa2bea0d0297292dca51688b535ead9acf8e46cb0e9e954c8b8e333618b5` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/same_scan_all_search_calls.tsv` |
| `fb903709852b39041f92a55beeac3cf899870d27688b78e2baefea275318f696` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/selected_cases.tsv` |
| `158abf414b43037e1560972bb97e31241aff2f85e44480b557327dd678570963` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/spectrum_audit.tsv` |

### `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/inputs/manifest_revision2.json`

| sha256 | file |
|---|---|
| `63f27def70f1f05a3ac688d56ede5301a437bf5d7c9ae3e93afcfbb784be8cfa` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/inputs/manifest_revision2.json` |

## Released code versus origin

| file | origin (LAB_ROOT relative) | sha256 at origin | sha256 here |
|---|---|---|---|
| `search_revision2.py` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/frozen/search_revision2.py` | `e8593f5b24caccaf7f6affb9126a3ac6f3ec0825b88de7250968cea75db3b7a2` | `e8593f5b24caccaf7f6affb9126a3ac6f3ec0825b88de7250968cea75db3b7a2` |
| `fragment_helpers.py` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/frozen/fragment_helpers.py` | `d0dca418ddbd95084a920ff5928d06dacd76341cd9530c8f9f33ed8358c4cd4e` | `2d8db5e3e9daba24c8b61d076e2d01028a72eb974667c31bbdff2872bf0782d0` |

Static checks: `python -m py_compile` on every `.py` file; `tools/scan_release.py`.

## Hashing code

The hashes above were produced by a Slurm job (4 to 8 CPUs) that walks each root once, hashes every file with SHA-256 in 8 MiB blocks and, for roots with more than 20 per-acquisition directories, hashes only the `COMPLETE.json` of each. The code is kept in the work directory, not in the repository.
