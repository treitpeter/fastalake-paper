# Receipts for cohort_database_accumulation

Paths are relative to LAB_ROOT.

## Jobs

| stage | Slurm job id |
|---|---|
| pilot (three specimens per cohort) | 6025421 |
| full PREDICT | 6025426 |
| full CAPSCAN | 6025427 |
| figures (afterok 6025426:6025427) | 6025428 |

worker.py sha256 331ad8e43660a0e4a25d9d27cbd80b577d296b46519aca02a22ad4ef997a7fe2 equals code_sha256 in both full/<cohort>/COMPLETE.json, and equals the copy in this folder.

## Outputs and inputs hashed

| path | bytes | sha256 |
|---|---|---|
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/manifest.json | 2658135 | 02b020362f3d65708d5a187effda21e9e83973fd20e07d4420ad6df9152f1249 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/CAPSCAN/accumulation.tsv | 33672 | 377b1c2f36023d9343cebf9a0c321c125fedbdc0a2a09111f23ff7640f8dc3d7 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/CAPSCAN/COMPLETE.json | 592 | 15446469f952e3e244ff574025323010f00230d61a7ea46fcb6b6d83872127b6 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/CAPSCAN/FL_public_sequence_incidence.npz | 13591181 | 3562a59f2a648123f0c1e2b34129b7bce373c49321927c4b72db34905e7fa50d |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/CAPSCAN/INPUT_RECEIPTS.json | 563018 | 19407699be6f5b56d258b7b1df622b1a50100ddce5fdccb6706248cbb84bede2 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/CAPSCAN/MG_sequence_incidence.npz | 899584820 | ae57177d8add9f6ea3c428e566d1dedea48b84d59f8e00524f0123c6480b6b96 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/CAPSCAN/MG_SEQUENCE_PRESERVATION.json | 31145 | be8af33dd2cb2988d39ac4cbec826571c4d5ac4576872f9bff606d504e94d09c |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/CAPSCAN/per_acquisition.tsv | 49063 | 71fe0e9a2c2b7960ef3845544aeeaa9e38152147e4d36f18521527b8f1689ac6 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/CAPSCAN/per_specimen.tsv | 18168 | e4dd948261ce2e06a9095490204fdc44b153588cef3e25ef06909e2e548495f5 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/PREDICT/accumulation.tsv | 45832 | 659958991f13680c169206930c230cb510908194e281b2ce7fd622a432fa44c9 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/PREDICT/COMPLETE.json | 592 | a61bd6ccd309843cb0869bde68373ee237f05bf1b6416c430c44ffb6080f25ee |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/PREDICT/FL_public_sequence_incidence.npz | 17430042 | 7f9b9ebbe3db273daa8168bad7ed2fa53ea5355ecf91b359ab6328ceda96c918 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/PREDICT/INPUT_RECEIPTS.json | 616775 | 56806d5190ffa3320407be45e2c6a6919675b916513fc718f3f59689ca5e2467 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/PREDICT/MG_sequence_incidence.npz | 1407755987 | 71dc3cc64e38927f634b520400047242d0212bcd3c707716df708a0e46d79d46 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/PREDICT/MG_SEQUENCE_PRESERVATION.json | 35050 | eafacec77b8bcdaa7416d03fcc4bd8e89f9edded29b07450ae8e86b598e1a5a0 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/PREDICT/per_acquisition.tsv | 62913 | aebb39d48b9098e78637e6a55a4a33518ac899fb0c5a5909e3c166d45c97cef5 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/PREDICT/per_specimen.tsv | 29996 | 461f167509940e01b3602d87db40404c8717af9a42e88caea5b9f605c53cf8d2 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/figures/SUMMARY.json | 1219 | a1fb1845b83d07e848151cd9215081c71b0266760c1463901970a5a3feb70251 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/figures/COMPLETE.json | 4445 | c58370463dc76aa0a69d4d81ec7871bb245d7c9f2f6b6cfffeee28d382978157 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/figures/sequence_accumulation_dual_axis_linear.pdf | 42078 | 43bc7660e08c21e96047501b00773a26812ff696b261ad0597e5ed686d96a68b |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/figures/sequence_accumulation_dual_axis_log.pdf | 42890 | a43a8f9c30ae516d024bee2c496db4441de9c04ea10a7370b658922dce422495 |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/figures/sequence_accumulation_order_averaged.pdf | 18784 | 0d9f5c402958fd43a582e9ceefd2b7b099d4399ff19acba6ac74c7117cc3716a |
| LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/worker.py | 8399 | 331ad8e43660a0e4a25d9d27cbd80b577d296b46519aca02a22ad4ef997a7fe2 |
