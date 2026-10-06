# Receipts

Original Slurm jobs, measured resources (from `sacct`, via `sacct -P -j`), and sha256 of every output file the figures and tables can draw on. Paths are relative to `LAB_ROOT`. Per-acquisition directories that carry a `COMPLETE.json` are summarised by the sha256 of the sorted list `<sha256 of COMPLETE.json>  <path>`; each `COMPLETE.json` itself records the sha256 of every output of its acquisition. Hashing was done on Slurm with `hashjob.py` (see the end of this file).

## Slurm jobs

| job id | stage | states | CPUs allocated | memory requested | wall time max (min) | peak RSS (GB) |
|---|---|---|---|---|---|---|
| 5929070 | pilot_verified_completed_inputs | COMPLETED: 1 | 2 | 32G | 0.1 | 0.0 |
| 5929071 | denovo | COMPLETED: 1 | 2 | 32G | 2.7 | 2.0 |
| 5929072 | capscan | COMPLETED: 1 | 1 | 32G | 6.8 | 2.1 |
| 5929073 | predict | CANCELLED: 1 | 0 | 32G | 0.0 | 0.0 |
| 6058101 | participant-15 run: denovo_threshold_routes_P15 | COMPLETED: 16 | 1 | 2G | 0.5 | 0.0 |
| 6058102 | participant-15 run: A27_available_P15 | COMPLETED: 16 | 1 | 6G | 2.9 | 3.6 |
| 6058103 | participant-15 run: terminal_route_counts_P15 | COMPLETED: 16 | 1 | 2G | 0.6 | 0.4 |
| 6058104 | participant-15 run: A27_unspecific_P15 | COMPLETED: 16 | 1 | 6G | 1.1 | 1.6 |
| 6058105 | participant-15 run: A24_baseline331 | FAILED: 1 | 1 | 3G | 0.1 | 0.0 |
| 6058106 | participant-15 run: A28_baseline331 | FAILED: 1 | 1 | 3G | 0.1 | 0.0 |
| 6058107 | participant-15 run: A24_fixed347 | CANCELLED: 1 | 0 | 3G | 0.0 | 0.0 |
| 6058108 | participant-15 run: A28_fixed347 | OUT_OF_MEMORY: 1 | 1 | 3G | 281.6 | 2.6 |
| 6058109 | participant-15 run: numbers_and_RESULTS_fig6 | CANCELLED: 1 | 0 | 8G | 0.0 | 0.0 |
| 6058147 | participant-15 run: A24_baseline331_resubmit(6058105_failed_missing_parent_dir) | COMPLETED: 1 | 1 | 3G | 6.8 | 2.2 |
| 6058148 | participant-15 run: A28_baseline331_resubmit(6058106_failed_missing_parent_dir) | OUT_OF_MEMORY: 1 | 1 | 3G | 308.1 | 2.6 |
| 6058149 | participant-15 run: A24_fixed347_resubmit(6058107_cancelled) | COMPLETED: 1 | 1 | 3G | 7.3 | 2.2 |
| 6058150 | participant-15 run: numbers_and_RESULTS_fig6_resubmit(6058109_cancelled) | CANCELLED: 1 | 0 | 8G | 0.0 | 0.0 |
| 6064030 | participant-15 run: A28_baseline331_resubmit2(6G) | COMPLETED: 1 | 2 | 6G | 149.9 | 2.9 |
| 6058023 | participant-15 run: extra | COMPLETED: 16 | 16 | 405G | 68.0 | 269.7 |
| 6064031 | participant-15 run: extra | COMPLETED: 1 | 1 | 6G | 271.3 | 3.0 |

## Output hashes

### `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results`

| sha256 | file |
|---|---|
| `b145d6eba201d157a32cfb474dc3afca40083996d181857372932d543114a1c2` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/CAPSCAN/COMPLETE.json` |
| `aec06b03b5a3703740ec80575fbb4d0d1118b4c03d9b29665a6a3473d803b1da` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/CAPSCAN/REPORT.md` |
| `cea2a987563c7388f033810471fa82b59ce7591ca35031cefaad09709b54edf3` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/CAPSCAN/acquisition_design.tsv` |
| `05ab47153209a4fc070b988031ec03ddd00f0d8f207d46901304171a85e964a5` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/CAPSCAN/complete_statistical_families.tsv` |
| `893647f9402f6b1af347471a88b8a936c1525db04869590e35941bdbf8df6de5` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/CAPSCAN/participant_paired_differences.tsv.gz` |
| `15b0d3363e06c451e05bce99b2946c946ac71a2e986854d14599d124f86d2f1b` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/CAPSCAN/per_acquisition_endpoints.tsv.gz` |
| `a87c470f38fea6f313b8c981aba61b1853d467bd661a1d4f210cb6923f6086e7` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/CAPSCAN/per_participant_region_endpoints.tsv.gz` |
| `eb21f3f041dbbea1dba0a770bb3d9cbddb4b34ed5ee96755dd8b833361f1ea28` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/CAPSCAN/per_specimen_endpoints.tsv.gz` |
| `024434e7747f944316d145e250dadc180d7237fa26626b0910fc224dbafdf46f` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/CAPSCAN/region_summary.tsv` |
| `c2675db4ed4b70fda2272f4565891fc4b05bf99bb8d3878ad950d6ba89d6133a` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/CAPSCAN/statistical_family_summary.tsv` |
| `1fb651036b219b2183daab05e7c43c7495a9f9853b451ea9b23065c78f17d5a7` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/CAPSCAN/threshold_pairs_historical_all_lengths.tsv.gz` |
| `891a267c9dfbb457746b550cddde3cff8b2929f0e0416bc6ef89d01130cbc4ee` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/CAPSCAN/threshold_pairs_same_9_50_lengths.tsv.gz` |
| `9eccf82c3812b7f00467577f1a3b1b3bc60e28fb845a7b6fbdef71f4d288b80d` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/CAPSCAN/threshold_pairs_same_production_eligibility.tsv.gz` |
| `f176eeb98c63aa01361f142a2e3a74ac99b1cc476222bd1f068bc1bfbae40f20` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/denovo/COMPLETE.json` |
| `aec06b03b5a3703740ec80575fbb4d0d1118b4c03d9b29665a6a3473d803b1da` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/denovo/REPORT.md` |
| `6b46fe011fb6f540f4740d901d00c8422b33a5a908a1f874221a659dcb49df31` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/denovo/acquisition_design.tsv` |
| `20904c728944a04aea0f949790e92e805fa940efe9dbcfe659eeee361604b7eb` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/denovo/complete_statistical_families.tsv` |
| `483d0791cc58834432a0710551e7dfafcf8be2862624e381a7cb577742c917ce` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/denovo/participant_paired_differences.tsv.gz` |
| `9b6df56d25065e0e8520411b66f7587b75db114841ba8808466fb8eace76e7d9` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/denovo/per_acquisition_endpoints.tsv.gz` |
| `627db02f807a3b0863326e03ab03a92f4d4814ebb6c6ad89891aad6fe808ed0f` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/denovo/per_participant_region_endpoints.tsv.gz` |
| `eedfc9c3cca927cac669965bbd3184c1fb39c5526c982e491e289cb1b14c40d2` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/denovo/per_specimen_endpoints.tsv.gz` |
| `9ef8278423ea5cf19cfbe09742d1b316737a8c97e98a4f22ba17dfb8df6bdc1e` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/denovo/region_summary.tsv` |
| `d7f622b21bb66245e8cd928e239fa366889012b8fd85a15142f19dc6c1307789` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/denovo/statistical_family_summary.tsv` |
| `94e679d8530a602cd4746471f94d28abf5ba2642db75437c53431ae96f99f156` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/denovo/threshold_pairs_historical_all_lengths.tsv.gz` |
| `d7b7dafb2a6724c4e5b5dc47db254e316e37901f08d89829825160a9b8981633` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/denovo/threshold_pairs_same_9_50_lengths.tsv.gz` |
| `de173a0dffdbdc2f4ef7686844d2e4bc803edd756d687de9ec67254a1c0f6b73` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/denovo/threshold_pairs_same_production_eligibility.tsv.gz` |
| `37f4b21012e5ffa663efe28342413109c0d00bda03c34f337d8ef29b606dc024` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/pilot/COMPLETE.json` |
| `aec06b03b5a3703740ec80575fbb4d0d1118b4c03d9b29665a6a3473d803b1da` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/pilot/REPORT.md` |
| `40d1aa56d970b21de06dfe68792998f8d8e66be183d820a49f100b8f77c08cf4` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/pilot/acquisition_design.tsv` |
| `aff33a45d50cd46f2106a5366cddf4e2695d254a7a3178631aa73dc30cd52fbe` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/pilot/per_acquisition_endpoints.tsv.gz` |
| `e15ccc6cb7514525db6b8dd74424d084fa92e2150c66bf034538b9c6bd62bc17` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/pilot/per_participant_region_endpoints.tsv.gz` |
| `640bbe6fd6fefaf03ab17ba04f9d809d94d3a360499e214d96aea780e8263d4f` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/pilot/per_specimen_endpoints.tsv.gz` |
| `5e183d3d33a998178392f955602c3146802d3acbcee94b76595d9951e958df93` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/results/pilot/region_summary.tsv` |

### `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/inputs/manifest.json`

| sha256 | file |
|---|---|
| `538ccef398311dec9113e0b844e8529003c09b3c8c2b4a95084dc77c4c27022b` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/inputs/manifest.json` |

### `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results_revision4`

Directory not present at the time of hashing.

### `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results_revision2`

Directory not present at the time of hashing.

### `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results`

| sha256 | file |
|---|---|
| `5a867ed25cb29275cca24fb0ad6ff1fa2910a4ee6385728dab13ac9f3b431999` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results/PREDICT/COMPLETE.json` |
| `aec06b03b5a3703740ec80575fbb4d0d1118b4c03d9b29665a6a3473d803b1da` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results/PREDICT/REPORT.md` |
| `0c47ea4db0fc5a9ee06b3161b264efd4126ba0f93fb289bc5ae9948c4976c796` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results/PREDICT/acquisition_design.tsv` |
| `dff7051d29fec43dbc250a35acf9a65588477a1bec0ed314be2f085c5ca9e1c3` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results/PREDICT/complete_statistical_families.tsv` |
| `d351e5716f9e991134cbc481297cc771681bee6cccb76b06a8190ab82e8d55b4` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results/PREDICT/participant_paired_differences.tsv.gz` |
| `0a2111fc504ca95404302bf4297468f7493a9b1a4dc79ab69f5ffbdbeda4fe41` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results/PREDICT/per_acquisition_endpoints.tsv.gz` |
| `e5e3b59e8ecddf8b864cc941d1a73fefa7ef9ed32bc73dbf955a6236d96ec1d9` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results/PREDICT/per_participant_region_endpoints.tsv.gz` |
| `de29c32becb0ea2dd28b492d232db0a3430a70b7476a80d57ff463fddaca0dc5` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results/PREDICT/per_specimen_endpoints.tsv.gz` |
| `c0d3d352eb92d2a8a22767492149f2615ba50768b5c8b78084293f2b0fc9efe2` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results/PREDICT/region_summary.tsv` |
| `629d2a5c819253736002a205698fe1ba55170af0b177fec4b270eee32c6d9aed` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results/PREDICT/statistical_family_summary.tsv` |
| `c171d7123550402be95864aa1800e0bf34d61f1ab023703798d8250e5c30b5c2` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results/PREDICT/threshold_pairs_historical_all_lengths.tsv.gz` |
| `719705146d06c5770ae7d20e16d758a42777de800927494f7dcd736daf9530a7` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results/PREDICT/threshold_pairs_same_9_50_lengths.tsv.gz` |
| `6d3682ae358b75b7d761daa4dec983ae1d7f5bf9a8d70dcd0b5c09c038fbadd7` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A24/results/PREDICT/threshold_pairs_same_production_eligibility.tsv.gz` |

### `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_baseline331/results`

| sha256 | file |
|---|---|
| `84406709236362c9da74ef4fb52f6885c9d8e162831dcc9ec418a5e09297e37c` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_baseline331/results/CAPSCAN/COMPLETE.json` |
| `aec06b03b5a3703740ec80575fbb4d0d1118b4c03d9b29665a6a3473d803b1da` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_baseline331/results/CAPSCAN/REPORT.md` |
| `cea2a987563c7388f033810471fa82b59ce7591ca35031cefaad09709b54edf3` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_baseline331/results/CAPSCAN/acquisition_design.tsv` |
| `05ab47153209a4fc070b988031ec03ddd00f0d8f207d46901304171a85e964a5` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_baseline331/results/CAPSCAN/complete_statistical_families.tsv` |
| `03c2de8d64e86010fe7d68ff6cb8cf90fd2c19c6d189d0fd25e34cfd6c64f6f6` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_baseline331/results/CAPSCAN/participant_paired_differences.tsv.gz` |
| `fb1dab0b8e2763d04065a4e8613abf97e1bcb1bc0f9d3b49946115d23403a3f5` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_baseline331/results/CAPSCAN/per_acquisition_endpoints.tsv.gz` |
| `bb6a48bffd0aa232b8470a0949ceae1bd6e8cb92929c1444efa1d02cd13633eb` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_baseline331/results/CAPSCAN/per_participant_region_endpoints.tsv.gz` |
| `1544e7eac9580ea73746be2e1ff8a03d81524ad5c05a0847cdfc1d14b94fbce0` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_baseline331/results/CAPSCAN/per_specimen_endpoints.tsv.gz` |
| `024434e7747f944316d145e250dadc180d7237fa26626b0910fc224dbafdf46f` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_baseline331/results/CAPSCAN/region_summary.tsv` |
| `c2675db4ed4b70fda2272f4565891fc4b05bf99bb8d3878ad950d6ba89d6133a` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_baseline331/results/CAPSCAN/statistical_family_summary.tsv` |
| `21d1948e104ea8c432eb94b12371215398fb5c1e4b67221b1673f7bee870247e` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_baseline331/results/CAPSCAN/threshold_pairs_historical_all_lengths.tsv.gz` |
| `505f02b4c00c761d849394164e61b0dff1139c813fd6ef5b8de3608540a62726` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_baseline331/results/CAPSCAN/threshold_pairs_same_9_50_lengths.tsv.gz` |
| `743e1fa00ed64fd665a433d5755fcc4fc01a061d8e2789a5151a6c84dbcbd8c6` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_baseline331/results/CAPSCAN/threshold_pairs_same_production_eligibility.tsv.gz` |

### `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_fixed347/results`

| sha256 | file |
|---|---|
| `efe314443e76a9a87224c9a79931146044830a6b2a952b7d6b84ad73e5b27951` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_fixed347/results/CAPSCAN/COMPLETE.json` |
| `aec06b03b5a3703740ec80575fbb4d0d1118b4c03d9b29665a6a3473d803b1da` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_fixed347/results/CAPSCAN/REPORT.md` |
| `4638d7894cdd3423e33a7416bfe2de11d9e761e6bbde4749d29e0ed4a03c7fc4` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_fixed347/results/CAPSCAN/acquisition_design.tsv` |
| `919360392724871c2cef8e99e5a0172cdb7d8ad8a364a476778aac7c8f10f70b` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_fixed347/results/CAPSCAN/complete_statistical_families.tsv` |
| `d0bcfb79b1c10ca3d719affc152bdfe8ffba1809067de9e8d358d47b637e9b90` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_fixed347/results/CAPSCAN/participant_paired_differences.tsv.gz` |
| `b214eb69ce8d243b201403f3599a65b1fcfd31ee767d9c4b8edca4d7018d2e65` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_fixed347/results/CAPSCAN/per_acquisition_endpoints.tsv.gz` |
| `6d998e903173e85cff0166be599d05fc0005d1b462bf29101a8a93fd3f746dac` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_fixed347/results/CAPSCAN/per_participant_region_endpoints.tsv.gz` |
| `e8a7c91163d7a927148fbf2f3825fff216503b79bed0f621d626639931c4cb6b` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_fixed347/results/CAPSCAN/per_specimen_endpoints.tsv.gz` |
| `6fe4adce9406b3b3f8673a3a663f2f665ec8c1291159ccfe1e9f8128c68d501a` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_fixed347/results/CAPSCAN/region_summary.tsv` |
| `e6b03cc3b7ed63ec663354d0e268c6e361c9fe5481aa290272f82a8ea6e16d52` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_fixed347/results/CAPSCAN/statistical_family_summary.tsv` |
| `3d57339c7e81580cf25ccabd8fabc9b2fb38acb18817d4110d36392426bdb738` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_fixed347/results/CAPSCAN/threshold_pairs_historical_all_lengths.tsv.gz` |
| `a4d6486289f6dac7f06588a3e2f2ace1b11b5ac137f93122178fd7f7c4461f1c` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_fixed347/results/CAPSCAN/threshold_pairs_same_9_50_lengths.tsv.gz` |
| `b2690873b4427aaccab0655407a018246ab6c1356bc46e1c1b0005e47a1bb9df` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A24_fixed347/results/CAPSCAN/threshold_pairs_same_production_eligibility.tsv.gz` |

## Released code versus origin

| file | origin (LAB_ROOT relative) | sha256 at origin | sha256 here |
|---|---|---|---|
| `search.py` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/frozen/search.py` | `c980eccbf1bc2989a7ea1ac4e3b7eb3baa3446eb33d5f3e869b1a54a19f34a8e` | `c980eccbf1bc2989a7ea1ac4e3b7eb3baa3446eb33d5f3e869b1a54a19f34a8e` |

Static checks: `python -m py_compile` on every `.py` file; `tools/scan_release.py`.

## Hashing code

The hashes above were produced by a Slurm job (4 to 8 CPUs) that walks each root once, hashes every file with SHA-256 in 8 MiB blocks and, for roots with more than 20 per-acquisition directories, hashes only the `COMPLETE.json` of each. The code is kept in the work directory, not in the repository.
