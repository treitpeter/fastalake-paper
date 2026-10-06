# Receipts

Original Slurm jobs, measured resources (from `sacct`, via `sacct -P -j`), and sha256 of every output file the figures and tables can draw on. Paths are relative to `LAB_ROOT`. Per-acquisition directories that carry a `COMPLETE.json` are summarised by the sha256 of the sorted list `<sha256 of COMPLETE.json>  <path>`; each `COMPLETE.json` itself records the sha256 of every output of its acquisition. Hashing was done on Slurm with `hashjob.py` (see the end of this file).

## Slurm jobs

| job id | stage | states | CPUs allocated | memory requested | wall time max (min) | peak RSS (GB) |
|---|---|---|---|---|---|---|
| 5929609 | pilot | COMPLETED: 1 | 2 | 32G | 1.2 | 0.7 |
| 5929610 | available_capscan | CANCELLED: 1 | 0 | 32G | 0.0 | 0.0 |
| 5929611 | available_predict | CANCELLED: 1 | 0 | 32G | 0.0 | 0.0 |
| 5929612 | full_capscan | CANCELLED: 1 | 0 | 32G | 0.0 | 0.0 |
| 5929613 | full_predict | CANCELLED: 1 | 0 | 32G | 0.0 | 0.0 |
| 5930574 | pilot_revision2 | COMPLETED: 1 | 2 | 32G | 1.2 | 0.6 |
| 5930575 | statistics_smoke | COMPLETED: 1 | 2 | 32G | 12.2 | 0.9 |
| 5930576 | available_capscan_revision2 | CANCELLED: 1 | 0 | 32G | 0.0 | 0.0 |
| 5930577 | available_predict_revision2 | CANCELLED: 1 | 0 | 32G | 0.0 | 0.0 |
| 5930578 | full_capscan_revision2 | CANCELLED: 1 | 0 | 32G | 0.0 | 0.0 |
| 5930579 | full_predict_revision2 | CANCELLED: 1 | 0 | 32G | 0.0 | 0.0 |
| 5931248 | pilot_revision3 | COMPLETED: 1 | 2 | 32G | 1.3 | 0.7 |
| 5931249 | statistics_balanced_revision3 | FAILED: 1 | 2 | 32G | 3.0 | 0.7 |
| 5931250 | available_capscan_revision3 | CANCELLED: 1 | 0 | 32G | 0.0 | 0.0 |
| 5931251 | available_predict_revision3 | CANCELLED: 1 | 0 | 32G | 0.0 | 0.0 |
| 5931252 | full_capscan_revision3 | CANCELLED: 1 | 0 | 32G | 0.0 | 0.0 |
| 5931253 | full_predict_revision3 | CANCELLED: 1 | 0 | 32G | 0.0 | 0.0 |
| 5931303 | independent_statistics_audit | CANCELLED: 1 | 0 | 8G | 0.0 | 0.0 |
| 5931567 | pilot_revision4 | COMPLETED: 1 | 2 | 32G | 1.3 | 0.5 |
| 5931568 | statistics_balanced_revision4 | COMPLETED: 1 | 2 | 32G | 12.1 | 0.8 |
| 5931569 | independent_statistics_audit_revision4 | COMPLETED: 1 | 2 | 8G | 0.2 | 0.0 |
| 5931570 | available_capscan_revision4 | COMPLETED: 1 | 1 | 32G | 223.7 | 2.2 |
| 5931571 | available_predict_revision4 | COMPLETED: 1 | 1 | 32G | 154.4 | 1.6 |
| 5931572 | full_capscan_revision4 | COMPLETED: 1 | 1 | 32G | 282.3 | 2.9 |
| 5931573 | full_predict_revision4 | CANCELLED: 1 | 0 | 32G | 0.0 | 0.0 |
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

### `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4`

| sha256 | file |
|---|---|
| `3304b7e1ce458800f8c9ad46f75bb32d6a50b2c928bef3a0d0e2e46f81247f8f` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/COMPLETE.json` |
| `fad2737bf7a636206d45f5e69f5e3fad419053043adb58d8a40e804ffed00fcd` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/P1_classification_coverage.tsv` |
| `046e0053b5ee11847c5d0d01413a81472aeb27522d46d3b6031551b5abaa3c2a` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/REPORT.md` |
| `f9e89b636179fa5d9b9942e3bd62b1485bb974a082ea01229fdf4f32a7a1475b` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/acquisition_enrichment.tsv.gz` |
| `027deb4ac4caa468a2baeebb74f30cf882b0c75b56369c1f5ad5ce68f53dcd7a` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/acquisition_motif_counts.tsv.gz` |
| `04768624f33db3e147640f69dd8d40223604c597c0b7d0de9c4f32aff311be0d` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/boundary_recurrence.sqlite` |
| `aba18c61ea96acdc0461eb3cdc6aa09e4b014ec42331972fc53951c42e9436a6` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/boundary_recurrence.tsv` |
| `f9ae06df1ad175b29537e9cf20e48e78fd21c29553cdce1d3da394440e04bccd` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/complete_statistical_families.tsv` |
| `21a0a0182e69f19a4251e6d54b3f1dc621b217584226feaecc59d947924346a9` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/mapping_coverage.tsv` |
| `1e9db01458c077ca7ac9aec65310864c6560f96ffbe1c7ce84eefb293b9dca36` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/matched_reciprocal_controls.tsv` |
| `cab4ecd597f3e7ed89c390bb5548bf999f6e18d7ebb6a0da9abc1eefde060690` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/method_overlap.tsv.gz` |
| `8a4e58b769bccfce288c529ff9b7340b0f4cbfb60ffa6bce6d98f21a086cb5a0` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/participant_concentration_sensitivity.tsv` |
| `a1d772868d553839d98c3a375c3beaf058c1bc840f9cec85dbc7633b2ab1e313` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/participant_motif_differences.tsv.gz` |
| `12db148fec005924f2b421076293e4acd3b1436c9acfc74e44de5183bad9c90e` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/participant_reciprocal_controls.tsv` |
| `f8f8ea597dbf43b75b263d78aa46677b7baab2bc0c1b258745f717f683bd97fe` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/participant_region_enrichment.tsv` |
| `7b93aa8d624e65c6f5ae932ac64d1d7c514e00cb7e121941ebdedabfdfb0a2b8` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/protein_concentration_sensitivities.tsv.gz` |
| `dcfd0a0a457e3e56d9422d0e8d4797cbe15b816763f6b1d3b0cffed6e0e42537` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/statistical_family_summary.tsv` |
| `8fe78376a874cc211b320e0ae190c8549ea1f0266b7268ac157b0d3f0292ccd3` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_CAPSCAN/unresolved_P1_sites.tsv.gz` |
| `7dd89d355154995d0bc301b5828aa1430a3b5c451d2b89f767f1f021d003e493` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/COMPLETE.json` |
| `a4b1f3d6434ed608664d20eb1a770056573549b79991bc8735137e3fe13cabb8` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/P1_classification_coverage.tsv` |
| `046e0053b5ee11847c5d0d01413a81472aeb27522d46d3b6031551b5abaa3c2a` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/REPORT.md` |
| `ce0db59e5dbc0e049e5b8d8352bb11341b9a144c9b92816e302a71d1200649f6` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/acquisition_enrichment.tsv.gz` |
| `bb1954d360cbae3d98beef93850c0e69acad8e8a9e6a9c3c46be132223248c3e` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/acquisition_motif_counts.tsv.gz` |
| `0d6dfa843671ff142cb0daaf56f7de8e4796159a59984fe5b17169c3c93b2d6a` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/boundary_recurrence.sqlite` |
| `fb746d840712dd744bf8af12db7b2c20693c3f0d83443b9b53f7c21b783cc634` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/boundary_recurrence.tsv` |
| `c552413d3f2719a646facac9c7bf06bd4425dde5dd5274222d0ad6dfed3e72cc` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/complete_statistical_families.tsv` |
| `36bf77ad94abe4a7caf54e8a9dc0643b425410ac35c1ace0880468941986c57e` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/mapping_coverage.tsv` |
| `a8e6b37aef0406919ec48d43c86781f77520aac5ba5988885699f09f6ea3dcd7` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/matched_reciprocal_controls.tsv` |
| `4f1c4759b803fea5b2fb846dd499ff89386a410f1dc7892c8427bd532357408b` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/method_overlap.tsv.gz` |
| `6ded39f2391a995916b5629ae4b22b1dfd4df555e645acdf64856b4b50735a7a` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/participant_concentration_sensitivity.tsv` |
| `ebdf1e54f944e7c0b0870cbf9cde82891bd6b0a2df83abc58fd01d2fc8bfc3f9` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/participant_reciprocal_controls.tsv` |
| `ef0236bdfa30908293b5f189bb6b98dee490a2c597922ee371826b39fcfe713c` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/participant_region_enrichment.tsv` |
| `60ecd61b453ed5d949318539777761a59e4703dc5648a81a74f58793fee5fc5b` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/protein_concentration_sensitivities.tsv.gz` |
| `72bec8d707f2a2e7b74188577c73116041f66b6d95154ae549e6c0e138f3f614` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/statistical_family_summary.tsv` |
| `476813ada32fca7e6ceda4fc964a2977f9b5e986606c1c8f45c1665e811d6a97` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/available_PREDICT/unresolved_P1_sites.tsv.gz` |
| `d5a466f5987cdf3e4d86e7d38f0a902733574cc722db762693a2c652eec7d765` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/COMPLETE.json` |
| `f3b5a5858c4e2ff1463019c02715d6bb4a810ff1fc379330bbfcf2a0bb7f8713` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/P1_classification_coverage.tsv` |
| `046e0053b5ee11847c5d0d01413a81472aeb27522d46d3b6031551b5abaa3c2a` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/REPORT.md` |
| `67d1d2bc68a20efba8b3d3ce5d3e405d3000dd1dc10593e55763281ea4e2d446` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/acquisition_enrichment.tsv.gz` |
| `7da4fc18a2b8470da19557b04d31252a5f8ad6615f47b88092134ee80290868a` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/acquisition_motif_counts.tsv.gz` |
| `4ef8ada0d9685ae46dd2f0be45d307d7eb04a2dbdfc5f2622736bee7a431ca33` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/boundary_recurrence.sqlite` |
| `7239ab71ec981e4048cc8760a37bc21dee4ff3df060f186b20cabc1dd25ffbd8` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/boundary_recurrence.tsv` |
| `90debaad4141a8b9ddece9fed1029612a103bd86e89c0bc1f67e3ee288c91a2f` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/complete_statistical_families.tsv` |
| `ec6caca80682bd34dc6ed357b7c4fda845101a8bd3bbeb3de7abd6f571e756c6` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/mapping_coverage.tsv` |
| `9f3191dc1a667cdb14197df5645e5a09104150ee7e27452e96c4ab6ca55d2c18` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/matched_reciprocal_controls.tsv` |
| `fbffda7a5134c94fe9f66bbb3f0e2f0cd12e2434bba91c14e967b5a313844ee6` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/method_overlap.tsv.gz` |
| `c43bf5bd8ac40ff031f4e6f0f6ae49028205786c0fb745638f02823fdb05787e` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/participant_concentration_sensitivity.tsv` |
| `0d89678b7650c7b673054c6a44618122a24e3f2e063ce33a251fa748916fd182` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/participant_motif_differences.tsv.gz` |
| `5c864175c2ab5c31f0d9f1913726475e44d22a3a1f612b893039662bdbb8dc6f` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/participant_reciprocal_controls.tsv` |
| `9d42a56a7d7dd9581b590c9b076b09317a69de215a1794d7bcaed65212e0665b` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/participant_region_enrichment.tsv` |
| `891bd0de10b8499d2f5ef433d07741338d22bdb787ecc229d8591582cb999c07` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/protein_concentration_sensitivities.tsv.gz` |
| `ac923bad3535f87c0b7830ba17c09f881aa7d5c77ede950b0cfcda902f297da0` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/statistical_family_summary.tsv` |
| `6eeefabf5130df4eeb1d6f73c5ceca38f0dd77b8380b512c43474684158f0c5f` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/full_CAPSCAN/unresolved_P1_sites.tsv.gz` |
| `b232cb0e9a64ad2ae82ad7226751ee3dbdcbde1f1fafd084a642f428b9343514` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/independent_audit/AUDIT.json` |
| `21fe37d74a13f1f6954e8de0fbe55bc6d6d933abea39ec48f1be190d895830f4` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/independent_audit/COMPLETE.json` |
| `528283f1baa29d65bcfc47d153c37df0a228ced0a1cbfd9c937554f4537ac7da` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/COMPLETE.json` |
| `5d20e5e1815dac7b9b37d1251f8e32bd19f9ffc9285ea774475cbee9d12f2a0c` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/P1_classification_coverage.tsv` |
| `046e0053b5ee11847c5d0d01413a81472aeb27522d46d3b6031551b5abaa3c2a` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/REPORT.md` |
| `77416d48941111b17d7644611f7be0233569793483f6fb8b067858d5611ba895` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/acquisition_enrichment.tsv.gz` |
| `7b0153d45428e59292d039b6129f07180e5983468e92fa2bd7033aa5eb1efbb4` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/acquisition_motif_counts.tsv.gz` |
| `afe0e8e9a5e6f4bab8defb2e6ad9f45f4fe2480ed9f6f797490d57bf5098b8b5` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/boundary_recurrence.sqlite` |
| `fb773690007d4dac3d31bf20cad87721fa274bfd00f54530a5200d71cbef1df5` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/boundary_recurrence.tsv` |
| `9eb7e08534487186249431b7ffb2e4fdff8710bed5dba776106869373190cb76` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/mapping_coverage.tsv` |
| `804d7ab2d04f850a7d746b251493f9077d4ced2b83f3c42e392dc2f0469edd34` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/matched_reciprocal_controls.tsv` |
| `6c6d5fb7d83589f7c7b79dd21a8630bb4c6174796b139a40a2c90f784e5b9c9e` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/method_overlap.tsv.gz` |
| `1f4f6d88ba23c4a7d8da5bb47995a633e927f40cc750371378918618889dbeb4` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/participant_concentration_sensitivity.tsv` |
| `94b54c5722587aaa4ab37011c34d20227e2d53829ea7d149c81f4786e345c015` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/participant_motif_differences.tsv.gz` |
| `b046bd3ad833062faaabb27a98ca6e84d1e7ba5c9fdbb6181c99c9cb1f89f6cb` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/participant_reciprocal_controls.tsv` |
| `4be35ada0a4f96359c251aab5e0fb030a55fe538a1151be0ffac8252df29074d` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/participant_region_enrichment.tsv` |
| `b93ae1e73e5559ee72120e6eb1ab497f1efa55d66ff9cc0a4d613cef06df0a46` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/protein_concentration_sensitivities.tsv.gz` |
| `dfff93c607f0c67c9e5877c0e3f9e1ed321bc231757e95795127ad14a0b03c69` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/pilot/unresolved_P1_sites.tsv.gz` |
| `7ef122e41380f9495bd4f8ea98281e5eb670d4007b5914f611715ea8c6ba7d45` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/COMPLETE.json` |
| `0925391ea68f014a7ed6b6f96d80e4c59d3ba63c0d095ce2f1cb916f7c6f2c0e` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/P1_classification_coverage.tsv` |
| `046e0053b5ee11847c5d0d01413a81472aeb27522d46d3b6031551b5abaa3c2a` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/REPORT.md` |
| `f15b3d384a7c1c4eea84b10f44ec135674d4ebcd4cb8f26aa23e6296374965dd` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/acquisition_enrichment.tsv.gz` |
| `e8e14e6586213466c44679688a5cbbdd5c1cf4e245120869559dfe8ad5bf4522` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/acquisition_motif_counts.tsv.gz` |
| `40dcc50d723742f94f69f87fce3567d6f133fc3de760c7b72f631fc7e7633d68` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/boundary_recurrence.sqlite` |
| `49670bf95f03303abfe5033203464bb8e235ec0f8a28758a4baf447032acf692` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/boundary_recurrence.tsv` |
| `32eed5a87df604fb07089f60bcc637d216d1f786f94d1608569fff41ff2ac0c3` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/complete_statistical_families.tsv` |
| `27abb501ba652b7789eddd0111227f2965efafc70fcafb8d842eb0124715e094` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/mapping_coverage.tsv` |
| `884d475cebfc47121cdfa8b99f2fb7b619933bed83ad264104570d4b19c9cac1` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/matched_reciprocal_controls.tsv` |
| `cb12ba1f76399bcab25430082d47d39e477b97221e240c1fb9055648164a0406` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/method_overlap.tsv.gz` |
| `c3d00de1376b970790ea865d8d890b630e198510db614d7712157773e515f9de` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/participant_concentration_sensitivity.tsv` |
| `01c46d2515bfe428a59660c38fe9e8fc739c3a90d1eec9308550035922dbff6d` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/participant_motif_differences.tsv.gz` |
| `0b03a0eabb3605eb8f2bf10170bd9f04a694778af715e7f09d3870984ec8b88d` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/participant_reciprocal_controls.tsv` |
| `e4cc3938d3f867d4ef735c555993d48de5718ae89a60523f8e448fbbeb835b54` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/participant_region_enrichment.tsv` |
| `94ec0391e852f1e947916df091e4c0de4d53fde16c55a987397c14e5d71b0c83` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/protein_concentration_sensitivities.tsv.gz` |
| `ece5ad32b62b39013c0fd5173551e72e698e5b51138009979f1661c90d8e2380` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/statistical_family_summary.tsv` |
| `6d10d2797732d10f08abeef404b4172fb04c5103909b298fde7c055d7ea66e69` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/results_revision4/statistics_smoke_CAPSCAN/unresolved_P1_sites.tsv.gz` |

### `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/inputs/manifest_revision4.json`

| sha256 | file |
|---|---|
| `a3002790cda61894a107e36f1432442df0df043b16ba65210662d145e9d836d2` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/inputs/manifest_revision4.json` |

### `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4`

| sha256 | file |
|---|---|
| `789da31f6b0d6d93012e196885e89e4758c50334a55b3b40f4d2ca01b9840104` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/COMPLETE.json` |
| `9c377a94415cffd6b6092f3bf566ace390b4aeb56ebcca4ce7b557af033f4993` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/P1_classification_coverage.tsv` |
| `046e0053b5ee11847c5d0d01413a81472aeb27522d46d3b6031551b5abaa3c2a` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/REPORT.md` |
| `8e49c01bac8e2af1f8f2570c4718fc4bf7dea5d648af1ff62f6f5d88cc14c1ab` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/acquisition_enrichment.tsv.gz` |
| `93d8b0e71f723e43b9c9fc3b41e516d016406f1c0f9e00bb334b491369fccc96` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/acquisition_motif_counts.tsv.gz` |
| `5cb2afa6ed34e8b03fbd9586268597035f35be84c89d2ba6e72ad859f7906d1a` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/boundary_recurrence.sqlite` |
| `1780d6e56187522d2b5fbe5407b01b64e4299b3acfc2224aa8de74a3da30e34b` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/boundary_recurrence.tsv` |
| `85c774cb604a1d9e6c9526b04750b10d6601a07855fe1552357f66be016d28d7` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/complete_statistical_families.tsv` |
| `978600f7498862d10364e7084c58568288cb09c58dff91563e4cdd026490e556` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/mapping_coverage.tsv` |
| `0b335476fbafee3a80d1fef1d99e849242bf9e6053dbfa3f4f227df3df689300` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/matched_reciprocal_controls.tsv` |
| `084ffc028a481499939dd3cea792c7f77ad1102596e7ca478190b657aebf7539` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/method_overlap.tsv.gz` |
| `7cd6606785968e34746ef3ffd29b775d09751e0480bb32a564a78609b412ec05` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/participant_concentration_sensitivity.tsv` |
| `56dccbad6c2f028b812e28703075fb959e8eed768509fd207b3a664830ba81d9` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/participant_reciprocal_controls.tsv` |
| `b94e2cf16b4122ea9783572fc0a429eb2e08a2932b4adbd3fa9b93bfe226d36b` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/participant_region_enrichment.tsv` |
| `421db17287db70bf1861c0d04a3c0ee067e395a0ffdb427fccad0e7903776654` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/protein_concentration_sensitivities.tsv.gz` |
| `57ec8f66616f596162672e98b1ea9f970166f441287476bc3cb236ddaa7e3395` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/statistical_family_summary.tsv` |
| `832182a7dda3949d23d6c8a9eaed079c4fd28a0f0f56309be34ef334a7ec1e92` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/unresolved_P1_sites.tsv.gz` |

### `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4`

| sha256 | file |
|---|---|
| `615150bd72f92a693dd9769895e4bc92b74cf533f4c0ba93017cc38a300ee878` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/COMPLETE.json` |
| `f3b5a5858c4e2ff1463019c02715d6bb4a810ff1fc379330bbfcf2a0bb7f8713` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/P1_classification_coverage.tsv` |
| `046e0053b5ee11847c5d0d01413a81472aeb27522d46d3b6031551b5abaa3c2a` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/REPORT.md` |
| `8b225d6f3162fecebd70e75e6930fc94f46a5e58c08fa85cd230543ea8767c12` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/acquisition_enrichment.tsv.gz` |
| `09a9111f7d93eb610ce6366242aa552efc27477035bf8deecf2ffc11b650bef9` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/acquisition_motif_counts.tsv.gz` |
| `4ef8ada0d9685ae46dd2f0be45d307d7eb04a2dbdfc5f2622736bee7a431ca33` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/boundary_recurrence.sqlite` |
| `7239ab71ec981e4048cc8760a37bc21dee4ff3df060f186b20cabc1dd25ffbd8` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/boundary_recurrence.tsv` |
| `90debaad4141a8b9ddece9fed1029612a103bd86e89c0bc1f67e3ee288c91a2f` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/complete_statistical_families.tsv` |
| `ec6caca80682bd34dc6ed357b7c4fda845101a8bd3bbeb3de7abd6f571e756c6` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/mapping_coverage.tsv` |
| `9f3191dc1a667cdb14197df5645e5a09104150ee7e27452e96c4ab6ca55d2c18` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/matched_reciprocal_controls.tsv` |
| `182859460b925d7cb0204174b5af2290e81d55c4646337ff2ac99979cef491ed` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/method_overlap.tsv.gz` |
| `c43bf5bd8ac40ff031f4e6f0f6ae49028205786c0fb745638f02823fdb05787e` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/participant_concentration_sensitivity.tsv` |
| `c04e03b51d875179576ea5dbdfac5d6fa793b055dad1cd73edd3addd62bb233f` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/participant_motif_differences.tsv.gz` |
| `5c864175c2ab5c31f0d9f1913726475e44d22a3a1f612b893039662bdbb8dc6f` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/participant_reciprocal_controls.tsv` |
| `9d42a56a7d7dd9581b590c9b076b09317a69de215a1794d7bcaed65212e0665b` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/participant_region_enrichment.tsv` |
| `11eeef87fcccfafa0cd278cde5752ce2974f6788a46aee941e3d87d7d7856d8b` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/protein_concentration_sensitivities.tsv.gz` |
| `ac923bad3535f87c0b7830ba17c09f881aa7d5c77ede950b0cfcda902f297da0` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/statistical_family_summary.tsv` |
| `aecb73a67189e20b6d0536a1427ac7aeee109ce0d3977aaeaabcf0c823964592` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN/unresolved_P1_sites.tsv.gz` |
| `f3b5a5858c4e2ff1463019c02715d6bb4a810ff1fc379330bbfcf2a0bb7f8713` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN.oom_6058148/P1_classification_coverage.tsv` |
| `9e02d6a67555febcd2a989ec715cec2b9bb91e8b4ed1f37e18c606798a52390d` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN.oom_6058148/acquisition_enrichment.tsv.gz` |
| `63805f87054c8ebdec0a245c7fe9b090eb6622039d23d15a59f7c31bc5bfa0e8` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN.oom_6058148/acquisition_motif_counts.tsv.gz` |
| `4ef8ada0d9685ae46dd2f0be45d307d7eb04a2dbdfc5f2622736bee7a431ca33` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN.oom_6058148/boundary_recurrence.sqlite` |
| `7239ab71ec981e4048cc8760a37bc21dee4ff3df060f186b20cabc1dd25ffbd8` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN.oom_6058148/boundary_recurrence.tsv` |
| `ec6caca80682bd34dc6ed357b7c4fda845101a8bd3bbeb3de7abd6f571e756c6` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN.oom_6058148/mapping_coverage.tsv` |
| `9f3191dc1a667cdb14197df5645e5a09104150ee7e27452e96c4ab6ca55d2c18` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN.oom_6058148/matched_reciprocal_controls.tsv` |
| `30193517796cd4c9bc86fd479730c0e1dd50a1092e8343c6d51f7ed46e13076e` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN.oom_6058148/method_overlap.tsv.gz` |
| `9d42a56a7d7dd9581b590c9b076b09317a69de215a1794d7bcaed65212e0665b` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN.oom_6058148/participant_region_enrichment.tsv` |
| `bd7646cc7b6f6e78cb18100fba0546cef588ca870858d1f708dc0fdb2161d298` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN.oom_6058148/protein_concentration_sensitivities.tsv.gz` |
| `cb06c6812e0b360590d5ee3ddd82a909576fa0a834507c3e0171f676255c8783` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_baseline331/results_revision4/full_CAPSCAN.oom_6058148/unresolved_P1_sites.tsv.gz` |

### `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4`

| sha256 | file |
|---|---|
| `fc9bfa1ae55426e3acd59c66789b0aa760e1f8e89005182ef4301f57e0500299` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/COMPLETE.json` |
| `6c509617e709b782e7cfbb10c2d308584e973bb7cc4310abac32a57c3e6cc1f5` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/P1_classification_coverage.tsv` |
| `046e0053b5ee11847c5d0d01413a81472aeb27522d46d3b6031551b5abaa3c2a` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/REPORT.md` |
| `a7268d744c6246ffe204fc513535f62ef0504b2ce744ba8feb911b2966abc01f` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/acquisition_enrichment.tsv.gz` |
| `2f26f4704994b09608afece1294ba7c76ab93e5248febadfd1d221c4390e52e0` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/acquisition_motif_counts.tsv.gz` |
| `4bc9a510ef71d22ae1a5a05ee92b53b6ccfda006039047e319517480bc2af537` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/boundary_recurrence.sqlite` |
| `55c2ef6b0ebb0a085b4c850e093d33ef83d96835a94f58553ec3beab2c5137e5` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/boundary_recurrence.tsv` |
| `3c24a42eec725e0b43bdeecca3928f971925c25395dde7478d92a1db9e1c28d7` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/complete_statistical_families.tsv` |
| `47e9cd65c8a27d4c9a767e350c168fba7d2d09bde3b587feb5d5d26b26eea711` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/mapping_coverage.tsv` |
| `9df850531fdadc6616e4bee5b845061ad9e8e02f7b4745a0fc1ad30e7806825d` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/matched_reciprocal_controls.tsv` |
| `10492f64c013f2315d5bb6b9a133f10806bfa3635c61b5c917ee20003fba3455` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/method_overlap.tsv.gz` |
| `b5c9a6615f698896b707b33d64ea761f1c706de02a29a97a69d89dfc2f3f6779` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/participant_concentration_sensitivity.tsv` |
| `df84cc2720b0ed63448048ef0676ad33c44ef9fe46eb05414a1d98a3ecb5ae77` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/participant_motif_differences.tsv.gz` |
| `74613f6da70994b5307cf95660f40f1eb3a22cc9757175c48d5b28161980263c` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/participant_reciprocal_controls.tsv` |
| `a830d42712865d23f6642060c782b2dd00b49db19b551dc3575cdb359f619722` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/participant_region_enrichment.tsv` |
| `0905cceac71aeb1af9b6157ccfc232a0b04edd18f10b5ca27a6388a7a88bc3f1` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/protein_concentration_sensitivities.tsv.gz` |
| `1f707e19f73834a674ac26a7dafb85d3bd4d9697629ed075c8a3d7325ba100f4` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/statistical_family_summary.tsv` |
| `09b87f0626f3adef4cfa4ee020ca11210dc47cb1e97bada4745a6d492605827b` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN/unresolved_P1_sites.tsv.gz` |
| `6c509617e709b782e7cfbb10c2d308584e973bb7cc4310abac32a57c3e6cc1f5` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN.failed_6058108/P1_classification_coverage.tsv` |
| `9ddee796f07efdc5e2023020fc2bcce24b987e10d7886044186f0a9af1bd907f` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN.failed_6058108/acquisition_enrichment.tsv.gz` |
| `42a1685cd5d8b49e5b463aba4dbabef56f54918bb27983b7afbf25aa58fd615f` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN.failed_6058108/acquisition_motif_counts.tsv.gz` |
| `4bc9a510ef71d22ae1a5a05ee92b53b6ccfda006039047e319517480bc2af537` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN.failed_6058108/boundary_recurrence.sqlite` |
| `55c2ef6b0ebb0a085b4c850e093d33ef83d96835a94f58553ec3beab2c5137e5` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN.failed_6058108/boundary_recurrence.tsv` |
| `47e9cd65c8a27d4c9a767e350c168fba7d2d09bde3b587feb5d5d26b26eea711` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN.failed_6058108/mapping_coverage.tsv` |
| `9df850531fdadc6616e4bee5b845061ad9e8e02f7b4745a0fc1ad30e7806825d` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN.failed_6058108/matched_reciprocal_controls.tsv` |
| `5d757a74e01d9be0e5bec9a975be43d67aa33521252adf5b0f6a05da1d5b8093` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN.failed_6058108/method_overlap.tsv.gz` |
| `a830d42712865d23f6642060c782b2dd00b49db19b551dc3575cdb359f619722` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN.failed_6058108/participant_region_enrichment.tsv` |
| `543523317bb69e4606e29a890244182abb141531a75aba54606930e9b648d41c` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN.failed_6058108/protein_concentration_sensitivities.tsv.gz` |
| `c34d8ef7b2db1eba929a4a28e28c645fe87c741db4905f8918b676275b0f43c3` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/A28_fixed347/results_revision4/full_CAPSCAN.failed_6058108/unresolved_P1_sites.tsv.gz` |

## Released code versus origin

| file | origin (LAB_ROOT relative) | sha256 at origin | sha256 here |
|---|---|---|---|
| `collect_revision4.py` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/frozen/collect_revision4.py` | `95e49d899ba7e22144287f467ad7a70d1f42e02751fa49ed706c4eb841bdeb74` | `95e49d899ba7e22144287f467ad7a70d1f42e02751fa49ed706c4eb841bdeb74` |
| `audit_inference_revision4.py` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/frozen/audit_inference_revision4.py` | `4c064fd5031fbb5ccfb2b541490ee679c97315a3ac8cc2ebdb8e08a148321cea` | `4c064fd5031fbb5ccfb2b541490ee679c97315a3ac8cc2ebdb8e08a148321cea` |
| `histidine_helpers.py` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/frozen/histidine_helpers.py` | `387460e01d28ff90ed17037f1e5ac7d08acb33287b04f503cacb7d1971b888f7` | `387460e01d28ff90ed17037f1e5ac7d08acb33287b04f503cacb7d1971b888f7` |
| `review_histidine_predict.py` | `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/code/review_histidine_predict.py` | `22f10374ac480ca1f9fb7f3ff2f76d4c0d2ec9a919ea935ef2d908a05e2e0b9c` | `7bad9c14166ed92d54797cc9a26ebc2055b271e605b0c1692cea9e253809d680` |
| `collect_repaired_predict_20260915.py` | `LAB_ROOT/MM_review_analysis/code/collect_repaired_predict_20260915.py` | `96289d6e9629c3fd1b00bb5097341ca51601c857758ca383f244808776144e97` | `96289d6e9629c3fd1b00bb5097341ca51601c857758ca383f244808776144e97` |

Static checks: `python -m py_compile` on every `.py` file; `tools/scan_release.py`.

## Hashing code

The hashes above were produced by a Slurm job (4 to 8 CPUs) that walks each root once, hashes every file with SHA-256 in 8 MiB blocks and, for roots with more than 20 per-acquisition directories, hashes only the `COMPLETE.json` of each. The code is kept in the work directory, not in the repository.
