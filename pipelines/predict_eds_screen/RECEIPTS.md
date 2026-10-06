# RECEIPTS

## Original Slurm jobs

| role | job |
|---|---|
| baseline replay | 6044262 |
| refit (all EDS families, separate adjustments, individual follow-up) | 6048399 |

## Code

| original (LAB_ROOT-relative) | file here | sha256 original | sha256 here |
|---|---|---|---|
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/common.py` | `common.py` | `a675dd91ad8da1751d9aeea6b531159d1570a9a92bcbfe4c3d5d3a1cf4e617cd` | `2dd3a2c8ebe07350f6486469e6f2b4f8f69ae85bc6a02e2d5da6a3d7639bb7c3` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_biochemistry.py` | `eds_biochemistry.py` | `988ec3220b9ee42efac383bcbf4b5931d424e52f472878c5e0c973900c6fb686` | `988ec3220b9ee42efac383bcbf4b5931d424e52f472878c5e0c973900c6fb686` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/analyze.py` | `analyze.py` | `9925977db928a866b73ef424b34de45a99d04cc571c160b83f099e5d496b1a1a` | `8b196a44fc5ffc9c87780c6dc75c6a10ed433a02b0ed570f67fd347d7f1a5de5` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/individual_followup.py` | `individual_followup.py` | `86ecade6b17e9291a38b57b72cf671af01a7dfa8fb627d5f2da0f52eb85fadfb` | `86ecade6b17e9291a38b57b72cf671af01a7dfa8fb627d5f2da0f52eb85fadfb` |

## Inputs (sha256 recorded by the original run)

| input (LAB_ROOT-relative) | sha256 |
|---|---|
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/inputs/matrices_repaired/PREDICT_KO_microbial_broad.tsv.gz` | `7184e3da856ef7e3a3eda180d21691dc4a8615f06fb5dba11df3c64238921474` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/inputs/matrices_repaired/PREDICT_Pathway_microbial_broad.tsv.gz` | `9b7a89360209f0803ed84693679873533ac13f4b40d6192b690fccbc2d9dbf26` |
| controlled input `heatmap_85_patient_order.tsv` (not in repository) | `4ee42bf87fdc5b73878ff70d20c2f0d995bacd8e35ea378baf5738a45256384d` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/PREDICT_repaired__original/source_data/Figure4/heatmap_85_values.npz` | `0cf55ba3f6ce26b4051ad33fb7b85fb066ed37c3a0291cb3282f8abfbf9be886` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/PREDICT_repaired__original/source_data/Figure4/heatmap_columns.tsv` | `60a773ff51861574930a4a17da3c55143f2c02a05872010c0eed2521c62490ad` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/PREDICT_repaired__original/source_data/Figure4/heatmap_features.tsv` | `a5acdfec0e270dc04097699c610d01dfae518a33bbb09e932138fd96588d9eb1` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/PREDICT_repaired__original/source_data/Figure4/heatmap_values.npz` | `862037a97826a3e729bb52d37e37d675fe9024e61d84c39c9e5ec269e0d097b5` |
| `analyses/enterotypes_2026-09-07/PREDICT_metabolite_annotations.tsv` | `14e044fbbf5f6dedcd68231a6b25053f38ec4e4b6fbc2d7ad953c38b580f0121` |
| `analyses/enterotypes_2026-09-07/PREDICT_metabolite_intensity.tsv` | `0663b71545e31e677368c1fcbfe7d2c503e987f4e5a4eff677da6d8848418cc9` |
| `denovo_vs_tpm/metab_reverse_2026-09-14/reference/compound_reaction.tsv` | `09084463946eeb27bcccb6871c489248935d676cd7d1d66772e0f10d9d2e5689` |
| `denovo_vs_tpm/metab_reverse_2026-09-14/reference/reaction_ko.tsv` | `539b237200129db56c91ab5308f0e76d7ea58500d4c098b415808ad6cf00a2e1` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/ADDITIONAL_LINKS.json` | `2bae424ab02b718a7a2b592928214be7ca103d7a373fce88c611519ca766374d` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/INPUTS.json` | `1d7d15d2a79bddedb6cacafd8d66c929b253c8fc049f1b805ab2d06a5d45c337` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/all_EDS_associations.tsv` | `a46342612d7f4bf90c538f39cd4854d641c8c8fe1db55f223c48a82cd1b1a4d2` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/all_linked_associations.tsv` | `33264126a194d835877b8c095ce85c9c00b6eeddc159a1e3152919f4c4095f00` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/analyze.py` | `9925977db928a866b73ef424b34de45a99d04cc571c160b83f099e5d496b1a1a` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/augmented_pair_screen.tsv` | `8c81ed1bded410d86ddc450c8e7e1833eed334a20d5d1e81bec27e6ba71d2f30` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/individual_followup.py` | `86ecade6b17e9291a38b57b72cf671af01a7dfa8fb627d5f2da0f52eb85fadfb` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/separate_adjustments/TARGETS.json` | `6270ac79a4ae241989a81df47e920316d3dfcb95b2925779a8b135de2d380a5d` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/separate_adjustments/all_models.tsv` | `8991f713a7597c77ff7cdf2faa7b764630e6355a3c1902b350b2d673d2ad461e` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/eds_reaction_check_2026-09-30/separate_adjustments/run_analysis.py` | `e27b62da1d938229a6cea7688d7ecbf335417fe0c3b62c16b533bb08832597a3` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/input_data/KO_profile.tsv.gz` | `164fe6cde20cec1e6cc6803c7b2ad50a3f373f02e411e211769670b0f39cec02` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/input_data/Pathway_profile.tsv.gz` | `8453f2b543726afcdb30a7b29159c7cf02f009a7598dd72e9f777dd2396c40f7` |
| controlled input `heatmap_85_patient_order.tsv` (not in repository) | `43dc72c4b31e4bcf44c2c6bfbcc16a0a713ea47a8f377e86c87280fc5327b9ba` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/input_data/heatmap_85_values.npz` | `11c171206d200b9051c0f785041288bbc621be61695baae1bfa8eaffb88b7b37` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/input_data/heatmap_features.tsv` | `b883e1659743f5437c3d1d3981abd07a36d1ae41ee5d87b6321596fe90968e51` |
| `reviews/ipath3_v6.2_MM_2026-09-13/source_data/Simple_PREDICT_metabolite_id_audit.tsv` | `5f6875b82b2b4baeaa0fa8ab993e026df4216a29c0ed12adca3dafb8e9aec927` |
| `simplified_grouping/final_revision/Fig5/kegg_cache/compound_to_pathway.tsv` | `692cc003e36f7753beae7fa45f33c4ea9cc28d2b4f540bc78a2d118d4500bfb0` |

## Outputs (sha256 recorded by the original run)

| output (LAB_ROOT-relative) | sha256 |
|---|---|
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_baseline/individual_followup_results.tsv` | `9c0353d103d4129d03422307dcd85fb338a440f5004b0ee4ad5abcb7ea440bbe` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_baseline/ADDITIONAL_LINKS.json` | `2bae424ab02b718a7a2b592928214be7ca103d7a373fce88c611519ca766374d` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_baseline/FOLLOWUP_COMPLETE.json` | `00dceec1183e5b88c0ebf64ae795301d72990d9348dc4bdd36285bff87d7c10a` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_baseline/analyze.py` | `9925977db928a866b73ef424b34de45a99d04cc571c160b83f099e5d496b1a1a` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_baseline/augmented_pair_screen.tsv` | `5eab9a9b6d057a58e148f725ed2e77a51e27fe22f37108436d07af61f88d45a8` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_baseline/joined_evidence.tsv` | `27400e2b3872253ab3fd1fcdb637d73af778c166d3c033085bc0110b8c2bae76` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_baseline/all_linked_associations.tsv` | `1d3a8a3fbc49e9f20cb1094d2bda834e91f59ba6a6c7002167be1e202e3c66c9` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_baseline/reaction_links.tsv` | `d161a5961ed033d9bd95c6e4c4e9aee2435426a7474f0282cc06c3e2918c9cdf` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_baseline/individual_followup.py` | `86ecade6b17e9291a38b57b72cf671af01a7dfa8fb627d5f2da0f52eb85fadfb` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_baseline/INPUTS.json` | `7b7bbf3283e4cee20897c5e5f33244eaf47402ad72ff959395fbd6b5367bcf9b` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_baseline/PILOT.json` | `d5c87acfd4fbe8ae8a67bb42a6dc55b00fc612cd7f31ac814f1afa3c6310caf1` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_baseline/all_EDS_associations.tsv` | `fd3c5d4ad21a0d7476f488350144b6bf4e8d6832514b0d0fea05c337d279a869` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_baseline/heatmap_metabolite_linked_evidence.tsv` | `243812daa09a955c7b0b761a2e043f4fcb0ebe782571a6504df8319254e3e42c` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/individual_followup_results.tsv` | `d80a2f061b01193d7fe3f534fa1b206bb69530cdd45ffeb16d77cb975162f93f` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/ADDITIONAL_LINKS.json` | `2bae424ab02b718a7a2b592928214be7ca103d7a373fce88c611519ca766374d` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/FOLLOWUP_COMPLETE.json` | `81bbdc84624bce07456c1832d42e627df333e51bdbbdee7295110550d1210c8c` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/analyze.py` | `9925977db928a866b73ef424b34de45a99d04cc571c160b83f099e5d496b1a1a` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/augmented_pair_screen.tsv` | `d7e59ab759d62aa5b2068ac5fbe0b774ece288ad3894bd633cfd219cac2a3967` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/joined_evidence.tsv` | `a782be33ab19269bd763943d1c710ad1f742540ffca2197f8ab8c9dc84ddefd5` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/all_linked_associations.tsv` | `6eeb08652c094290626707deaad0ab0fda5dafdf98f9218a723130b803d23ced` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/all_models_comparison.tsv` | `7ad49b24ef1027722dd732f505df46822cf45f91e03e47cdd6c8c7309863711d` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/reaction_links.tsv` | `b6c332d2d884b13f2077410c870eebe09971fc5277d35a9490e026429e13dad2` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/augmented_pair_screen_comparison.tsv` | `55ccca7a1d9f5bb06a626a07335b54a1e553a9f2ae651fa438208ce2f4e81157` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/individual_followup.py` | `86ecade6b17e9291a38b57b72cf671af01a7dfa8fb627d5f2da0f52eb85fadfb` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/INPUTS.json` | `52d52c22779dacec81206fa8a613b09f02d6a6b50e45a5f7fd05b30ce3f0c22b` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/all_EDS_associations_comparison.tsv` | `22115899b85b943e9c63eed223a9a75eb3900aa60639f68e64e8a1fca0acc7a9` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/PILOT.json` | `97786d41f252223c0a13d7cb10a7ca6d2f7855d123ab7e36a9846550dc78db59` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/all_EDS_associations.tsv` | `dcd15919dde69cce09d3fdd3e96c5517e35f1786b824bb5575ca7dc0131008dd` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/eds_refit/heatmap_metabolite_linked_evidence.tsv` | `92f51edc07b7c39d661248503c5a2d7bfc90d7821b139bd8bbdca1d25c8fbeff` |

## Path aliases

Lab path segments are shortened in code and tables (`jess_integration_2026-09-29` to `integration_2026-09-29`, `v10_mm_philip_comment_audit_2026-09-22` to `v10_comment_audit_2026-09-22`, `coauthor_comment_check_v25_2026-10-05` to `comment_check_v25_2026-10-05`). `FASTALAKE_DATA_DIR` must hold the frozen inputs under the shortened names; the original names are the lab directory names.
