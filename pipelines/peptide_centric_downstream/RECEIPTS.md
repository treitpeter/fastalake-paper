# Receipts: peptide_centric_downstream

Original locations are relative to LAB_ROOT.

## Slurm jobs

5973723, 5973693, 5973716, 5973758, 6058096, 6058097, 6058098

## Output files

| file | bytes | sha256 |
|---|---|---|
| `LAB_ROOT/MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/AGGREGATE_VALIDATION.json` | 6925 | `9d3bf2f0616f301b69e98398372b501f1e710bda245bde87659a1025814d45b5` |
| `LAB_ROOT/MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/CAPSCAN_own_stool_donors.tsv` | 1876 | `d018052d7af40fde4c6f3c2164e54b793672f24d7b228d44d4d04eb47e4d591c` |
| `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig5e_f/run347/NUMBERS.json` | 1012 | `55c5ccf007cc3aedb21593d18851ed01cbc39f42084f08279872ccd3bb705486` |
| `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig5e_f/gated_joint_347/CAPSCAN_gated_joint_stats.json` | 1105 | `12ca54fa3202a2d7bf71b39c9d2be296b6659d2d43163c867771a170e3beb543` |
| `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig5e_f/design_347.tsv` | 96708 | `abe0fd5290c2ec2f6db4b94ba15c7171508caeddc5c82e34aab32b3bc2440c09` |

## Code files in this folder

| file | sha256 |
|---|---|
| `build_design_347.py` | `5ccf811f7547428565caa08f9e80416b8e6cfe769fc32bd67308fc4fa1e3aad7` |
| `build_gated_joint.py` | `c5735479f00b00fc03614dc1ffc4696d6079cde94a4d4acd358ec8c1edb72d3b` |
| `build_supp_9_10_19.py` | `fd9b6bf9502ebc5715fe8c5055a2b850baf8446cd12e515842e703e9d465575e` |
| `collector.py` | `f5d6e0ba3765bfa931b69e1631462cb4023b56e8d1a1c6e6fa1060d5439a6c9e` |
| `common_pepcentric.py` | `bc1f0a68b3e92c416ee08f35a5758444a1da84987b434cdc49b8fe552103d9a5` |
| `compare_control.py` | `6d1f0c6ac4eed06a0da26c508c1ca11765052b75124b8410b0f5535fa9d54ccb` |
| `derive_hstool_pc.py` | `cafe0e11f4a55bdf1cf6ebb95b6834ec8d5591a3b9484b47f8a6b40978b74783` |
| `fig5ef_run.py` | `14d857cd68fd06b88a353eabd862150d5d258f6b2ce5d7a6d4d074c446fe00b8` |
| `hstool_pc.py` | `20cf0b930b2efac12bfcbb9dfa75237507a01b04a2547398b71c90349b400e9f` |
| `make_leader_tables.py` | `5fffdec2b169ca633ba11e10947affc3d435dc0f70741ca7f0ade41896bfe616` |
| `scripts__build_design_347.py` | `5ccf811f7547428565caa08f9e80416b8e6cfe769fc32bd67308fc4fa1e3aad7` |
| `search_pc.py` | `60295063c2e7e395acdba6eb2c7d87e53c934a4024f3a45a973a7488562c1bff` |
| `slurm/hstool.sbatch` | `6873b6d5c137703ef1ef2930b01fb809f13faf7a83fac83deca0e1c977e7c0ab` |

## Direct-assignment recount of Fig. 5f (v30)

Slurm job 6074394 (p.hpcl92, 2 cores, 8G). Code: `strict_subset_direct_recount.py` = `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/grouping_recount/recount.py` (sha256 `02ef747a1bda5cfa5d95bfb4cecc933cdd7df9bd6734e6b099a5d54bda42737e`) with docstrings added.

| file | role | sha256 |
|---|---|---|
| LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig5e_f/run347/results/contributor_method_cells.tsv | input | 42105a1a22c4cbb6bb658c8b3ad8d02e476fa5b200d74c580be8a77c24d076f6 |
| LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig5e_f/run347/results/contributor_eligibility.tsv | input | 9aac5aa7fcc04153dbd2ab6e7f8681fc4324e5db6ec4b4d49f0ce8022f9ecdfb |
| LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/grouping_recount/recount.json | output | 6bc4ee17b3fc7daa7bdfc3f9797a408f73c380701d695687bee2172eeeca9100 |
| LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/grouping_recount/direct_resolved_per_participant.tsv | output | 4c1becc57677dfd48cc0678a1e518129d7938a529c174b48b77dfed683ed045f |
| LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/grouping_recount/direct_strict_per_participant.tsv | output | 1efa6eee231ec63231ebe09911956637ff1321d9a0149cb47a872e3260b6aa3f |
| LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/grouping_recount/ctrl_resolved_both_per_participant.tsv | output (control) | 0cd0c2aad5bba1be3fe4b6043de2b368a337a63e8fd9e2d3324c90121d967708 |
| LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/grouping_recount/ctrl_strict_both_per_participant.tsv | output (control) | 411bfe9e22b09bdeebe00515a79a44d8c2ba14b7350858d33d970d6c354d4fb1 |
