# RECEIPTS

## Original Slurm jobs

| role | job |
|---|---|
| refit | 6042666 |

## Code

| original (LAB_ROOT-relative) | file here | sha256 original | sha256 here |
|---|---|---|---|
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/clinical_refit.py` | `clinical_refit.py` | `b496050501da5075efdf4aff75127bd357b2d4aceb277cb2209399f5e3735072` | `6e5a4c63edf0c8df3c2b4f81c29a1930b3ba410443cc48dc8d83a18401f76a4d` |

## Inputs (sha256 recorded by the original run)

| input (LAB_ROOT-relative) | sha256 |
|---|---|
| controlled input `acquisition_accounting.tsv` (not in repository) | `1a2dac477320d9d7362b572db683369adab168879df0dab40ce901fe621a43a5` |
| controlled input `clinical_group_global.json` (not in repository) | `dc82211971a6a540b56de76bf27559fd0c3fa1fc8f308b528c242ca110249e15` |
| controlled input `clinical_group_mean_intervals.tsv` (not in repository) | `419e8bea0b8e5e5854be96f544433abca6ca09c10e3575f4f1c801adce8229b4` |
| controlled input `clinical_group_tests.tsv` (not in repository) | `23f8aea5f57402d7346bafb0357578ca205068c68c19a4882de17f1cdd2a66b5` |
| controlled input `clinical_paired_statistics.tsv` (not in repository) | `2745ab1c21a9ced6c299e6a7827e607a839faa9a487847a2df9c2860f1a36770` |
| controlled input `paired_antibiotic_profiles.tsv` (not in repository) | `c41f5582e2ae54328c24d31b0414c85263c5d904f795e10de9a0fc85feda9563` |
| controlled input `visit1_patient_profiles.tsv` (not in repository) | `3099f78b6de461cb154cd9a831d1c5166edea688dc354291d1eab9910f50d993` |
| controlled input `derive_clinical_pc.py` (not in repository) | `90c5c4a1c15c4e869f80dc7a7456d4ec751e628caff156045f8204b206841ac8` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/profiles/PREDICT_noclust/COMPLETE.json` | `a0334cb944224790b3b28a7055e8d9c6a5988f65988029c3f11dd6a04e91e903` |
| `review/v15_integration/integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/profiles/PREDICT_noclust/taxonomy_signal.tsv` | `5a659ce13c30dd949c5b753898e34d2484fe2adbca794795c57a7a2e31081aab` |

## Outputs (sha256 recorded by the original run)

| output (LAB_ROOT-relative) | sha256 |
|---|---|
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/clinical_refit/clinical_group_mean_intervals.tsv` | `35a9c5421cc08b30dffccebe6ee2246b5eb41f508a43061cd17371c4b696e2be` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/clinical_refit/clinical_group_tests.tsv` | `20c3c52a74c406fdeff9f0681f5f133f3eb43bdc200d8e650ec978e3d529d007` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/clinical_refit/clinical_group_global.json` | `f4a2bf9d34a9eedbd59bd0f8c4a42eb9a6a06bf1937e0bbe5a531404e30d9f38` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/clinical_refit/acquisition_accounting.tsv` | `3880642bef93ccb390ff7be86b30a4aec715134b4f96c63ce0d5c6a4f5deaa23` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/clinical_refit/clinical_paired_statistics.tsv` | `9518b0117267eb3198b1aafe785dab47b53b020bd4deb0d9d2fb24f67ccc3bc7` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/clinical_refit/visit1_patient_profiles.tsv` | `68a3f3de774084016f61e4bad05a13f5d7e6c6ff610822847c8905d3e208f3bb` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/clinical_refit/group_method_comparison.tsv` | `1c644fb50453332db9ea44d5abc66e6b60f9018d94658a5820d7d0358bb662b1` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/clinical_refit/paired_method_comparison.tsv` | `e2282b692c28808eed07db0fb68ceb281048fbeeaedc57bffe96c0a1720c2e57` |
| `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/clinical_refit/paired_antibiotic_profiles.tsv` | `b99dd07a5cf75150f7cb6520eedcc816cefcd6e562dea9b457c309dfa95f27ec` |

## Path aliases

Lab path segments are shortened in code and tables (`jess_integration_2026-09-29` to `integration_2026-09-29`, `v10_mm_philip_comment_audit_2026-09-22` to `v10_comment_audit_2026-09-22`, `coauthor_comment_check_v25_2026-10-05` to `comment_check_v25_2026-10-05`). `FASTALAKE_DATA_DIR` must hold the frozen inputs under the shortened names; the original names are the lab directory names.
