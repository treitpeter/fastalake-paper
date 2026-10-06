# predict_clinical_refit

Internal ID: A91 (previous folder name A91_predict_clinical_refit)

After the profile repair, how does the human-assigned fraction and the quantified protein-group count of PREDICT stool change between antibiotic off and on within the same patient (42 paired patients, 84 patient-states) and across clinical groups (91 patients at the first visit)? Used by Fig4, STab3 and the Methods; the public output is the aggregate `clinical_paired_statistics.tsv`.

Status: code copied as it ran (path edits only), not re-run in this repository. See RECEIPTS.md for original job ids and hashes.

## Inputs

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

Patient-level inputs are controlled (see `controlled/CONTROLLED_INPUTS.md`): `acquisition_accounting.tsv` is read from `FASTALAKE_CONTROLLED_DIR`. The earlier pipeline (`derive_clinical.py` of the Figure 4 workshop) is the source of the transcribed `derive()` and is not copied; it chains broader figure inputs.

## Outputs

`clinical_refit/clinical_paired_statistics.tsv` (public aggregate); controlled and not published: `visit1_patient_profiles.tsv`, `paired_antibiotic_profiles.tsv`, `clinical_group_tests.tsv`, `clinical_group_mean_intervals.tsv`, `acquisition_accounting.tsv`. Hashes in RECEIPTS.md.

## Slurm resources and run time (original logs, `sacct`)

Job 6042666 (refit): 2 cores, 4G requested, peak RSS 2.4 MB (sampled), wall 2 s.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| paired population | patients with both antibiotic states; 42 patients, 84 patient-states | within-patient design | none |
| bootstrap | 50,000 resamples of patient differences, seed 20260907 | percentile 95% interval | none |
| sign-flip test | 100 x 2,000 random sign flips, seeds 20260907 + endpoint index, plus-one p over 200,001 | exact enough for p near 1e-5 | none |
| endpoints | human-assigned fraction (percentage points) and quantified simplified groups | the two claims made | BH over the two endpoints |
| baseline | first visit, 91 patients | one value per patient | none |
| baseline replay | the pre-repair outputs are replayed first (mode baseline) | proves the code reproduces the earlier numbers | none |

## Notes

Code edits relative to the original: lab absolute paths replaced by `FASTALAKE_DATA_DIR`; lab-tree roots located by a marker file replaced by `FASTALAKE_DATA_DIR`. No statistics changed. The scripts keep their original integrity gates (`verify()`, `gate()`, hash checks against `*_CODE.json` and `*_INPUTS.json` manifests of the frozen run). Those manifests contain absolute lab paths and are not in the repository; the hashes they held are listed in RECEIPTS.md. Because the code here carries the path edits described below, its own sha256 differs from the original (both are in RECEIPTS.md), so a re-run needs the manifests regenerated for this copy.`ROOT` (the lab tree) is read as `FASTALAKE_DATA_DIR`; baseline mode reads further per-patient tables from that tree and is controlled; only refit mode is the published path.
