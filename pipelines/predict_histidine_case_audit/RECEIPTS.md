# Receipts

Original Slurm jobs, measured resources (from `sacct`, via `sacct -P -j`), and sha256 of every output file the figures and tables can draw on. Paths are relative to `LAB_ROOT`. Per-acquisition directories that carry a `COMPLETE.json` are summarised by the sha256 of the sorted list `<sha256 of COMPLETE.json>  <path>`; each `COMPLETE.json` itself records the sha256 of every output of its acquisition. Hashing was done on Slurm with `hashjob.py` (see the end of this file).

## Slurm jobs

| job id | stage | states | CPUs allocated | memory requested | wall time max (min) | peak RSS (GB) |
|---|---|---|---|---|---|---|
| 5940795 | case_audit_pilot | COMPLETED: 1 | 2 | 16G | 0.3 | 0.0 |
| 5940796 | case_audit_full | COMPLETED: 1 | 2 | 16G | 0.5 | 0.0 |

## Output hashes

### `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/spectrum_case_review/full/independent_case_review.tsv`

| sha256 | file |
|---|---|
| `c09af4794aa26f91aa213b67bada4e3709a24166b58050d0c64056ed380e7881` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/spectrum_case_review/full/independent_case_review.tsv` |

### `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/predict_spectrum_review/full`

| sha256 | file |
|---|---|
| `029d9eecf5c29d88acb7a40e2a7940fa198a8b1f907725c3f009e14c5991b53d` | `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/predict_spectrum_review/full/COMPLETE.json` |
| `26d1626714ba453dc783b0ba9d091075c251dc79777d6f12207f140b518f8517` | `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/predict_spectrum_review/full/independent_case_review.tsv` |
| `217e7b05b3dd2c7065e0c3f9a701ac92e04a389f7100726474b16cc896e5f5ef` | `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/predict_spectrum_review/full/measured_case_review.pdf` |
| `4725845f0a4b866d39d2c9eb622f855c034aa73a65cd475c8ef58999ab4444a6` | `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/predict_spectrum_review/full/measured_case_review.png` |

### `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/predict_spectrum_review/MANIFEST.json`

| sha256 | file |
|---|---|
| `45bcd329d2f8cb8964919f20c93d885e6fe12c4498181bb9ba54e1ad5f7e08a1` | `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/predict_spectrum_review/MANIFEST.json` |

## Released code versus origin

| file | origin (LAB_ROOT relative) | sha256 at origin | sha256 here |
|---|---|---|---|
| `audit_predict_spectrum_cases.py` | `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/predict_spectrum_review/frozen/audit_predict_spectrum_cases.py` | `ae0edcea5f876f4394d50ab36a775769fad0863759fa99712d8a9176bc1d4320` | `ae0edcea5f876f4394d50ab36a775769fad0863759fa99712d8a9176bc1d4320` |

Static checks: `python -m py_compile` on every `.py` file; `tools/scan_release.py`.

## Hashing code

The hashes above were produced by a Slurm job (4 to 8 CPUs) that walks each root once, hashes every file with SHA-256 in 8 MiB blocks and, for roots with more than 20 per-acquisition directories, hashes only the `COMPLETE.json` of each. The code is kept in the work directory, not in the repository.
