# Receipts: capscan_crossfit_10fold

Original locations are relative to LAB_ROOT.

## Slurm jobs

6051943, 6065335, 6065336, 6065344

## Output files

| file | bytes | sha256 |
|---|---|---|
| `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration4/crossfit/run/PER_ACQUISITION.tsv` | 464180 | `f7616394a3ffd829c92c2256244f5ce174116f7596b579f5db06b7b56a95a6e3` |
| `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration4/crossfit/run/SUMMARY.tsv` | 13984 | `abd33b20058b9cec3374d5b63f2afc9fd8718eeccf79e92421e17e54b070d85a` |
| `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/public347/crossfit_capscan347/run/SUMMARY.tsv` | 13986 | `4449fbe156c1abc2177db2e9c3d7d3bac378a5802a948b9037b150508e87fe0f` |

## Code files in this folder

| file | sha256 |
|---|---|
| `analyze_final.py` | `835e0c93c04c41d81a60dc6fbbe1a61ffe25055b3a918af175d9cf1abdd0b2da` |
| `common.py` | `b51505b0712bf90b81f8a1474ff2ae705cd1282cfd07aeb0da169b1796dc35cf` |
| `recovery.py` | `b2d904fbed9547337ff0aac00a4e76a9c01d6985dc3bf1aec0550a115811c5b8` |
| `run__analyze_final.py` | `e681881b3102b09556d1078656aece62cbc33cb486e15fd825750213ba4e7be4` |
| `sage_lda.py` | `33c5e3d0d278c54cd70e38fe9ed1972f7968272dde93843a885d0f3d489a1aa1` |
| `sage_q.py` | `829a7e59bbd9dc86c8ea8647cccadfd12aa96f442b5db4ad81d108cb59548c10` |
| `worker.py` | `3cd642c5b2d0f0bff1eec298745f3eab98a88ec38ec6ccbdfe7f174900f2d6ba` |
| `slurm/final.sbatch` | `3ef323bedc3c99efb1334ade1ae06cfb4bb1815e9c171542714c2774ff9d0c9c` |
| `slurm/run.sbatch` | `e516147025e56efc5f23ae501750fffbc58915ef75516e9bbc864c011ebaa282` |

## Cross-fitted PEP sensitivity (added 2026-10-06)

Slurm jobs: 6057983 (validate), 6057984 (pilot), 6057988 (negative control), 6057990 (array, 20 tasks), 6058010 (analysis); all COMPLETED.

| file | bytes | sha256 |
|---|---|---|
| `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C1_pep_primary/crossfit/CF_PEP_PER_ACQUISITION.tsv` | 230937 | `41e3cfd156a40cd34a318ea9ba40fc1773d08840428f852f7ae57e2d16e2de48` |
| `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C1_pep_primary/crossfit/CF_PEP_SUMMARY.tsv` | 3353 | `e43dac3c0bb212ac762e0fd5f44fe2bf65328a129c2f41abe77ebb741500e7e2` |
| `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C1_pep_primary/crossfit/CF_PEP_SUMMARY_BY_REGION.tsv` | 4170 | `5ccc94fe30375b8da581aa1f251ec1aebf5df8cef98c1e6e38e2e00054af0f72` |
| `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C1_pep_primary/crossfit/validate_0.tsv` | 8161 | `027108347eafc53fc40fca94101a9ba69f18ee8820bbc4603259e71b2725202c` |
| `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C1_pep_primary/crossfit/PEP_MODEL.json` | 183 | `d1c92336f9c92ffc8569357f62f3d96d6a921cca755733cebcf8fc89a95c061f` |
| `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C1_pep_primary/crossfit/per_acq/*.json` | 578 files | not hashed individually |

Source Data written by export_pep_source_data.py: crossfit_pep_per_acquisition.tsv `d46e9712efd7a4dc0e4d8393fc0ab1d916c34863da9e0bc41f7bf66f03790839` (summaries are byte copies, hashes above).
