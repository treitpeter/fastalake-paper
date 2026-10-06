# prediction_parser_audit

Internal ID: A102 (previous folder name A102_parser_audit)

## Question

This pipeline produces: Per-file parser comparison for 726 prediction files (parser_audit_acquisitions_726.tsv), line 573; Slurm 6064887, 6064888 (original 6000982)

Used by: text_statistics/supplementary_methods_1_5.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `M07_parser_audit__collect.py` | Collect exact M07 parser comparisons; fail on incomplete or changed inputs. |
| `fixed726__collect.py` | Collect exact M07 parser comparisons; fail on incomplete or changed inputs. |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `parser_audit_acquisitions_726.tsv` | not found on disk |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

| job | name | cores | requested memory | peak RSS (GB) | elapsed | state |
|---|---|---|---|---|---|---|
| 6000982_0 | M07_parser_full | 1 | 4G | 0.1 | 00:11:13 | COMPLETED |
| 6000982_1 | M07_parser_full | 1 | 4G | 0.1 | 00:12:50 | COMPLETED |
| 6064887 | A7pa_new16 | 2 | 80M | 0.0 | 00:00:12 | COMPLETED |
| 6064888 | A7pa_collect | 2 | 500M | 0.0 | 00:00:06 | COMPLETED |

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| none found as top-level constants | | | |
