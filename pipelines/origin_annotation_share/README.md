# origin_annotation_share

Internal ID: A103 (previous folder name A103_partial_origin_share)

## Question

This pipeline produces: Annotated-origin rows and pooled intensity of the annotated-origin consensus (origin_rows_*.tsv), line 651; Slurm 6064869, 6064870

Used by: text_statistics/supplementary_methods_1_5.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `verify_proposals_methods__partial_human.py` | no module docstring in original |
| `fixed347__partial_human.py` | no module docstring in original |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `origin_rows_*.tsv` | not found on disk |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

| job | name | cores | requested memory | peak RSS (GB) | elapsed | state |
|---|---|---|---|---|---|---|
| 6064869 | A7p_347 | 2 | 100M | 0.0 | 00:00:01 | COMPLETED |
| 6064870 | A7p_collect | 2 | 100M | 0.0 | 00:00:01 | COMPLETED |

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| none found as top-level constants | | | |
