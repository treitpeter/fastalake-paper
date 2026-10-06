# cohort_selected_fasta_unions

Internal ID: A97 (previous folder name A97_selected_fasta_unions)

## Question

This pipeline produces: Selected-FASTA unions of HSTOOL58, PREDICT306 and CapScan 331/347 rosters (union_contributions.tsv: PREDICT 555,987; CapScan 452,255 (331) or 460,977 (347)), line 751

Used by: text_statistics/methods, text_statistics/supplementary_methods_6_7.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `workflow.py` | Exact FASTA unions and a single joint Sage search per retained cohort/arm. |
| `run_union.py` | Build the exact-sequence union with the frozen workflow.build_union (byte-identical copy in ../code). |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `union_contributions.tsv` | not found on disk |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `TERMINAL` | `{'COMPLETED', 'FAILED', 'OUT_OF_MEMORY', 'TIMEOUT', 'CANCELLED',` | no comment in original code | not stated in code; see workflow.py |
