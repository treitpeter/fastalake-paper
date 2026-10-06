# bag_mass_paired_nulls

Internal ID: A04 (previous folder name A04_bag_mass_paired_nulls)

## Question

This pipeline produces: per-acquisition NULL_COMPARISON.json receipts (exact, bag, mass, combined routes, two shuffle seeds), frozen 60-receipt snapshot; all sheets of STab6

Used by: tables/STab6.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `run_bag_nulls.py` | Compare matched exact/bag/mass nulls from an existing full, route-tagged scan. |
| `guard.py` | no module docstring in original |
| `finish.py` | no module docstring in original |
| `submit.py` | no module docstring in original |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

See RECEIPTS.md.

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `NULL_ARMS` | `{` | no comment in original code | not stated in code; see run_bag_nulls.py |
| `ALL_ARMS` | `dict(ARMS, **NULL_ARMS)` | no comment in original code | not stated in code; see run_bag_nulls.py |
