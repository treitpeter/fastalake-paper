# functional_allocation_sensitivity

Internal ID: A107 (previous folder name A107_functional_allocation_sensitivity)

## Question

This pipeline produces: Per-acquisition allocated intensity, CapScan own-stool distances, PREDICT distance matrices and healthy-stool per-pair correlations under four allocation rules (calls A84 transform/evaluate and A90 measure/pairs)

Used by: text_statistics/supplementary_methods_1_5.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `alloc_sensitivity.py` | GH-029: functional-allocation sensitivity on frozen annotated-origin consensus outputs. |
| `alloc_capscan347.py` | A3: CapScan part of the GH-029 functional-allocation sensitivity (iteration1/proposals_methods/alloc_sensitivity, |
| `collect.py` | A3 collector: allocation sensitivity, old (331, production job 6048208) vs new (347, CapScan rerun; HSTOOL/PREDICT unchanged). |

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
| `CONS` | `A / 'annotated_origin_consensus_2026-10-01/full'` | no comment in original code | not stated in code; see alloc_sensitivity.py |
| `FROZEN` | `C / 'profiles'` | no comment in original code | not stated in code; see alloc_sensitivity.py |
| `OUT` | `HERE / 'out'` | no comment in original code | not stated in code; see alloc_sensitivity.py |
| `LEVELS` | `['KO', 'OG', 'Pathway']` | no comment in original code | not stated in code; see alloc_sensitivity.py |
| `VARIANTS` | `['equal', 'unique', 'complete', 'unique_complete']` | no comment in original code | not stated in code; see alloc_sensitivity.py |
| `UNASSIGNED` | `'Unassigned_at_this_resolution'` | no comment in original code | not stated in code; see alloc_sensitivity.py |
| `HUMAN` | `'Human_assigned_signal_withheld'` | no comment in original code | not stated in code; see alloc_sensitivity.py |
