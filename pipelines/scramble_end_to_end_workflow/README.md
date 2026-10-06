# scramble_end_to_end_workflow

Internal ID: none; previous folder name UNASSIGNED_terminal_scramble_workflow

## Question

This pipeline produces: scramble_end_to_end.tsv (6.67 M, 12,606, 29,663, 5,841, 525, 0) for Fig1 panel a and Results lines 53, 59, 61; scramble_lookup_comparison.tsv (selected-FASTA sizes, accepted peptides, query-match fractions)

Used by: figures/Fig1, text_statistics/supplementary_methods_6_7.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `search_worker.py` | no module docstring in original |
| `submit_search.py` | no module docstring in original |
| `search_manifest.json` | launcher |
| `build_evidence.py` | no module docstring in original |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `scramble_end_to_end.tsv` | not found on disk |
| `scramble_lookup_comparison.tsv` | not found on disk |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| none found as top-level constants | | | |
