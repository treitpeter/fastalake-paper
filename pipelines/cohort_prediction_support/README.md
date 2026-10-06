# cohort_prediction_support

Internal ID: A09 (previous folder name A09_prediction_support)

## Question

This pipeline produces: acquisition_support_summary_726.tsv for Fig1 panel d and Results line 63; worked_acquisition_prediction_support.json for 4,258 of 29,663

Used by: figures/Fig1.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `analysis.py` | no module docstring in original |
| `collect.py` | Combine validated acquisition outputs and independently verify displayed counts. |
| `summarize.py` | Acquisition support rows from A09 protein_support.tsv, using the record logic of |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `acquisition_support_summary_726.tsv` | `6f81da0a03893bcdf2ed5406ddb0e52eb16791cfd46f55a3253b31c65ac7dec2` |
| `worked_acquisition_prediction_support.json` | not found on disk |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `OUT` | `pathlib.Path(__file__).resolve().parents[1]; R = OUT / 'a09_run'` | no comment in original code | not stated in code; see summarize.py |
| `ORIG` | `A9 / 'results/acquisition_support_summary.tsv'` | no comment in original code | not stated in code; see summarize.py |
