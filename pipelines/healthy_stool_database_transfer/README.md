# healthy_stool_database_transfer

Internal ID: A15 (previous folder name A15_transfer_pool)

## Question

This pipeline produces: results/all_comparisons.tsv (transfer_pool_comparisons.tsv): Fig3 b-e, Results lines 99, 103, legend 107, Methods 739; 58 targets and donor assignments

Used by: figures/Fig3, text_statistics/supplementary_methods_6_7.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `digest.py` | Frozen CAMPI stage and reference-member audit; never perform a search. |
| `lib.py` | no module docstring in original |
| `parameters.py` | Matched public/cohort/metagenome identity summaries with declared weighting. |
| `peptide.py` | The canonical peptide definition. Import this; do not re-implement it. |
| `report.py` | Collect validated §4.4/§4.5 results without changing their frozen outputs. |
| `search.py` | Recount transfer/pooling raw results and verify database availability. |
| `slurm/audit.sbatch` | launcher |
| `slurm/pilot.sbatch` | launcher |
| `slurm/report.sbatch` | launcher |
| `setup_transfer_pool.py` | no module docstring in original |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `all_comparisons.tsv` | `a7618438c577a20b6037bf801b4daa89c181dcafb724066cfada229b06035947` |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

Requests in the launchers (`#SBATCH` lines):

- `slurm/audit.sbatch`: no resource lines
- `slurm/pilot.sbatch`: no resource lines
- `slurm/report.sbatch`: no resource lines

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `COL` | `'peptide_prediction_detokenized_unmodified'` | no comment in original code | not stated in code; see lib.py |
| `PYTHON` | `'/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python'` | no comment in original code | not stated in code; see lib.py |
| `Q_FIELD` | `"peptide_q"` | no comment in original code | not stated in code; see peptide.py |
| `Q_THRESHOLD` | `0.01` | no comment in original code | not stated in code; see peptide.py |
