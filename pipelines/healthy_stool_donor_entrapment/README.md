# healthy_stool_donor_entrapment

Internal ID: A11 (previous folder name A11_control_reconciliation)

## Question

This pipeline produces: donor_recount.tsv and donor_(entrapment_)seed_means.tsv: donor-selected (58) and self-selected (5) explicit-peptide entrapment recounts, three seeds; Fig1 panel e, Results line 65; SFig05 panel d, lines 65, 425, 727-729

Used by: figures/Fig1, figures/SFig05, text_statistics/supplementary_methods_6_7.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `analysis.py` | Recount frozen donor and fusion controls; preserve their different estimands. |
| `donor_original.py` | Run a frozen-database peptide entrapment assay on a recorded HSTOOL transfer. |
| `slurm/donor_array.sbatch` | launcher |
| `slurm/report.sbatch` | launcher |
| `audit_fusion.py` | Independent accounting audit of frozen fusion searches; no FDR certification. |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `donor_recount.tsv` | `4e51dd75c37deed3f13c84ad2407032b01ae16d877efac9d1b60a174aafb5ea8` |
| `donor_seed_means.tsv` | `5c6b0af7c7224b7f1968b82a8a2c22ecae4d25679f691ef04c0957644c27ea66` |
| `donor_entrapment_seed_means.tsv` | not found on disk |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

Requests in the launchers (`#SBATCH` lines):

- `slurm/donor_array.sbatch`: no resource lines
- `slurm/report.sbatch`: no resource lines

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `THRESHOLDS` | `(0.001, 0.005, 0.01, 0.02, 0.05)` | no comment in original code | not stated in code; see donor_original.py |
| `THRESHOLDS` | `(0.001, 0.005, 0.01, 0.02, 0.05)` | no comment in original code | not stated in code; see audit_fusion.py |
