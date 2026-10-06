# clinical_common_fasta

Internal ID: A23 (previous folder name A23_clinical_common_fasta_synthesis)

## Question

This pipeline produces: PREDICT/per_participant.tsv and CAPSCAN/per_participant.tsv, concatenated and pseudonymised as common_search_participants.tsv (SFig07 b)

Used by: figures/SFig07.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `summarize_v2.py` | no module docstring in original |
| `slurm/CAPSCAN.sbatch` | launcher |
| `slurm/CAPSCAN_v2_control.sbatch` | launcher |
| `slurm/PREDICT.sbatch` | launcher |
| `slurm/PREDICT_v2.sbatch` | launcher |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

See RECEIPTS.md.

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

Requests in the launchers (`#SBATCH` lines):

- `slurm/CAPSCAN.sbatch`: no resource lines
- `slurm/CAPSCAN_v2_control.sbatch`: no resource lines
- `slurm/PREDICT.sbatch`: no resource lines
- `slurm/PREDICT_v2.sbatch`: no resource lines

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| none found as top-level constants | | | |
