# cohort_database_sizes

Internal ID: A01 (previous folder name A01_database_sizes)

## Question

This pipeline produces: acquisition_sizes_726.tsv (710 rows from A01 audit + 16 participant-15 CapScan rows) and human_reference (20,431 Swiss-Prot entries): Fig1 panel b; Fig2 CAMPI reservoir 375,333,654; SFig01 panel a; SFig15 233-375 M range; abstract/discussion lines 35,181,193; Methods line 215; Supp Methods lines 665-713

Used by: figures/Fig1, figures/Fig2, figures/SFig01, figures/SFig15, text_statistics/abstract_discussion, text_statistics/methods, text_statistics/supplementary_methods_6_7.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `audit_v2.py` | Audit actual FASTA entries and render the approved four-cohort comparison. |
| `slurm/audit_v2.sbatch` | launcher |
| `sizes.py` | Per-acquisition FASTA size rows with the A01 producer's own scan() (frozen/audit_v2.py, imported unchanged). |
| `assemble.py` | Insert the 16 participant-15 rows into the Source Data tables (CapScan block re-sorted by sample, as in the 710 |
| `audit.py` | Audit actual FASTA entries and render the approved four-cohort comparison. |
| `slurm/00_lake.sbatch` | NO --partition DIRECTIVE: directives are parsed before the shell runs, so they |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `acquisition_sizes_726.tsv` | `fd9ade3bb4ebd2ee6c129f69819b9697cf20abd7c10c77a12644bef25396e386` |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

Requests in the launchers (`#SBATCH` lines):

- `slurm/audit_v2.sbatch`: no resource lines
- `slurm/00_lake.sbatch`: #SBATCH --nodes=1 --ntasks=1 --cpus-per-task=64 --mem=400G --time=12:00:00

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `OUT` | `pathlib.Path(__file__).resolve().parents[1]` | no comment in original code | not stated in code; see sizes.py |
| `SRC` | `A01 / 'frozen/audit_v2.py'` | no comment in original code | not stated in code; see sizes.py |
| `OUT` | `pathlib.Path(__file__).resolve().parents[1]` | no comment in original code | not stated in code; see assemble.py |
| `HUMAN` | `20431` | no comment in original code | not stated in code; see assemble.py |
