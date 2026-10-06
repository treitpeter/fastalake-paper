# healthy_stool_selection_sweep

Internal ID: A70 (previous folder name A70_selection_sensitivity_sweep)

## Question

This pipeline produces: Merged with the former selection-grid entry (same table, selection_grid_acquisitions.tsv). 19-arm x 58-acquisition selection sweep (SupplementaryTable_selection_sensitivity.xlsx Per_acquisition and Arm_summary) for SFig02 c,d and line 679; 24,560 queries in Fig1 panel a and Results line 53

Used by: figures/Fig1, figures/SFig02, text_statistics/supplementary_methods_6_7.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `sweep_tables.py` | Per-sample and per-arm tables for the HSTOOL length / score / k-mer sweeps. |
| `submit_arm.sh` | One HSTOOL sweep arm through the REAL chain driver (gate -> evidence -> persample -> |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `selection_grid_acquisitions.tsv` | not found on disk |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `ARMS_ROOT` | `f"{CANON}/publication_data/arms"` | no comment in original code | not stated in code; see sweep_tables.py |
| `MOD_RE` | `re.compile(r"\[[^\]]*\]")` | no comment in original code | not stated in code; see sweep_tables.py |
