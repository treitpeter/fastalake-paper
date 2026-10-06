# cohort_spectrum_split_selection

Internal ID: A73 (previous folder name A73_spectrum_split_selection)

## Question

This pipeline produces: Per-acquisition FULL / JOINT (self) / DISJOINT accepted-peptide counts and identity overlaps on odd-scan spectra, 710 acquisitions + 726 re-run (ACQUISITION_REVIEW.tsv, SEARCH_MODEL_REVIEW.json); SFig06, Methods line 267, lines 755-761

Used by: figures/SFig06, text_statistics/methods, text_statistics/supplementary_methods_6_7.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `pipeline.py` | no module docstring in original |
| `submit.py` | no module docstring in original |
| `slurm/run.sbatch` | launcher |
| `slurm/gate.sbatch` | launcher |
| `review.py` | no module docstring in original |
| `collect347.py` | no module docstring in original |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

See RECEIPTS.md.

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

Requests in the launchers (`#SBATCH` lines):

- `slurm/run.sbatch`: no resource lines
- `slurm/gate.sbatch`: no resource lines

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `ROOT` | `next(p for p in W.parents if (p/'00_START_HERE_V7.md').exists())` | no comment in original code | not stated in code; see pipeline.py |
| `CAN` | `ROOT/'canonical'; ARMS=CAN/'publication_data/arms'` | no comment in original code | not stated in code; see pipeline.py |
| `PAR` | `CAN/'rust/parsimony_engine/target/release/parsimony_engine'` | no comment in original code | not stated in code; see pipeline.py |
| `SAGE` | `ROOT/'sample_specific_merge_all_raw/2026-09-11/bin/sage'` | no comment in original code | not stated in code; see pipeline.py |
| `COHORTS` | `['HSTOOL','CAMPI_MOCK','PREDICT','CAPSCAN']` | no comment in original code | not stated in code; see pipeline.py |
| `OLD` | `PRODREV.parents[1];BASE=OLD.parents[2]` | C347: production locations unchanged | not stated in code; see collect347.py |
