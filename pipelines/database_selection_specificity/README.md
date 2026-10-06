# database_selection_specificity

Internal ID: A02 (previous folder name A02_lookup_specificity)

## Question

This pipeline produces: length_by_arm.tsv, five_seed_summary.tsv, query_source data: five-seed composition-preserving scramble and exact reservoir lookup of one healthy-stool acquisition (Slurm 5920360); Fig1 panel c, Results line 59, 36,588 in panel a; SFig02

Used by: figures/Fig1, figures/SFig02.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `worker.py` | no module docstring in original |
| `lib.py` | no module docstring in original |
| `report.py` | no module docstring in original |
| `composition.cpp` | include <fstream> |
| `slurm/lookup.sbatch` | launcher |
| `prepare.py` | no module docstring in original |
| `report_figure_v2.py` | no module docstring in original |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `length_by_arm.tsv` | `63594d955d35bd564a32fe68317274537519376858d3bc0daee1b09cc3b1cc29` |
| `five_seed_summary.tsv` | `ebb3de14cd36ab997a1e0e2245e882de30f822b3b06d78cbf42e8af79bf73cd6` |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

| job | name | cores | requested memory | peak RSS (GB) | elapsed | state |
|---|---|---|---|---|---|---|
| 5920360 | MM_A02_report | 2 | 8G | 0.0 | 00:01:22 | COMPLETED |

Requests in the launchers (`#SBATCH` lines):

- `slurm/lookup.sbatch`: no resource lines

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `RUN` | `pathlib.Path(sys.argv[1]).resolve()` | no comment in original code | not stated in code; see worker.py |
| `STAGE` | `sys.argv[2]` | no comment in original code | not stated in code; see worker.py |
| `COL` | `'peptide_prediction_detokenized_unmodified'` | no comment in original code | not stated in code; see lib.py |
| `PYTHON` | `'/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python'` | no comment in original code | not stated in code; see lib.py |
| `RUN` | `pathlib.Path(__file__).resolve().parent` | no comment in original code | not stated in code; see prepare.py |
| `ROOT` | `RUN.parents[1]` | no comment in original code | not stated in code; see prepare.py |
