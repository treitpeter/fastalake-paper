# healthy_stool_pooling_controls

Internal ID: A94 (previous folder name A94_hstool_pooling_controls)

## Question

This pipeline produces: per_acquisition.tsv (arms own, pool_all12, pool_other11; hstool_pooling_per_acquisition.tsv) for the two counts at v28d line 743; search array 5879093, reduce 5879104

Used by: tables/STab2, text_statistics/supplementary_methods_6_7.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `prepare_record.py` | no module docstring in original |
| `run_task.py` | no module docstring in original |
| `slurm/search.sbatch` | launcher |
| `reduce.py` | no module docstring in original |
| `slurm/reduce.sbatch` | launcher |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `per_acquisition.tsv` | `474de2e78bbfc37fcf132930b688069df06909699b6f7cb9a0c12bf0c42222f8` |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

| job | name | cores | requested memory | peak RSS (GB) | elapsed | state |
|---|---|---|---|---|---|---|
| 5879093_0 | fig3_pool_controls | 16 | 8G | 5.7 | 00:01:12 | COMPLETED |
| 5879093_1 | fig3_pool_controls | 16 | 8G | 5.7 | 00:01:16 | COMPLETED |
| 5879093_2 | fig3_pool_controls | 16 | 8G | 5.8 | 00:01:17 | COMPLETED |
| 5879093_3 | fig3_pool_controls | 16 | 8G | 5.8 | 00:01:17 | COMPLETED |
| 5879093_4 | fig3_pool_controls | 16 | 8G | 5.8 | 00:01:18 | COMPLETED |
| 5879093_5 | fig3_pool_controls | 16 | 8G | 5.7 | 00:01:16 | COMPLETED |
| 5879093_6 | fig3_pool_controls | 16 | 8G | 5.8 | 00:01:14 | COMPLETED |
| 5879093_7 | fig3_pool_controls | 16 | 8G | 5.6 | 00:01:15 | COMPLETED |
| 5879093_8 | fig3_pool_controls | 16 | 8G | 5.8 | 00:01:17 | COMPLETED |
| 5879093_9 | fig3_pool_controls | 16 | 8G | 5.7 | 00:01:17 | COMPLETED |
| 5879093_10 | fig3_pool_controls | 16 | 8G | 5.8 | 00:01:19 | COMPLETED |
| 5879093_11 | fig3_pool_controls | 16 | 8G | 5.6 | 00:01:17 | COMPLETED |
| 5879093_12 | fig3_pool_controls | 16 | 8G | 5.9 | 00:01:13 | COMPLETED |
| 5879093_13 | fig3_pool_controls | 16 | 8G | 5.7 | 00:01:16 | COMPLETED |
| 5879093_14 | fig3_pool_controls | 16 | 8G | 5.7 | 00:01:13 | COMPLETED |
| 5879093_15 | fig3_pool_controls | 16 | 8G | 5.7 | 00:01:14 | COMPLETED |
| 5879093_16 | fig3_pool_controls | 16 | 8G | 5.7 | 00:01:19 | COMPLETED |
| 5879093_17 | fig3_pool_controls | 16 | 8G | 5.7 | 00:01:19 | COMPLETED |
| 5879093_18 | fig3_pool_controls | 16 | 8G | 5.9 | 00:01:12 | COMPLETED |
| 5879093_19 | fig3_pool_controls | 16 | 8G | 5.7 | 00:01:16 | COMPLETED |
| 5879093_20 | fig3_pool_controls | 16 | 8G | 5.8 | 00:01:17 | COMPLETED |
| 5879093_21 | fig3_pool_controls | 16 | 8G | 5.8 | 00:01:17 | COMPLETED |
| 5879093_22 | fig3_pool_controls | 16 | 8G | 5.8 | 00:01:17 | COMPLETED |
| 5879093_23 | fig3_pool_controls | 16 | 8G | 5.7 | 00:01:17 | COMPLETED |
| 5879104 | fig3_control_reduce | 32 | 64G | 47.2 | 00:01:18 | COMPLETED |

Requests in the launchers (`#SBATCH` lines):

- `slurm/search.sbatch`: #SBATCH --partition=p.hpcl93; #SBATCH --cpus-per-task=8; #SBATCH --mem=8G; #SBATCH --time=01:00:00; #SBATCH --array=0-23%6
- `slurm/reduce.sbatch`: #SBATCH --partition=p.hpcl93; #SBATCH --cpus-per-task=16; #SBATCH --mem=64G; #SBATCH --time=01:00:00

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| none found as top-level constants | | | |
