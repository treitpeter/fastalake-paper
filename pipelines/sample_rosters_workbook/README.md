# sample_rosters_workbook

Internal ID: A87 (previous folder name A87_sample_rosters_workbook)

Which acquisitions, specimens and participants enter each cohort and analysis?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `build_histidine_scope.py` | call chain of the used outputs |
| `build_rosters.py` | call chain of the used outputs |
| `build_workbook.py` | call chain of the used outputs |

Launchers: none, the step ran interactively or inside another launcher.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: roster tables and workbook layout (all 34 sheets of the v28d sample-roster workbook) | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: tables/SData1. What is used: roster tables and workbook layout (all 34 sheets of the v28d sample-roster workbook); sheet writer and final checks ported into panels.py

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

No Slurm job ids are recorded for this step (see OPEN_ITEMS.md). It ran without a recorded accounting entry.

## Choices

Approach: the 34-sheet sample-roster workbook.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| sheets | 34 sheets | one per roster and scope | none |

## Known limitations

See `OPEN_ITEMS.md`.
