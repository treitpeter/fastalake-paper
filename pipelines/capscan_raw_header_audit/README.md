# capscan_raw_header_audit

Internal ID: A89 (previous folder name A89_capscan_raw_header_audit)

What are the instrument-record file names inside the 368 CapScan RAW files?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `audit_headers.py` | call chain of the used outputs |

Launchers: none, the step ran interactively or inside another launcher.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: all_RAW_header_names.tsv (368 CapScan RAW files, instrument-record file names) for the specimen-number discrepancy census | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: tables/SData1. What is used: all_RAW_header_names.tsv (368 CapScan RAW files, instrument-record file names) for the specimen-number discrepancy census

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

No Slurm job ids are recorded for this step (see OPEN_ITEMS.md). It ran without a recorded accounting entry.

## Choices

Approach: all_RAW_header_names.tsv.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| scope | all 368 RAW files | specimen-number discrepancy census needs every file | none |

## Known limitations

See `OPEN_ITEMS.md`.
