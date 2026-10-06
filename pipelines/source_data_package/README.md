# source_data_package

Internal ID: A88 (previous folder name A88_source_data_package)

How is the Source Data package (STab4, STab5, SData1) built from the verified tables, with pseudonymised identifiers?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `build_package.py` | call chain of the used outputs |
| `derive_verified.py` | call chain of the used outputs |
| `specials.py` | call chain of the used outputs |
| `tables_cr003.py` | call chain of the used outputs |

Launchers: none, the step ran interactively or inside another launcher.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: Merged STab4_workbook_build and STab5_workbook_build (same Source Data package builder build_package.py/specials.py/tables_cr003.py). PREDICT pseudonyms (PREDIC | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/Fig4, tables/SData1, tables/STab4, tables/STab5. What is used: Merged STab4_workbook_build and STab5_workbook_build (same Source Data package builder build_package.py/specials.py/tables_cr003.py). PREDICT pseudonyms (PREDICT_P###) for Fig4 participant_id; pseudonymisation and path cleaning of SData1, specimen-number discrepancy census (derive_verified.py); STab4 and STab5 workbook assembly, Read_me text, column/label renames; logic ported to panels.build_sheets

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

No Slurm job ids are recorded for this step (see OPEN_ITEMS.md). It ran without a recorded accounting entry.

## Choices

Approach: Source Data package builder, pseudonymisation, path cleaning and specimen-number discrepancy census.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| pseudonyms | PREDICT_P### for PREDICT participants | no raw identifiers released | none |
| path cleaning | regular expressions replace lab and home paths | released tables must carry no lab paths | none |

## Known limitations

See `OPEN_ITEMS.md`.
