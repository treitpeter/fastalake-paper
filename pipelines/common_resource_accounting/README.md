# common_resource_accounting

Internal ID: A72 (previous folder name A72_compute_resource_accounting)

What did the FastaLake stages cost in CPU hours and memory, from the Slurm accounting of the original runs?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `build_revision.py` | call chain of the used outputs |
| `collect.py` | call chain of the used outputs |
| `compute_block.py` | call chain of the used outputs |
| `h10_compute_provenance_v2.py` | call chain of the used outputs |
| `original_collect_sacct.py` | call chain of the used outputs |

Launchers: none, the step ran interactively or inside another launcher.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: Slurm accounting (sacct) of the production chain restricted by array index to manuscript acquisitions, incl. CAMPI chain and AlphaNovo inference: wall times per | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/SFig15, text_statistics/abstract_discussion, text_statistics/methods. What is used: Slurm accounting (sacct) of the production chain restricted by array index to manuscript acquisitions, incl. CAMPI chain and AlphaNovo inference: wall times per stage and task (SFig15, line 193, Methods line 299)

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

| jobs | cores | requested memory | peak RSS | longest task | summed task time | states |
|---|---|---|---|---|---|---|
| 2 tasks | 2 | 120M, 200M | 0.0 GB | 0 min 01 s | 0 min 02 s | COMPLETED 2 |

Partitions: p.hpcl92. Job ids: 6064882, 6064883.

## Choices

Approach: per-task Slurm records and per-stage summaries for the baseline 331 and the fixed 347 CapScan runs.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| source | sacct records of the original jobs | accounting is measured, not estimated | none |
| blocks | baseline331 and fixed347 compute blocks | the 331-to-347 rerun changed the roster | both reported |

## Known limitations

See `OPEN_ITEMS.md`.
