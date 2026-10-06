# peptide_centric_functional_allocation

Internal ID: A33 (previous folder name A33_peptide_centric_functional_allocation)

How is each acquisition's signal allocated over KOs, and how much remains FunctionallyUnassigned, for HSTOOL, repaired PREDICT and the 347-acquisition CapScan rerun?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `build_variants.py` | call chain of the used outputs |

Launchers: `slurm/repaired_predict.sbatch`, `slurm/run.sh`, `slurm/stage.sh`.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: KO x acquisition allocation matrices with the FunctionallyUnassigned remainder (HSTOOL, repaired PREDICT, CapScan rerun on the 347 roster, Slurm 5973734 and lat | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/SFig12. What is used: KO x acquisition allocation matrices with the FunctionallyUnassigned remainder (HSTOOL, repaired PREDICT, CapScan rerun on the 347 roster, Slurm 5973734 and later)

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

| jobs | cores | requested memory | peak RSS | longest task | summed task time | states |
|---|---|---|---|---|---|---|
| 3 tasks | 2 | 10G, 5000M, 6G | 4.4 GB | 5 min 22 s | 10 min 38 s | COMPLETED 3 |

Partitions: p.hpcl92. Job ids: 6065988, 6066031, 6066033.

## Choices

Approach: KO x acquisition allocation matrices with the FunctionallyUnassigned remainder.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| remainder | kept as an explicit FunctionallyUnassigned row | allocation must sum to the whole signal | none |
| CapScan roster | 347 acquisitions | fixed by A80 | 331 reproduction run (job 6065988) checks the producer |

## Known limitations

See `OPEN_ITEMS.md`.
