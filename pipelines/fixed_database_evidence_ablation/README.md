# fixed_database_evidence_ablation

Internal ID: A61 (previous folder name A61_fixed_dictionary_evidence_ablation)

With the database held fixed, how much of the identification depends on each kind of evidence (de novo support partitions)?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `h01_eligibility.py` | call chain of the used outputs |
| `h09_evidence_ablation.py` | call chain of the used outputs |
| `partition.py` | call chain of the used outputs |
| `plot_corrected_ablation.py` | call chain of the used outputs |
| `prepare.py` | call chain of the used outputs |

Launchers: `slurm/H09_full.sbatch`, `slurm/H09_pilot.sbatch`, `slurm/full.sbatch`, `slurm/pilot.sbatch`, `slurm/prepare.sbatch`.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: acquisition_summary.tsv of PREDICT (296) and CapScan (282), combined into plot_values.tsv, SFig08 b-d | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/SFig08. What is used: acquisition_summary.tsv of PREDICT (296) and CapScan (282), combined into plot_values.tsv, SFig08 b-d

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

| jobs | cores | requested memory | peak RSS | longest task | summed task time | states |
|---|---|---|---|---|---|---|
| 7 tasks | 2 | 24G, 8G | 2.8 GB | 2 min 54 s | 11 min 35 s | COMPLETED 7 |

Partitions: p.hpcl92. Job ids: 5940770, 5940771, 5940772, 5928242, 5928243.

## Choices

Approach: per-acquisition ablation summaries for PREDICT (296) and CapScan (282), combined into plot_values.tsv.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| database | held fixed per acquisition | separates evidence from dictionary size | none |
| eligibility | h01_eligibility rule | frozen eligibility of the H09 analysis | none |
| units | acquisition | per-acquisition summaries | participant level in the figure |

## Known limitations

See `OPEN_ITEMS.md`.
