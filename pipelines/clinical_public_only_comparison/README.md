# clinical_public_only_comparison

Internal ID: A16 (previous folder name A16_public_comparison)

Per acquisition, how many peptides does the public-only FastaLake search identify compared with the matched metagenome search on the same spectra, and how many metagenome peptides does it recover?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `lib.py` | call chain of the used outputs |
| `report.py` | call chain of the used outputs |
| `search.py` | call chain of the used outputs |
| `summarize_ordinary.py` | call chain of the used outputs |

Launchers: `slurm/pilot.sbatch`, `slurm/report.sbatch`, `slurm/summary.sbatch`.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: same-spectrum FL_public vs MG per-acquisition peptide counts and identity recovery, 296 PREDICT and 282 CapScan acquisitions (Source Data Figure4 benchmark_{0,1 | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/Fig4, figures/SFig06, text_statistics/supplementary_methods_6_7. What is used: same-spectrum FL_public vs MG per-acquisition peptide counts and identity recovery, 296 PREDICT and 282 CapScan acquisitions (Source Data Figure4 benchmark_{0,1}_per_acquisition.tsv; SFig06; confidence_ordinary_all_pairs.tsv, lines 633, 777)

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

| jobs | cores | requested memory | peak RSS | longest task | summed task time | states |
|---|---|---|---|---|---|---|
| 5 tasks | 2 | 8G | 0.6 GB | 1 min 40 s | 3 min 23 s | COMPLETED 5 |

Partitions: p.hpcl92. Job ids: 5927617, 5927618, 5928360.

## Choices

Approach: same-spectrum counts (peptide_q01, 5-50 residues) and identity recovery for 296 PREDICT and 282 CapScan acquisitions.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| peptide count | accepted peptides, 5-50 residues, q <= 0.01 | production acceptance rule | 9-50 also reported |
| pairing | same acquisition, same spectra | isolates the database | none |
| unit | acquisition | statistics are reported per acquisition and per participant | participant medians in the figure tables |

## Known limitations

See `OPEN_ITEMS.md`.
