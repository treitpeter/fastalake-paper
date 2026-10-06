# metagenome_split_spectra_comparison

Internal ID: A96 (previous folder name A96_split_spectra_mg_comparison)

Does the gain over matched metagenome searches hold when only odd-scan spectra are searched?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `analyze.py` | call chain of the used outputs |
| `worker.py` | call chain of the used outputs |

Launchers: none, the step ran interactively or inside another launcher.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: Per-acquisition counts and gains of the matched MG searches on odd-scan spectra (Source Data package sheet) and the 10,000-draw participant bootstrap | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: text_statistics/methods. What is used: Per-acquisition counts and gains of the matched MG searches on odd-scan spectra (Source Data package sheet) and the 10,000-draw participant bootstrap

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

| jobs | cores | requested memory | peak RSS | longest task | summed task time | states |
|---|---|---|---|---|---|---|
| 512 tasks | 16, 32 | 16G | 16.2 GB | 6 min 48 s | 29 h 52 min | COMPLETED 510, OUT_OF_MEMORY 2 |

Partitions: p.hpcl92, p.hpcl93. Job ids: 6048981, 6065298.

## Choices

Approach: per-acquisition counts and gains, participant bootstrap.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| bootstrap | 10,000 draws, seed 20261003, percentile 95% CI | participant-level resampling keeps acquisitions of a participant together | none |
| split | odd scans | disjoint from the evidence spectra | 10-fold in A95 |

## Known limitations

See `OPEN_ITEMS.md`.
