# capscan_dominance_evenness

Internal ID: A83 (previous folder name A83_capscan_dominance_evenness_347)

How dominant is the leading species and how even is the community, per acquisition and region?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `dominance_evenness.py` | call chain of the used outputs |
| `dominance_evenness_347.py` | call chain of the used outputs |

Launchers: none, the step ran interactively or inside another launcher.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: per-acquisition leading-species share and Pielou evenness and the paired Wilcoxon tests | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/Fig5, text_statistics/supplementary_methods_8_11. What is used: per-acquisition leading-species share and Pielou evenness and the paired Wilcoxon tests; fully ported into panels.py (estimator mismatch recorded in OPEN_ITEMS.md)

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

| jobs | cores | requested memory | peak RSS | longest task | summed task time | states |
|---|---|---|---|---|---|---|
| 1 tasks | 2 | 16G | 0.0 GB | 0 min 06 s | 0 min 06 s | COMPLETED 1 |

Partitions: p.hpcl92. Job ids: 6063821.

## Choices

Approach: leading-species share and Pielou evenness, paired Wilcoxon tests.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| evenness | Pielou | standard | none |
| test | paired Wilcoxon on participant medians | paired by participant | none |

## Known limitations

See `OPEN_ITEMS.md`.
