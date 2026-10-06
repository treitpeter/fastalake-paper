# capscan_region_mixed_model

Internal ID: A79 (previous folder name A79_capscan_region_mixed_model)

Does the FastaLake-to-metagenome peptide ratio differ by sampling site along the gut?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `region_test.py` | call chain of the used outputs |

Launchers: `slurm/run.sbatch`.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: linear mixed model of ln(FL/MG) by sampling site with participant random intercept (job 6057982) | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/Fig4. What is used: linear mixed model of ln(FL/MG) by sampling site with participant random intercept (job 6057982); ported into panels.py:stat_mixed_model and recomputed

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

| jobs | cores | requested memory | peak RSS | longest task | summed task time | states |
|---|---|---|---|---|---|---|
| 1 tasks | 2 | 4G | 0.0 GB | 0 min 05 s | 0 min 05 s | COMPLETED 1 |

Partitions: p.hpcl92. Job ids: 6057982.

## Choices

Approach: linear mixed model of ln(FL/MG) on sampling site with a participant random intercept.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| response | ln(FL/MG) per acquisition | ratio is multiplicative | participant-level table (PARTICIPANT_LEVEL.tsv) |
| random effect | participant intercept | acquisitions nest in participants | none |
| unit | acquisition nested in participant | per CONVENTIONS section 3 | none |

## Known limitations

See `OPEN_ITEMS.md`.
