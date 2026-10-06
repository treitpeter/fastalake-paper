# predict_curated_mg_comparator

Internal ID: A60 (previous folder name A60_curated_specimen_mg_comparator)

Does curating the specimen metagenome database change the peptide counts relative to FastaLake, in the 296 PREDICT acquisitions?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `collect_parallel.py` | call chain of the used outputs |
| `search.py` | call chain of the used outputs |

Launchers: `slurm/collect_full.sbatch`, `slurm/collect_parallel_full.sbatch`, `slurm/collect_parallel_pilot.sbatch`, `slurm/collect_pilot.sbatch`, `slurm/search.sbatch`.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: curation_counts.tsv (accepted 5-50-residue peptides per acquisition in FL_cohortMG, MG, FL_public and MG_own_curated arms, 296 PREDICT acquisitions), SFig08 a | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/SFig08. What is used: curation_counts.tsv (accepted 5-50-residue peptides per acquisition in FL_cohortMG, MG, FL_public and MG_own_curated arms, 296 PREDICT acquisitions), SFig08 a

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

| jobs | cores | requested memory | peak RSS | longest task | summed task time | states |
|---|---|---|---|---|---|---|
| 2 tasks | 2, 4 | 48G, 64G | 1.2 GB | 15 min 48 s | 23 min 04 s | CANCELLED 1, COMPLETED 1 |

Partitions: p.hpcl8. Job ids: 5942784, 5942965.

## Choices

Approach: accepted peptides (5-50 residues) in four arms: FL_cohortMG, MG, FL_public and MG_own_curated.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| peptide window | 5-50 residues | production rule | none |
| arms | cohort-metagenome FastaLake, specimen metagenome, public-only FastaLake, curated own metagenome | the four comparators of Supplementary Figure 8a | none |
| replay | five affected curated-own-metagenome searches replayed from archived FASTAs | their original work directories had been deleted by the archived job script | limitation kept in the original PROTOCOL; see OPEN_ITEMS |

## Known limitations

See `OPEN_ITEMS.md`.
