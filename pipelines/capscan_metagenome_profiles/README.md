# capscan_metagenome_profiles

Internal ID: A109 (previous folder name A109_capscan_metagenome_aggregate)

What is the DNA-based taxonomic and functional profile of the 301 CapScan runs with public paired WGS?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `aggregate.py` | call chain of the used outputs |
| `prepare.py` | call chain of the used outputs |
| `run_sample.py` | call chain of the used outputs |
| `run_selection.py` | call chain of the used outputs |
| `submit.py` | call chain of the used outputs |

Launchers: `slurm/environment.sh`, `slurm/profile.sbatch`.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: 301 runs of the CapScan functional DNA aggregate (capscan_dna_runs.tsv) | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: text_statistics/supplementary_methods_8_11. What is used: 301 runs of the CapScan functional DNA aggregate (capscan_dna_runs.tsv)

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

| jobs | cores | requested memory | peak RSS | longest task | summed task time | states |
|---|---|---|---|---|---|---|
| 351 tasks | 16, 2, 8 | 128M, 32G, 48G | 18.7 GB | 7 h 11 min | 1057 h 57 min | COMPLETED 205, TIMEOUT 145, FAILED 1 |

Partitions: p.hpcl8, p.hpcl92. Job ids: 5907383, 5907384, 5949070.

## Choices

Approach: 301-run DNA aggregate (mOTUs 3.0.3 taxonomy, MEGAHIT 1.2.9 assembly, Prodigal 2.6.3, Salmon 1.10.0, eggNOG-mapper 2.1.12).

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| assembler | MEGAHIT 1.2.9 defaults | per-sample assembly | none |
| quantification | Salmon 1.10.0, IU, --meta | relative DNA gene abundance | none |
| annotation | eggNOG-mapper 2.1.12, eggNOG 5.0.2, DIAMOND 2.0.11, prokaryota_broad | as the protein arm | none |
| exclusion | SRR18585253 excluded | recorded in the original selection | original 302 runs frozen |
| time limit | 4 h per task | 145 tasks timed out; recovery runs used longer limits | see OPEN_ITEMS |

## Known limitations

See `OPEN_ITEMS.md`.
