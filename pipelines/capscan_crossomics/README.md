# capscan_crossomics

Internal ID: A82 (previous folder name A82_capscan_crossomics_347)

How do protein-based, DNA-based and metabolite profiles agree across regions in the 13 participants with complete DNA?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `build_groups_capscan.py` | call chain of the used outputs |
| `crossomics_capscan.py` | call chain of the used outputs |
| `crossomics_src/check_prevalence.py` | call chain of the used outputs |
| `crossomics_src/native_pyramids.py` | call chain of the used outputs |

Launchers: `slurm/crossomics.sbatch`.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: CapScan native_correlations.tsv, functional_correlations.tsv and distance matrices (native_distances.npz, functional_distances.npz) at 347 acquisitions, 13 DNA- | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/Fig5, figures/SFig11, tables/SData1, text_statistics/supplementary_methods_8_11. What is used: CapScan native_correlations.tsv, functional_correlations.tsv and distance matrices (native_distances.npz, functional_distances.npz) at 347 acquisitions, 13 DNA-complete participants (Slurm 6063849): Fig5 g, SFig11 b,d; eligibility_audit.tsv (346 of 347) for SData1 rows

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

| jobs | cores | requested memory | peak RSS | longest task | summed task time | states |
|---|---|---|---|---|---|---|
| 4 tasks | 2 | 2800M, 2G, 4G | 0.4 GB | 1 min 19 s | 1 min 56 s | FAILED 1, COMPLETED 2, OUT_OF_MEMORY 1 |

Partitions: p.hpcl92. Job ids: 6063599, 6063836, 6063837, 6063849.

## Choices

Approach: native and functional correlations and distance matrices at 347 acquisitions.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| seed | 20260907 | fixed for reproducibility | none |
| permutations | 1999 | production setting of crossomics.py | none |
| participants | 13 DNA-complete | require protein, metabolite and DNA cells in all five groups | 12 in the earlier design |
| multiple testing | Benjamini-Hochberg per cohort/scope/prevalence family | as production | none |

## Known limitations

See `OPEN_ITEMS.md`.
