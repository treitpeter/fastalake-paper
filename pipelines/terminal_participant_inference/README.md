# Participant-level terminal inference and top-30 sensitivity (A24)

Internal ID: A24 (previous folder name A24_terminal_inference)

Pipeline id `terminal_participant_inference`. Used by: Fig6 (panels a to c); text_statistics supplementary_methods_8_11 (families terminal_fraction__length_9 and all18_nonKR_composition__length_9; statistics re-run in `panels.stat_terminal_families`).

## Question

Do capsule regions and stool differ in the terminal residues of their peptides, per participant, and does the de novo selection threshold matter? Answer: planned paired contrasts on the A20 counts and the A26 threshold counts, with every null kept. No search is repeated.

## Inputs

- A20 context, origin and route counts and A26 revision 2 threshold counts (manifest sha256 in `RECEIPTS.md`).
- Frozen specimen and participant mapping (CapScan 331 acquisitions; participant 15 added for the 347-acquisition run).

## Outputs

- `results/<cohort>/`: `acquisition_design.tsv`, `per_acquisition_endpoints.tsv.gz`, `per_specimen_endpoints.tsv.gz`, `per_participant_region_endpoints.tsv.gz`, `participant_paired_differences.tsv.gz`, `region_summary.tsv`, `complete_statistical_families.tsv`, `statistical_family_summary.tsv`, `threshold_pairs_*.tsv.gz`.

## Call chain and code

One script, `search.py`, is both the original and the participant-15 `search.py` (byte-identical). Run-directory builders for the 331 and 347 runs: `terminal_route_counts/prep_chain.py`.

| file | origin (LAB_ROOT relative) | sha256 at origin | sha256 here |
|---|---|---|---|
| `search.py` | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/frozen/search.py` | `c980eccbf1bc2989a7ea1ac4e3b7eb3baa3446eb33d5f3e869b1a54a19f34a8e` | `c980eccbf1bc2989a7ea1ac4e3b7eb3baa3446eb33d5f3e869b1a54a19f34a8e` |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| endpoints | non-K/R, all, H/all and H/non-K/R for N-first, actual preceding P1 and C-last; all 18 non-K/R residues as bare-sequence composition at 9 to 50 | predeclared | residue selection post hoc is not done |
| primary origins | All, Human, Microbial; other origin classes descriptive | conservative classes | mixed and unknown reported only |
| unit | participant (acquisitions averaged within specimen, specimens within participant and region); coarse compartments need all constituent regions; PREDICT averages specimens per patient | no pseudo-replication | row units reported separately from canonical units |
| contrasts | ten region pairs and four predefined compartment contrasts; threshold contrasts form separate cohort families | predeclared | none |
| test | exact sign-flip test for at most 16 participants, otherwise 200000 sign-flip draws with plus-one correction; at least 8 paired participants, otherwise p = 1 (non-estimable) | small paired samples | assumes sign symmetry; exploratory |
| bootstrap | participant bootstrap, 50000 draws for compact endpoints and 2000 for residue composition | interval for the mean paired difference | pointwise intervals |
| multiplicity | Benjamini-Hochberg over all tested routes, origins, termini and contrasts, separately for compact endpoints versus composition and for each length | one family per kind | none |
| random streams | seeded from sha256 of the stream key and 20260914 (`seed()` in `search.py`) | reproducible draws | none |
| resources | 1 CPU, 32 GB, 4 h, one task at a time | measured peak 2.1 GB | none |

## Resources and run time (measured)

De novo collection (job 5929071): 2.7 min, 2.0 GB; CapScan collection (5929072): 6.8 min, 2.1 GB; PREDICT collection of the first wave (5929073) was cancelled and replaced by the recovery job 5937596 (32 GB, 5.0 min, 2.0 GB). Participant-15 baseline (331) and fixed (347) runs: jobs 6058147 (6.8 min, 2.2 GB) and 6058149 (7.3 min, 2.2 GB), 3 GB requested.

Per-job accounting is in `RECEIPTS.md`.

## How to run

Every script takes a run directory as first argument and, for array stages, a task index as second argument. The run directory holds `inputs/manifest*.json` (task list with the path, size, modification time and sha256 of every input), `frozen/` (this pipeline's code, assembled by `capscan_nonspecific_search/slurm/assemble_run_dir.sh`) and `results/`.

    export FASTALAKE_PIPELINES_DIR=<repo>/pipelines
    export FASTALAKE_OUT_DIR=<run directory>      # contains inputs/manifest*.json
    export FASTALAKE_PYTHON=<python with the packages named below>
    sbatch --array=<first>-<last>%<concurrency> <flags from slurm/submit_example.sh> slurm/run_task.sbatch

The manifests list original instrument-derived inputs by absolute path and are not part of the repository (they name lab locations; their sha256 is in `RECEIPTS.md`). A re-run therefore needs a manifest rebuilt for the new location; the frozen code verifies every input against the manifest (path, size, modification time, sha256) and stops on any difference. Per the repository conventions this pipeline is not re-run for the paper; the code is the code that produced the outputs.

## Receipts and open items

`RECEIPTS.md` lists the original job ids and the sha256 of every output. `OPEN_ITEMS.md` lists what could not be traced.
