# De novo score >= 0.9 versus production top 30 percent (A26)

Internal ID: A26 (previous folder name A26_denovo_thresholds)

Pipeline id `denovo_threshold_routes`. Used by: Fig6 (panel a, top-30 percent route).

## Question

Do the historical high-confidence de novo route (score >= 0.9) and the production route (top 30 percent after eligibility) give different terminal counts, and is the difference the threshold or the length and complexity filter? Answer: both routes are counted for every acquisition and kept separate, with high-confidence 5 to 50 and 9 to 50 strata that isolate the threshold on a common length window. No search is repeated.

## Inputs

- De novo prediction tables of every acquisition (CapScan 331, PREDICT 306) and each acquisition searched FASTA, as listed in `inputs/manifest_revision2.json` (sha256 in `RECEIPTS.md`).

## Outputs

- `results_revision2/<cohort>/<acquisition>/` counts per route, residues, fractional parent P1 context and origin ambiguity classes, with `COMPLETE.json`.

## Call chain and code

`context_helpers.py` is shared with A25, A27 and A28 (one copy kept here). `search.py` is the first-round worker, `search_revision2.py` the worker of the results used; `origin.py` is in A20.

| file | origin (LAB_ROOT relative) | sha256 at origin | sha256 here |
|---|---|---|---|
| `search.py` | `LAB_ROOT/MM_review_analysis/A26_4_7_denovo_thresholds/frozen/search.py` | `8dce0be3bdd74fbaa5cdf1765bcf6e12083694691b6e1ce49a7cd26ccbc7aac5` | `8dce0be3bdd74fbaa5cdf1765bcf6e12083694691b6e1ce49a7cd26ccbc7aac5` |
| `search_revision2.py` | `LAB_ROOT/MM_review_analysis/A26_4_7_denovo_thresholds/frozen/search_revision2.py` | `30b46e73cad086fa9b9525f7285f3603dc050601f2f1cf5e0d8e19fb5cbc783c` | `30b46e73cad086fa9b9525f7285f3603dc050601f2f1cf5e0d8e19fb5cbc783c` |
| `historical_denovo.py` | `LAB_ROOT/MM_review_analysis/A26_4_7_denovo_thresholds/frozen/historical_denovo.py` | `a7adc4b0a0b6f27b382efeae9aca38913a2374ae9c2f84ab513cd14c7aa3eb2c` | `0032978983f2d11a54fc337e9182fc95a75d5d8028000e0525a915d312e6134d` |
| `context_helpers.py` | `LAB_ROOT/MM_review_analysis/A26_4_7_denovo_thresholds/frozen/context_helpers.py` | `1a4d4a8746e4ddc30fd25be5c224c829bb752f7710b9bbd3a82d0032ab0a7038` | `1174aea1f9d264abb45080687167e1db5ae89f931a886c185ac4793ebe7db662` |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| historical route | score >= 0.9, all non-empty positive-score exported predictions, all lengths (`historical_denovo.py`) | reproduces the manuscript counts; every CapScan acquisition must reproduce its original counts | kept as a separate route |
| production route | top 30 percent of positive-score, 9 to 50, non-low-complexity predictions, including threshold ties; cutoff = ceil(0.3 n)-th highest score | exact production selector (`lib.load_queries`) | none |
| added strata (revision 2) | high-confidence 5 to 50 and 9 to 50 | separates score selection from length and complexity | requested to isolate the top-30 effect |
| mapping scope | context mapping restricted to 5 to 50 residues; longer predictions stay in historical counts and are labelled outside mapping scope | not absent from the FASTA | none |
| origin labels | sequence-compatible source, not unique molecular origin | ambiguity kept explicit | see A20 |
| pilots | three: one CapScan, the zero-quantification CapScan acquisition, one PREDICT | empty prediction routes keep explicit zero or undefined counts | none |
| resources | 1 CPU, 16 GB, 1 h, one task at a time | measured peak below 1 GB | none |

## Resources and run time (measured)

Revision 2 pilots (job 5928961, 3 tasks) and full array (job 5928962, 634 tasks): 1 CPU, 16 GB requested, per-task wall time median 0.2 min, maximum 0.4 min. The first-round full array (5928927) was cancelled after 14 tasks when revision 2 was requested. Participant-15 stage (job 6058101, 16 tasks): 2 GB, maximum 0.5 min.

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
