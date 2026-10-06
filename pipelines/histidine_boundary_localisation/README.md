# Histidine boundary localisation and source expectations (A27)

Internal ID: A27 (previous folder name A27_histidine_defence)

Pipeline id `histidine_boundary_localisation`. Used by: Fig6 (panels d to f); SFig14 (site_evidence.tsv.gz, route sage_unspecific_bothq, defines the candidate peptides).

## Question

Where in the source protein does each observed internal cleavage boundary fall, and what P1 residue distribution would the same protein predict? Answer: every accepted sequence occurrence is mapped onto the acquisition FASTA, per-site evidence and source-matched expected distributions are saved for all 331 CapScan and 306 PREDICT acquisitions in two phases (available routes, then fully non-specific). Statistics are in A28.

## Inputs

- Existing tryptic, semi-enzymatic and non-specific accepted peptide sets (A05, A13), de novo tables, each acquisition FASTA (`inputs/manifest.json`, `inputs/predict_submission_manifest.json`; sha256 in `RECEIPTS.md`).

## Outputs

- `results/{available,unspecific}/<cohort>/<acquisition>/`: `site_evidence.tsv.gz`, source-expectation and complementary-product tables, `COMPLETE.json`.

## Call chain and code

`histidine_helpers.py` here differs from the A28 copy (different sha256, see table); both are kept because each was frozen with its stage.

| file | origin (LAB_ROOT relative) | sha256 at origin | sha256 here |
|---|---|---|---|
| `search.py` | `LAB_ROOT/MM_review_analysis/A27_4_7_histidine_defence/frozen/search.py` | `801337c61cd4ca5682a22501ea2d93a5cabaf21712e0b0af405363b17267fae5` | `801337c61cd4ca5682a22501ea2d93a5cabaf21712e0b0af405363b17267fae5` |
| `histidine_helpers.py` | `LAB_ROOT/MM_review_analysis/A27_4_7_histidine_defence/frozen/histidine_helpers.py` | `42c85103d2f0185ec096604e5eb5a55e63cd2d3acc236b85ddff6c3b7cf97521` | `42c85103d2f0185ec096604e5eb5a55e63cd2d3acc236b85ddff6c3b7cf97521` |
| `predict_index_wrapper.py` | `LAB_ROOT/MM_review_analysis/A27_4_7_histidine_defence/frozen/predict_index_wrapper.py` | `0ef337badf51716a511be51ceaed8ff939d65d8cf5e2ad39d5743fd7f4944f74` | `366868b6215e9ed1ff5c16219eaf65847575dce75f58e31973ff974731814607` |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| window | 9 to 50 residues, canonical I/L-collapsed identities | as the de novo comparison | 5 to 50 in A20 |
| de novo bins | production top 30 percent; score >= 0.9 with and without the complexity gate; six score bins (0,.5), [.5,.7), [.7,.9), [.9,.95), [.95,.99), [.99,inf) | retained so route and score effects stay separable | prediction scores are not q-values |
| Sage gate | target rank 1 with peptide q <= 0.01; nested spectrum q <= 0.01 | as A05 | none |
| compatible parents | recorded enzyme and P2 restriction, at most 2 missed cleavages, no initiator-M clipping; protein ends explicit and excluded from internal-site claims | internal-site claims need an internal site | ambiguous sites kept as descriptive sensitivities |
| primary unit | distinct internal sites supported by peptides with a single compatible full-protein locus; identical full proteins and I/L variants collapsed by sequence hash | avoids multiple counting of one site | all-ambiguous sensitivity saved |
| expected distribution | for each observed non-K/R internal site, the eligible internal non-K/R positions of the same source protein (needs a 9 to 50 peptide on that side, the route enzyme rule, mass 500 to 5000, C +57.0215, M +15.9949 up to 2) | controls source composition, not ionisation or detectability | all-internal-position background saved; no generic proteome background |
| validation | vectorised opportunity calculation checked by brute-force enumeration on real pilot proteins (`brute_validate`) | correctness check | none |
| resources | 1 CPU, 32 GB, 2 h | measured peak 5.6 GB (available phase), 2.4 GB (non-specific phase) | none |

## Resources and run time (measured)

Available phase: pilot job 5929425 (4 tasks, max 2.1 min, 4.6 GB), full job 5929426 (633 tasks, median 1.1 min, max 3.0 min, 5.6 GB). Non-specific phase: pilot 5929427 (4 tasks), CapScan 5929428 (328 tasks, max 1.4 min, 2.4 GB), PREDICT after an array-index repair 5929446 (302 tasks, max 1.6 min, 2.4 GB). All 1 CPU and 32 GB requested. Participant-15 stage: jobs 6058102 (available, 6 GB requested, max 2.9 min, 3.6 GB) and 6058104 (non-specific, max 1.1 min, 1.6 GB).

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
