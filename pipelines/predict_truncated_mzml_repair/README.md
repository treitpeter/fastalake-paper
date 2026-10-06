# Re-search of the truncated PREDICT acquisitions from reconverted RAW (A25 repair)

Internal ID: A25 (previous folder name A25_repaired_predict_search)

Pipeline id `predict_truncated_mzml_repair`. Used by: SFig14.

## Question

Five PREDICT acquisitions had truncated mzML files. Can they be recovered from the original RAW files and re-searched so that the PREDICT results use complete data? Answer: each corrupt mzML was reconverted from its hash-identified RAW file, gated on full scan counts and an intact-spectrum comparison, then searched twice (tryptic with the production enzyme, non-specific with the A13 configuration) with the unchanged Sage worker.

## Inputs

- The five original RAW files (hash-identified in `raw_recovery/manifest.json`) and the historical mzML, the A13 manifest and frozen worker, the original receipts (`raw_recovery/manifest.json` sha256 in `RECEIPTS.md`).
- Sage 0.14.6 binary, sha256 `d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625`.

## Outputs

- `raw_recovery/<sample>/{<sample>.mzML,COMPLETE.json,CONVERSION.json}`: reconverted mzML and its validation receipt (`ALL_THREE_VALIDATED.json` covers the three samples of the first manifest).
- `search_recovery/<sample>/{tryptic,unspecific}/results/<sample>__<arm>/`: complete Sage results and `COMPLETE.json`.

## Call chain and code

The search script calls `capscan_nonspecific_search/search.py` (identical to the A13 worker) with a run directory it builds under `search_recovery/`. The collector that re-ran A20, A24, A25 and A28 on these acquisitions is `terminal_histidine_inference/collect_repaired_predict_20260915.py`.

| file | origin (LAB_ROOT relative) | sha256 at origin | sha256 here |
|---|---|---|---|
| `repaired_predict_searches_20260915.py` | `LAB_ROOT/MM_review_analysis/code/repaired_predict_searches_20260915.py` | `b31b2e94fae08ed9a46173917c412e04a8711444edfd9595188682e2fdbba1ef` | `b31b2e94fae08ed9a46173917c412e04a8711444edfd9595188682e2fdbba1ef` |
| `recover_predict_raw_20260915.py` | `LAB_ROOT/MM_review_analysis/code/recover_predict_raw_20260915.py` | `906329f707cdf9c5e9bc1cf45276e7751b41848aa39d20a512ef37c6728b5047` | `906329f707cdf9c5e9bc1cf45276e7751b41848aa39d20a512ef37c6728b5047` |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| Sage engine | Sage 0.14.6 binary, sha256 `d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625`, version verified from the program output before each run | one engine for all searches of the review analyses | other engines: campi_second_engine |
| reconversion gate | new mzML is reused only after full scan-count agreement and a prespecified intact-spectrum comparison against the historical vendor-centroided data (`recover_predict_raw_20260915.py`) | a damaged input must not change results silently | historical files are never altered |
| tryptic arm | enzyme block copied from the original production configuration (`cleave_at` KR, not semi-enzymatic) | a matched baseline on the repaired data | original search kept as the historical comparison |
| non-specific arm | A13 configuration with only the raw path and output directory changed; must show at least one fully non-tryptic compatible peptide | same experiment as A13 | none |
| search code | `search.py` of A05/A13 run unchanged (called with the new run directory and task 0) | no new code path | none |
| resources | 8 threads, 480 GB, array of two tasks | measured peak 257.8 GB, 36.5 min maximum | none |

## Resources and run time (measured)

Search array job 5937565 (2 tasks), 8 threads (16 allocated hardware threads), 480 GB requested, peak 257.8 GB, 19.1 and 36.5 min. The reconversion jobs ran before (converter and gate: `recover_predict_raw_20260915.py`, launcher `convert.sbatch`, a one-line array wrapper that only names lab paths and is not copied).

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
