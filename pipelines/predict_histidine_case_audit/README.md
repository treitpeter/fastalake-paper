# Independent audit of the four PREDICT spectrum cases (A25 audit)

Internal ID: A25 (previous folder name A25_predict_spectrum_case_audit)

Pipeline id `predict_histidine_case_audit`. Used by: SFig14.

## Question

Are the four PREDICT histidine-boundary cases of SFig14 correctly selected and supported by their measured spectra? Answer: an independent reconstruction of selection, indexed RAW peaks, Sage preprocessing, local source loci and float64 b/y fragment matching against all 306 acquisition receipts; all four cases pass (13, 8, 11 and 7 matched Sage ions; backbone coverage 81.8, 60.0, 75.0 and 75.0 percent; peptide q and spectrum q <= 0.01).

## Inputs

- Receipts and outputs of the corrected PREDICT spectrum selection (case table, candidate counts for 4 x 306 strata, eligible selection PSMs, same-scan calls; `MANIFEST.json`, sha256 in `RECEIPTS.md`), the mzML of each case, each case FASTA and non-specific configuration, the baseline tryptic peptide sets.

## Outputs

- `full/independent_case_review.tsv` (one row per case), per-case plots and quantitative files; the pilot phase covers cases HIS_002 and HIS_004 and must pass before the full phase.

## Call chain and code

One script. It reads receipts through `MM_REVIEW_ROOT` (renamed `FASTALAKE_DATA_DIR` at run time in the launcher).

| file | origin (LAB_ROOT relative) | sha256 at origin | sha256 here |
|---|---|---|---|
| `audit_predict_spectrum_cases.py` | `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/predict_spectrum_review/frozen/audit_predict_spectrum_cases.py` | `ae0edcea5f876f4394d50ab36a775769fad0863759fa99712d8a9176bc1d4320` | `ae0edcea5f876f4394d50ab36a775769fad0863759fa99712d8a9176bc1d4320` |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| case rule re-check | median candidate-count acquisition of the stratum (index len//2 of the sorted candidates), median matched-intensity eligible peptide, rank-1 target, peptide q and spectrum q <= 0.01, 9 to 50 residues, peptide absent from the tryptic baseline | the same rule as A25, re-derived independently | none |
| fragment tolerance | +/- 20 ppm on b and y ions; precursor offset 0, -1 and +1 checked | as the search | offset checks are descriptive |
| phases | pilot (2 cases) then full (4 cases); the full phase requires the pilot receipt and the identical worker hash | guards against a changed worker | none |
| resources | 2 CPUs allocated, 16 GB | measured under 1 min | none |

## Resources and run time (measured)

Pilot job 5940795 (0.3 min) and full job 5940796 (0.5 min), 16 GB requested, 2 CPUs allocated; peak memory not recorded in the accounting rows (below the sampling resolution).

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
