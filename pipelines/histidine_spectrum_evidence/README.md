# Measured-spectrum evidence for histidine boundary peptides (A25)

Internal ID: A25 (previous folder name A25_histidine_spectra)

Pipeline id `histidine_spectrum_evidence`. Used by: SFig14; text_statistics supplementary_methods_8_11.

## Question

Do the measured spectra of histidine-boundary peptides carry fragment support on the cleavage side? Answer: one PSM per origin and boundary side is selected in a rule fixed before spectrum inspection, and its observed and processed peaks and matched b and y fragments are saved (stratum table `spectra/all_strata.tsv`).

## Inputs

- A27 non-specific site evidence, the acquisition mzML files and Sage results (`inputs/manifest_revision2.json`; sha256 in `RECEIPTS.md`).

## Outputs

- `results_revision2/{pilot,CAPSCAN,PREDICT}/all_strata.tsv`, per-stratum selected PSM, observed and processed peaks, matched fragment annotations, `COMPLETE.json`.

## Call chain and code

`lib.py`, `peptide.py`, `context_helpers.py` are the shared copies (A05, A26); the A25 copies are byte-identical. The `case audit` and the repaired search are the sibling pipelines predict_histidine_case_audit and predict_truncated_mzml_repair.

| file | origin (LAB_ROOT relative) | sha256 at origin | sha256 here |
|---|---|---|---|
| `search_revision2.py` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/frozen/search_revision2.py` | `e8593f5b24caccaf7f6affb9126a3ac6f3ec0825b88de7250968cea75db3b7a2` | `e8593f5b24caccaf7f6affb9126a3ac6f3ec0825b88de7250968cea75db3b7a2` |
| `fragment_helpers.py` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/frozen/fragment_helpers.py` | `d0dca418ddbd95084a920ff5928d06dacd76341cd9530c8f9f33ed8358c4cd4e` | `2d8db5e3e9daba24c8b61d076e2d01028a72eb974667c31bbdff2872bf0782d0` |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| strata | 12: CapScan capsule and stool x Human and Microbial x N-P1H and C-terminal H, plus PREDICT stool x origin x side; empty strata kept | fixed before inspection | none |
| eligibility | 9 to 50 residues, absent from the same-acquisition tryptic search, unique internal source locus, target rank 1, peptide q <= 0.01 and spectrum q <= 0.01 | strict set | none |
| acquisition choice | median positive candidate-count acquisition per stratum; deterministic nearest-rank fallback only if its strict PSM set is empty | selection independent of the spectrum | none |
| PSM choice | per canonical peptide the best peptide q, then hyperscore and scan; across peptides the median matched-intensity peptide | no cherry-picking | all selection denominators saved |
| fragment matching | b and y ions, charge 1 and 2, +/- 20 ppm, at most 150 peaks (Sage `max_peaks` 150, asserted) | as the search | none |
| resources | 1 CPU, 32 GB, 4 h, one task at a time | measured peak 0.7 GB | none |

## Resources and run time (measured)

Revision 2 pilot (job 5930491, 1.6 min, 2 CPUs allocated), CapScan (5930492, 16.1 min, peak 0.72 GB), 32 GB requested; the PREDICT task of this wave (5930493) was cancelled and the PREDICT strata were produced by the recovery job 5937598 (32 GB, 13.7 min, 0.7 GB; `terminal_histidine_inference/collect_repaired_predict_20260915.py`).

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
