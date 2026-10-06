# Histidine seven-point participant inference (A28)

Internal ID: A28 (previous folder name A28_histidine_inference)

Pipeline id `terminal_histidine_inference`. Used by: Fig6 (panels d to f); text_statistics supplementary_methods_8_11.

## Question

Is histidine over-represented at observed cleavage boundaries relative to the source-protein expectation, per participant region, and are the motifs and complementary products consistent across methods? Answer: participant-level sign-flip inference on the A27 site evidence, with matched controls and concentration sensitivities.

## Inputs

- A27 acquisition receipts (CapScan 331, PREDICT 306, participant 15 added for the 347-acquisition run), `inputs/manifest_revision4.json` (sha256 in `RECEIPTS.md`).

## Outputs

- `results_revision4/{full_CAPSCAN,available_CAPSCAN,available_PREDICT,pilot,statistics_smoke_CAPSCAN,independent_audit}/`: `acquisition_enrichment.tsv.gz`, `acquisition_motif_counts.tsv.gz`, `participant_region_enrichment.tsv`, `participant_motif_differences.tsv.gz`, `boundary_recurrence.tsv`, `matched_reciprocal_controls.tsv`, `participant_reciprocal_controls.tsv`, `participant_concentration_sensitivity.tsv`, `complete_statistical_families.tsv`, `statistical_family_summary.tsv`, `mapping_coverage.tsv`, `P1_classification_coverage.tsv`, `REPORT.md`.
- PREDICT full collection by the recovery job (repaired PREDICT acquisitions), code `collect_repaired_predict_20260915.py`.

## Call chain and code

`collect_revision4.py` imports `lib`, `histidine_helpers` and `context_helpers` from the run directory `frozen/`. `review_histidine_predict.py` independently recomputes the PREDICT rows; `collect_repaired_predict_20260915.py` runs the A20, A24, A25 and A28 collectors on the repaired PREDICT acquisitions (it is the single copy used by those four pipelines).

| file | origin (LAB_ROOT relative) | sha256 at origin | sha256 here |
|---|---|---|---|
| `collect_revision4.py` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/frozen/collect_revision4.py` | `95e49d899ba7e22144287f467ad7a70d1f42e02751fa49ed706c4eb841bdeb74` | `95e49d899ba7e22144287f467ad7a70d1f42e02751fa49ed706c4eb841bdeb74` |
| `audit_inference_revision4.py` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/frozen/audit_inference_revision4.py` | `4c064fd5031fbb5ccfb2b541490ee679c97315a3ac8cc2ebdb8e08a148321cea` | `4c064fd5031fbb5ccfb2b541490ee679c97315a3ac8cc2ebdb8e08a148321cea` |
| `histidine_helpers.py` | `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/frozen/histidine_helpers.py` | `387460e01d28ff90ed17037f1e5ac7d08acb33287b04f503cacb7d1971b888f7` | `387460e01d28ff90ed17037f1e5ac7d08acb33287b04f503cacb7d1971b888f7` |
| `review_histidine_predict.py` | `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/code/review_histidine_predict.py` | `22f10374ac480ca1f9fb7f3ff2f76d4c0d2ec9a919ea935ef2d908a05e2e0b9c` | `7bad9c14166ed92d54797cc9a26ebc2055b271e605b0c1692cea9e253809d680` |
| `collect_repaired_predict_20260915.py` | `LAB_ROOT/MM_review_analysis/code/collect_repaired_predict_20260915.py` | `96289d6e9629c3fd1b00bb5097341ca51601c857758ca383f244808776144e97` | `96289d6e9629c3fd1b00bb5097341ca51601c857758ca383f244808776144e97` |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| test | exact sign-flip test for at least 8 donors and at most 16, otherwise 200000 draws with plus-one correction; observed minus source-expected fractions | small paired samples | assumes sign symmetry; exploratory |
| bootstrap | 20000 draws; separate columns for mean-difference and ratio intervals; seeds 20261014+n (difference) and 20261114+n (ratio), 20260914+n (sign flips) | reproducible | none |
| families | all P1 residues in every available route and side, three primary origins and each CapScan region (PREDICT: all observed specimens) form one source-enrichment family; planned motif contrasts form a second family; missing cells p = 1 | one planned family each | Benjamini-Hochberg |
| motifs | eight positions (P4 to P4-prime), complete denominators and protein-edge missingness, H-conditioned background so P1 = H is not rediscovered as a motif | avoids a circular motif | no pseudocount for zero cells |
| spatial contrasts | ten region pairs plus proximal/stool, distal/stool, proximal/distal; ordinal slopes require all four fluid regions | predeclared | none |
| boundary agreement | same-protein boundary and canonical-peptide agreement are separate outputs; complementary support needs different spectra | distinguishes site from peptide agreement | matching of H against non-H sites within source protein, side opportunity and spectrum-depth bin is descriptive, not a cleavage probability |
| concentration sensitivity | exclude zero, one or five proteins ranked by all non-K/R site counts, with excluded ids saved | checks single-protein dominance | none |
| PREDICT weighting | specimens averaged within patient; no clinical group weighting | no clinical covariates are used | none |
| independent audit | `audit_inference_revision4.py` recomputes donor aggregation and the exact tests from the saved tables and must agree | guards the statistics | job 5931569 |
| resources | 1 CPU, 32 GB, 8 h, one task at a time | measured peak 2.9 GB | participant-15 runs needed 6 GB (3 GB ran out of memory at 2.6 GB sampled peak) |

## Resources and run time (measured)

Revision 4: pilot 5931567 (1.3 min, 0.5 GB), statistics smoke 5931568 (12.1 min, 0.8 GB), independent audit 5931569 (0.2 min, 8 GB requested), available CapScan 5931570 (223.7 min, 2.2 GB), available PREDICT 5931571 (154.4 min, 1.6 GB), full CapScan 5931572 (282.3 min, 2.9 GB); full PREDICT 5931573 was cancelled and replaced by recovery job 5937597 (64 GB requested, 199.6 min, 1.3 GB). All 1 CPU, 32 GB unless stated. Participant-15 runs: 6058108 and 6058148 ran out of memory at 3 GB; completed as 6064030 (baseline 331, 6 GB, 149.9 min, 2.9 GB) and 6064031 (347, 6 GB, 271.3 min, 3.0 GB).

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
