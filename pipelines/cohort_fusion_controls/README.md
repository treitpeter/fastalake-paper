# cohort_fusion_controls

Internal ID: A92 (previous folder name A92_fusion_completed_2026_09_08)

Question. Completed fusion-control calibration on four cohorts (CAMPI mock, HSTOOL, PREDICT, CAPSCAN) with two seeds. Sheets of STab2.

Code in this folder: SPEC.md, audit_completed.py, build_fused.cpp, config.json, run_experiment.py.
Launchers (slurm/): job.sh, prepare_array.sh, sample_array.sh.
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

config.json (frozen binaries and sample lists, cohort lakes) and SOURCE_SHA256.json of the run.

## Outputs

cohort_q001.tsv, cohort_threshold_summary.tsv, matched_seed_summary.tsv, short_5_to_8_summary.tsv, acquisition_outliers_q001.tsv, completeness.tsv, snapshot/per_acquisition.tsv, VALIDATION.json (hashes below). Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

prepare v2: 8 tasks, 32 CPU, 96 GB request, peak RSS 95.6 GB, 61.8 min. searches: 16 CPU, 64 GB request, peak RSS 32.7 GB, at most 8 min per task; 1,237 search tasks in total.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| seeds | 20260907 and 20260908 | two independent control seeds | n/a |
| FDR level | q 0.01 | reporting threshold | threshold summary table |
| prepare v1 | failed (5881251) and was rerun as v2 (5881273) | see RECEIPTS | n/a |

## Notes

run_experiment.py, job.sh, sample_array.sh, prepare_array.sh, build_fused.cpp, SPEC.md, config.json are the experiment code of analyses/fusion_entrapment_2026-09-07. The compiled binaries (fasta_extractor, parsimony_engine, sage, build_fused_v2) are not in the repository (open item).

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
