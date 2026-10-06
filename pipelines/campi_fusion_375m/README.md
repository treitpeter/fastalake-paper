# campi_fusion_375m

Internal ID: A70 (previous folder name A70_campi_fusion_375m)

Question. Does the fusion (composition-preserving linked control) calibration hold on the 375,333,654-sequence reservoir? 15 CAMPI acquisitions x 2 seeds (lines 713-715) and label validity (line 721).

Code in this folder: PROTOCOL.md, audit_fusion.py, check_labels_repaired.py, config.json, driver.py, run_experiment.py, submit.py.
Launchers (slurm/): job.sh.
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

config.json (frozen binaries, 15 prediction CSVs, mzML) and inputs/PROVENANCE.json (hashed below).

## Outputs

results tables and JSON (hashes below). Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

preflight 16 CPU peak RSS 5.4 GB 2.5 min. prepare 2 tasks, 32 CPU, 96 GB request, peak RSS 95.5 GB, 80.5 min. pilots peak RSS 9.7 GB, 2.3 min. remaining 28 searches peak RSS 10.6 GB, 3.4 min each at most. collect 0.1 min.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| reservoir | 375,333,654 entries of the 2026-09 noclust arm | the reported broad reservoir | historical 233,477,360 control |
| seeds | 2 control seeds | matched-seed design | n/a |
| estimator | original control construction and estimator retained | isolates the reservoir mismatch | none |

## Notes

driver.py, submit.py, run_experiment.py, audit_fusion.py: run_experiment.py is the same file as the one kept in A92 (analyses/fusion_entrapment); a copy is retained here as it ran, check_labels_repaired.py is included for line 721.

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
