# campi_clustering_recount

Internal ID: A82 (previous folder name A82_campi_clustering_independent_recount)

Question. Do independent recounts reproduce the identification and KO tables of the clustering benchmark? Used only by the STab1 consistency check.

Code in this folder: recount.py.
Launchers (slurm/): none (no Slurm launcher: the code runs on one core in seconds, or the original launch lines were not recovered).
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

Outputs of A81 and H08 identification_comparison.tsv.

## Outputs

summary.tsv, per_acquisition_recount.tsv, KO_recount.tsv, config_and_entries.tsv, recount.json (hashes below). Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

NOT RECOVERED: no Slurm job id for recount.py was found.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| method | recount from raw Sage results independent of collect.py | independent check | none |

## Notes

Output names differ from the plan (independent_recount_*.tsv); the files hashed here are the ones in the check folder.

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
