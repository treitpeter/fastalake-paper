# campi_database_metrics

Internal ID: A41 (previous folder name A41_campi_database_metrics)

Question. How many database entries does each of the five arms hold per CAMPI acquisition? Fig2 a and line 77.

Code in this folder: benchmark.py, derive.py, derive_fig2_revision.py.
Launchers (slurm/): none (no Slurm launcher: the code runs on one core in seconds, or the original launch lines were not recovered).
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

Arm directories of the five named database arms under FASTALAKE_DATA_DIR/canonical/publication_data/arms (frozen FASTAs and Sage configurations).

## Outputs

db_records per acquisition and arm, written by derive.py to --out. Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

NOT RECOVERED: no Slurm job id or resource record for derive.py was found.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| arms | blind UniProt noclust, +CAMPI REF, +CAMPI REF+MG, and the two comparators named in derive.py | the five arms of Fig2 | derive.py names dict |

## Notes

derive.py imports derive_fig2_revision.py (included). benchmark.py is the identification benchmark module of the same folder (called through derive_fig2_revision).

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
