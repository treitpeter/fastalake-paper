# campi_independent_selection

Internal ID: A93 (previous folder name A93_campi_independent_selection)

Question. Common-reference CAMPI comparison with independent selection: ordinary counts, fusion diagnostics. indep_* sheets of STab2.

Code in this folder: collect_campi.py, run_campi_reference.py.
Launchers (slurm/): none (no Slurm launcher: the code runs on one core in seconds, or the original launch lines were not recovered).
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

CAMPI_SOURCE_SHA256.json, reference and linked reservoir FASTAs, the 15 acquisitions.

## Outputs

ordinary_counts.tsv, fusion_diagnostics.tsv, VALIDATION.json (hashes below). Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

8 CPU, 24 GB request (12 GB for the second), peak RSS 4.25 GB and 8.98 GB, 2.6 and 4.7 min. Further seed jobs not enumerated.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| selection | independent of the production selection | not stated in the files read | none |

## Notes

Imports run_experiment from A92 (one copy kept there) and peptide from A03.

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
