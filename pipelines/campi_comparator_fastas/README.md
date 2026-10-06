# campi_comparator_fastas

Internal ID: A99 (previous folder name A99_campi_comparator_fastas)

Question. Which frozen CAMPI reference and MG target FASTAs were supplied to the comparator arms? Methods line 227.

Code in this folder: SPEC_ref_append.md, SPEC_ref_metag_append.md.
Launchers (slurm/): none (no Slurm launcher: the code runs on one core in seconds, or the original launch lines were not recovered).
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

SPEC files (copied) describing the two arms.

## Outputs

campi_ref_targets.faa and campi_mo_targets_fixed.faa (hashes below). Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

No computation: these are supplied files; hashing took seconds.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| content | as supplied, no processing | frozen inputs | n/a |

## Notes

No code; the SPEC files are the specification of how each arm used the FASTAs.

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
