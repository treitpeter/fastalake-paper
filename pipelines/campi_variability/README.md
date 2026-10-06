# campi_variability

Internal ID: A43 (previous folder name A43_campi_variability)

Question. How much of the add-back gain is available in the predictions, how many peptides are newly quantified, and what does equal-MS2 depth standardisation give? 98.4%, 99.9%, 3.38-fold, 1.98-fold (line 77); 38,254 MS2 spectra and 76,127 pooled reference peptides.

Code in this folder: analyze.py, build_mechanistic_review.py, mechanism_v2.py.
Launchers (slurm/): none (no Slurm launcher: the code runs on one core in seconds, or the original launch lines were not recovered).
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

INPUTS.json (hashed below) and inputs/figure_metrics.tsv, workflow_metadata.tsv; results of the five completed arms.

## Outputs

mechanism/gain_availability.tsv, new_quantified_peptides.tsv, gain_spectrum_fates.tsv, upstream_prediction_stages.tsv, group_support/paired_group_support.tsv (hashes below). Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

all jobs 2 CPU, 8 GB request; full run 2.1 min, peak RSS 0.51 GB; others under one minute.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| depth standardisation | equal MS2 count per acquisition | removes spectra-number differences between arms | reported in the same table |
| mechanism version | mechanism_v2.py | v1 label error corrected (MECHANISM_LABEL_CORRECTION.json) | n/a |

## Notes

group_support.py and mechanism.py (v1) are not in the plan list and not included; campi_spectral_inventory.tsv and campi_reference_peptides.tsv were not found by name (open item).

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
