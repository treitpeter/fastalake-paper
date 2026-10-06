# campi_alphanovo_accuracy

Internal ID: A44 (previous folder name A44_campi_alphanovo_accuracy)

Question. How often is the top AlphaNovo prediction identical to the REF-search peptide of the same spectrum? 51.8% of 801,613; median 67.7%, range 26.3-79.2%; 85.3% of 271,469 (line 81); permuted control 0.006% (line 553).

Code in this folder: accuracy.py.
Launchers (slurm/): run.sbatch.
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

The 15 prediction CSVs and archived REF results of A03 (hashed in the A03 manifest).

## Outputs

results/per_acquisition.tsv and pooled.tsv (hashes below); the manuscript file alphanovo_accuracy_per_acquisition.tsv is per_acquisition.tsv. Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

1 CPU, 640 MB request, peak RSS 0.52 GB, 1.0 min. Request = largest input scaled from pilot x 1.05.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| query filter | production lib.load_queries (score > 0, 9-50, complexity, top 30%) | same as lookup | all-rows window also reported |
| control | prediction sequences permuted across spectra within acquisition | exact-match background | one permutation |
| windows | 9-50 and 5-50 | as elsewhere | both |

## Notes

Uses lib.py and peptide.py of A03.

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
