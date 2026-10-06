# campi_candidate_clustering

Internal ID: A81 (previous folder name A81_campi_candidate_clustering)

Question. Does clustering the candidate database at 95/90/75/50% identity change identifications and KO coverage on the 15 CAMPI acquisitions? Table STab1 (105 acquisition-arm configurations) and line 695.

Code in this folder: PROTOCOL.md, collect.py, h04_cluster_retention_v2.py, h08_cluster_searches_v2.py, prepare.py, search.py, selection.py.
Launchers (slurm/): collect.sh, full.sh, pilot.sh.
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

manifest.json of the run (hashed below), PROTOCOL.md, fixed 97/99 results from H04/H08.

## Outputs

figures/per_acquisition_identifications.tsv, per_acquisition_KO_coverage.tsv, paired_medians.tsv (hashes below). Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

pilot 4 CPU, 32 GB request, peak RSS 6.56 GB, 142.5 min. full 13 tasks, peak RSS 5.72 GB, longest 47.1 min. collect 1.2 min, 0.38 GB.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| identity settings | 95, 90, 75, 50% | new settings beyond existing 97 and 99 | 97/99 from H04/H08 |
| FDR | whole-workflow FDR not claimed | functional evidence uses a fixed sequence/KO dictionary | n/a |

## Notes

h04_cluster_retention_v2.py and h08_cluster_searches_v2.py are copies from the earlier paper repository; checked 6 October 2026: no lab paths, e-mails, personal-genome references, review comments, AlphaNovo code or clinical data.

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
