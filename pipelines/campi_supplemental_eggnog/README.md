# campi_supplemental_eggnog

Internal ID: A106 (previous folder name A106_campi_supplemental_eggnog)

Question. Annotate the 17,935 previously unannotated CAMPI reference proteins with eggNOG-mapper (line 657).

Code in this folder: prepare_supplement.py, run_supplement.py.
Launchers (slurm/): none (no Slurm launcher: the code runs on one core in seconds, or the original launch lines were not recovered).
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

prepare_supplement.py builds the query FASTA from the CAMPI reference proteome and the crosswalk audit.

## Outputs

campi_supplemental_annotation_queries.tsv (not found by name; the closest hashed files are the crosswalk audit, completion check and collect log below). Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

5 tasks, 8 CPU, 18 GB request, peak RSS 3.51 GB, longest 10.8 min. collect 1 s.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| annotator | eggNOG-mapper | same annotation system as the rest of the paper | n/a |

## Notes

The eggNOG-mapper installation and version are not recorded in the files read (open item).

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
