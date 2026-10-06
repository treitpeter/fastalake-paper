# campi_simplified_grouping

Internal ID: A40 (previous folder name A40_campi_simplified_grouping)

Question. How many protein groups, and groups with two or more assigned peptides, does the simplified grouping give per CAMPI acquisition in the five arms? Fig2 c, line 77.

Code in this folder: collect.py, hstool_full_comparison.py, prepare.py, run_cohort.py.
Launchers (slurm/): collect_tasks.sh, run_task.sh.
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

manifest.json of the run (cohort CAMPI, five arms, 15 acquisitions each) and rosters/CAMPI.txt (sha256 below); Sage results and annotations of each arm.

## Outputs

cohorts/CAMPI_*/ per-acquisition tables and collection/ tables (hashes below). Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

array: 8 tasks, 2 CPU, 16 GB request, peak RSS 3.58 GB, longest 35.3 min (these tasks cover HSTOOL, PREDICT, CAPSCAN and CAMPI; the CAMPI tasks are the shortest-running subset, not separated here). collection 1.3 min, 4 GB.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| methods compared | original and simplified grouping | comparison named in the run manifest | METHODS in run_cohort.py |
| annotation levels | OG, KEGG_ko, KEGG_Pathway and others listed in LEVELS | as in run_cohort.py | LEVELS constant |
| code commit | 2791c207f11c844fb24e403320d44493de2fa776 of the FastaLake package | recorded in the run manifest | none |

## Notes

run_cohort.py imports fasta_lake.study from the FastaLake package and sits at publication/analyses/<name>/ inside it (ROOT = parents[3]); the package is the repository code, not duplicated here. FASTALAKE_CODE_DIR must point to a checkout at the commit above with these files in place.

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
