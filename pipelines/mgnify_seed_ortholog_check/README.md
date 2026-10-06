# mgnify_seed_ortholog_check

Internal ID: A104 (previous folder name A104_mgnify_seed_ortholog_check)

Question. Do native MGnify eggNOG annotations survive seed-ortholog verification by DIAMOND? 9,368,708 rows aggregated per genome (line 657).

Code in this folder: prep_native.py, summarize_verification.py, verify_shard.py.
Launchers (slurm/): none (no Slurm launcher: the code runs on one core in seconds, or the original launch lines were not recovered).
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

Native annotation rows prepared by prep_native.py; eggNOG DIAMOND database (FASTALAKE_EGGNOG_DMND) and diamond binary (FASTALAKE_DIAMOND).

## Outputs

verify/VERIFICATION_SUMMARY.json, native_row_verification.tsv.gz, INDEPENDENT_MAPPER_CONCORDANCE.json (hashes below); mgnify_seed_ortholog_check_by_genome.tsv is written by summarize_verification.py. Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

prep 4 CPU, 24 GB request, peak RSS 1.56 GB, 7.6 min. verify 64 tasks, 8 CPU, 6 GB request, peak RSS 0.56 GB, longest 5.3 min. control 4 tasks 4.4 min. summary peak RSS 6.09 GB, 3.5 min.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| mapper | DIAMOND against the eggNOG protein database | not stated in the files read | independent mapper concordance file |
| shards | 64 | parallelism only | n/a |

## Notes

summarize_verification.py was copied from the plan; there is no Slurm launcher file for these jobs in the sources (launch lines are in RECEIPTS as not recovered).

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
