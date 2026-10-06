# campi_reference_gap

Internal ID: A14 (previous folder name A14_campi_reference_gap)

Question. How many CAMPI member reference peptides are recovered by each arm, where are the others lost, and what are the UpSet intersections and length tables? Fig2 b,d,e,f; SFig03; line 35 and lines 683-687.

Code in this folder: aggregate.py, build_revision.py, derive_upset_a14.py, search.py.
Launchers (slurm/): aggregate.sbatch, audit.sbatch, stage_frozen.sh.
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

inputs/manifest.json and inputs/aggregate_manifest.json of the run: accepted 5-50 and 9-50 I/L-canonical peptide sets of the five completed Sage arms (hashes saved), cleaned predictions, candidate and final FASTAs, CAMPI OS-to-member map.

## Outputs

results/member_counts_all.tsv, member_summary.tsv, paired_reference_recovery.tsv, recovery_summary.tsv, AGGREGATE_VALIDATION.json (hashes below); derive_upset_a14.py derives the UpSet tables; build_revision.py is not included (see open items). Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

pilot 2 CPU, peak RSS 0.13 GB, 2.2 min. audit 14 tasks, peak RSS 0.18 GB, longest 3.1 min. aggregate 0.1 min. Requests: 2 CPU, 8 GB.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| digestion | Sage 0.14.6 KR/P, two missed cleavages, no initiator-methionine clipping | matches the production search | none |
| peptide windows | 5-50 and nested 9-50, I/L canonical | as in the main search | both reported |
| beams | rank zero only | exports contain rank zero; unavailable beams are not claimed | none |
| searches | none rerun | recovery uses the five completed arms | n/a |

## Notes

Shares lib.py and peptide.py with A03.

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
