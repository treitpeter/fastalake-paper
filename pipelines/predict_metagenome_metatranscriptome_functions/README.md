# predict_metagenome_metatranscriptome_functions

Internal ID: A89 (previous folder name A89_predict_metagenome_metatranscriptome_functions)

PREDICT metagenome (DNA) and metatranscriptome (RNA) function tables: OG (TPM), KO and pathway abundances per patient specimen, used by STab5 and as the metaG and metaT layers of A84 and A87.

Status: code copied as it ran (path edits only), not re-run in this repository. See RECEIPTS.md for original job ids and hashes.

## Inputs

Not traced. The producing code and its inputs could not be located.

## Outputs

The six tables under `analyses/enterotypes_2026-09-07/` (hashes below). Consumers in that folder (`fullstack.py`, `fullstack_sensitivity.py`, `matched_TPM_FKM.py`, `compare_TPM_matched.py`) only read them.

## Slurm resources and run time (original logs, `sacct`)

Not recorded (producer untraced). Hashing of the six files ran on Slurm.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| producer | not traced | none | none |

## Notes

No code is shipped because none could be traced to the producer of these tables. See OPEN_ITEMS.md.
