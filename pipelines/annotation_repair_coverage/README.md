# annotation_repair_coverage

Internal ID: A105 (previous folder name A105_annotation_repair_coverage)

Question. What share of compatible proteins in each study arm carries a gene-sequence-verified MGnify annotation? annotation_mgnify_coverage_by_arm.tsv (8 arms plus the 347-roster CapScan arm).

Code in this folder: collect.py, repair_projection.py.
Launchers (slurm/): none (no Slurm launcher: the code runs on one core in seconds, or the original launch lines were not recovered).
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

evidence_repaired/INPUTS.json and CODE_MANIFEST.json (hashed below).

## Outputs

full_collection tables (hashes below). Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

8 tasks, 2 CPU, 40 GB request, peak RSS 24.05 GB, longest 26.1 min. collect under 1 min.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| annotation rule | gene-sequence-verified MGnify rows only | avoids name-only transfer | repair_projection.py |

## Notes

The CapScan 347-roster arm was produced in a separate follow-up (capscan_followup_2026-10-01) that is not in the plan; open item. The job ids above are the reproject jobs, matched by name to repair_projection.py, not read from a manifest.

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
