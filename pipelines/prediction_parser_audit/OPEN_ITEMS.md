# Open items: prediction_parser_audit

- PIPELINE: input files are read from LAB_ROOT paths and their sha256 were not recorded by the original runs; no input hash is given.
- UNTRACED: code path given only as prose in the plan and not copied: 6-10-01/biology_refit/claude_review_2026-10-03/final/reruns_347/A7_other/parser_audit/bin (Rust binary m07_parser_audit)
- PIPELINE: the original has no Slurm launcher in the listed paths; no slurm/ folder.
- UNTRACED: output file parser_audit_acquisitions_726.tsv was not found under the expected original directory; no sha256.
- PIPELINE: the Rust binary m07_parser_audit is a compiled artefact; its source is not in the plan and the binary is not copied.
