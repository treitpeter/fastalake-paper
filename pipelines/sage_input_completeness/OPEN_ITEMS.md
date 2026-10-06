# Open items: sage_input_completeness

- PIPELINE: input files are read from LAB_ROOT paths and their sha256 were not recorded by the original runs; no input hash is given.
- UNTRACED: no Slurm job id could be traced; run time and memory are not given.
- PIPELINE: the original has no Slurm launcher in the listed paths; no slurm/ folder.
- UNTRACED: output file search_input_completeness.tsv was not found under the expected original directory; no sha256.
- PIPELINE: reconverted-flag and end-of-file checks of the 306 PREDICT mzMLs (lines 233, 303) are not traced to a script in the plan.
