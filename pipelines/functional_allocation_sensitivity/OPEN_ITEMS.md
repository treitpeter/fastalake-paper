# Open items: functional_allocation_sensitivity

- PIPELINE: input files are read from LAB_ROOT paths and their sha256 were not recorded by the original runs; no input hash is given.
- UNTRACED: no Slurm job id could be traced; run time and memory are not given.
- UNTRACED: output files are not located; no output sha256.
- PIPELINE: the original has no Slurm launcher in the listed paths; no slurm/ folder.
- PIPELINE: healthy_stool_ko_precision_refit: the measure() and pairs() functions called here live there (profiles.py, hstool_precision.py); the A84 transform/evaluate functions it also calls are not in the plan and are not copied.
