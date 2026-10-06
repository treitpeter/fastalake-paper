# Open items: bag_mass_paired_nulls

- PIPELINE: input files are read from LAB_ROOT paths and their sha256 were not recorded by the original runs; no input hash is given.
- UNTRACED: no Slurm job id could be traced; run time and memory are not given.
- UNTRACED: output files are not located; no output sha256.
- PIPELINE: the original has no Slurm launcher in the listed paths; no slurm/ folder.
- UNTRACED: manifest.json (2 MB, the task manifest read by submit.py) was not copied; the 60-receipt NULL_COMPARISON.json snapshot was not hashed.
