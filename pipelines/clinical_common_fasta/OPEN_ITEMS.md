# Open items: clinical_common_fasta

- PIPELINE: input files are read from LAB_ROOT paths and their sha256 were not recorded by the original runs; no input hash is given.
- UNTRACED: no Slurm job id could be traced; run time and memory are not given.
- UNTRACED: output files are not located; no output sha256.
- CONTROLLED: per_participant.tsv of PREDICT carries clinical columns; it is not hashed or copied here, and the pseudonymised common_search_participants.tsv is produced outside this folder.
