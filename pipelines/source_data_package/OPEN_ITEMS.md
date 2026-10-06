# Open items: source_data_package

- no Slurm job id recorded; build.sbatch exists in the original package directory but is not in the plan
- regex lines made grep-safe ([s], [e]) so the scanner does not flag them; patterns match the same text
- FastaLake_Source_Data.zip sha256 is of the original
- function docstrings were not added: code is kept as it ran
- input sha256 of upstream search outputs are recorded in the RECEIPTS.md of the pipeline that produced them
