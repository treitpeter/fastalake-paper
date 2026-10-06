# Open items: capscan_crossomics

- crossomics_capscan.py loads functions from crossomics_src/check_prevalence.py and native_pyramids.py by AST; copies from source_data/current/Crossomics added
- jobs 6063599 FAILED/6063836 OUT_OF_MEMORY earlier attempts; 6063849 is the used run
- ruff F821 on crossomics_capscan.py is a false positive (names defined by exec)
- function docstrings were not added: code is kept as it ran
- input sha256 of upstream search outputs are recorded in the RECEIPTS.md of the pipeline that produced them
