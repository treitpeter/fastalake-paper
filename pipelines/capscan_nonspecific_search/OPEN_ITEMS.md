# Open items

- MANIFEST: the builder of the original 331-task and 306-task manifests (`inputs/manifest.json`) was not found in the analysis directories; only the participant-15 builder (`prep_A05_p15.py`) is present. The manifests are hashed in RECEIPTS.md.
- DOCSTRINGS: frozen code is released as it ran, plus a module docstring where the original had none and the `PYTHON` constant read from `FASTALAKE_PYTHON`. Function-level docstrings were not added to the compact original functions.
- SUMMARY TABLE: `nonspecific_search_summary.tsv` named in the Fig6 upstream entry was not found next to the A05 results; it is presumably written by a Fig6 step outside this pipeline.
