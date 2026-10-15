# Cross-fitted entrapment controls

Internal ID: A-16 / D-R03.

The producer measures paired entrapment errors after each search database is
fixed, across all 591 matched acquisitions. It compares cross-fitted FastaLake
(FL), the matched metagenome (MG), the all-prediction FastaLake database (FULL),
and peptides accepted only by FastaLake (GAINED). Both q-only and q-plus-PEP
filters and all three peptide-length ranges are retained.

- `src/ecv.py`: search setup, pooled rescoring, peptide/PSM error estimators and
  participant-bootstrap summaries.
- `src/make_entrapment.cpp`: paired tryptic-block sequence shuffling.
- [Public inputs and summaries](../../source_data/v42_updates/README.md).
- [Source identities](RECEIPTS.md).

The upstream producer requires the original study workspace: raw spectra,
prediction tables, fixed fold FASTAs, Sage configuration and executable, and the
frozen cross-fit producer whose hashes are checked in the code. `LAB_ROOT` marks
the project root to configure for that workspace. Build the C++ program as
`bin/make_entrapment` before running `src/ecv.py run <index>`; `collect` combines
completed acquisitions. Clinical manifests are obtained through the
[controlled access route](../../docs/DATA_AVAILABILITY.md).

For reproduction from the included public tables, use
`python tools/check_updated_analyses.py --bootstrap` from the repository root.
That command repeats summary calculations, not searches.
