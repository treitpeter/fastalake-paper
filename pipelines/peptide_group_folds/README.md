# Peptide/precursor-group held-out testing

Internal ID: D-R04.

`src/dr04.py` keeps repeated observations together when defining held-out folds.
Spectra are connected by an identical de novo peptide prediction (I/L collapsed)
or a precursor match with the same charge, within 10 ppm and one minute.
Components are assigned to ten balanced folds before search outcomes are known.
All 591 matched acquisitions are retained, including negative FastaLake gains.

The script builds group maps and cohort caches, searches the held-out spectra,
and collects comparisons with native-scan folds and paired-entrapment controls.
It reuses [crossfit_entrapment](../crossfit_entrapment/) and the frozen upstream
cross-fit design. [RECEIPTS.md](RECEIPTS.md) records source identities;
`ECV_SHA256.txt` pins the adjacent public entrapment producer.

Raw execution requires the original study workspace described in the entrapment
README, with `LAB_ROOT` configured consistently. Modes are `map <index>`,
`cache <index>`, `run <index>`, `negctrl <index>` and `collect`. Prediction and
clinical manifests are not distributed here.

The [public tables](../../source_data/v42_updates/README.md) include acquisition
counts, gains, fold diagnostics, all strata and bootstrap intervals. Run
`python tools/check_updated_analyses.py --bootstrap` to check those summaries
without rerunning the searches.
