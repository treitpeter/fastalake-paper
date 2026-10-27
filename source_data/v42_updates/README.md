# Updated analyses

These tables supply the checked-annotation sensitivities, full fold-fraction
census, 20-perturbation rescoring analysis, 1,000-permutation confidence control
and all-seed nested EDS validation adopted in v42. Source identities and staging
transformations are in `PROVENANCE.json`.

`tools/v42_values.py` recalculates the registered ranges, counts and medians from
these tables. Allocation and EDS tables are outputs of the upstream fits; this
check does not repeat those fits. The numerical-sensitivity table retains all
591 acquisitions and all 20 perturbations. The label-control table retains all
58,000 acquisition–permutation outcomes, including zero results.

PREDICT acquisition pseudonyms use the same 306-acquisition key as Figure 4.
Patient identifiers and clinical-group fields were excluded from newly staged
per-acquisition inputs. Original source tables remain in the controlled analysis
workspace.

## Held-out comparisons and error controls

`group_folds_per_acquisition.tsv` and `group_folds_summary.tsv` contain the
peptide/precursor-group held-out sensitivity analysis. Gains are
`100 × (FastaLake / matched-metagenome − 1)`; 5 and 9 denote the lower peptide
length limit, with an upper limit of 50 residues. Both group-fold and native-scan
results, including negative gains, are retained.

`entrapment_per_acquisition.tsv` and `entrapment_summary.tsv` report the paired
entrapment controls for FL, MG, FULL and GAINED (identifications exclusive to FL).
`q` uses peptide q ≤ 0.01; `q_pep` also requires PEP ≤ 0.01. Length ranges are
5–50, 5–8 and 9–50 residues. The paired estimate is `(E + E0 + 2B) / N`, with the
unresolved contribution reported separately; all numerator and denominator
counts are included. Acquisition-level estimates are fractions; summary
estimates and intervals are percentages. `CAPSCAN_all` combines stool and
intestinal fluid.

These tables retain 591 acquisitions (296 PREDICT and 295 CapScan). PREDICT
acquisitions use the existing public pseudonyms, clinical groups are omitted,
and `cluster_order` replaces participant identifiers while preserving the order
needed for participant-cluster bootstrap resampling. The source summaries use
10,000 draws with seed 20261004. Run
`python tools/check_updated_analyses.py --bootstrap` to repeat their medians,
counts, estimator arithmetic and intervals from these inputs. Searches and
peptide assignments are not repeated by this command.

`crossfit_per_acquisition.tsv`, `crossfit_summary.tsv` and
`crossfit_summary_by_region.tsv` contain the current native-scan comparison:
all 295 CapScan acquisitions use the 347-acquisition cache roster, alongside
296 PREDICT acquisitions. Earlier cross-fit tables in the Supplementary Methods
unit are retained for the historical comparisons in that notebook. Use these
updated tables for the current comparison.

The Supplementary Methods 6–7 notebook uses `tools/current_crossfit.py` to
recompute 30 registered current claims from this acquisition table, including
the participant-bootstrap bounds. It retains the historical calculations in
the preceding cells. The archived fallback-control table remains an additional
input for the two PREDICT medians calculated after replacing fallback results.
