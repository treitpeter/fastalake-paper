# Changes since the first preprint

The [posted bioRxiv paper](https://www.biorxiv.org/content/10.64898/2026.10.08.757395v1)
will be updated to include these improvements.

- **Stricter held-out testing.** Holding all observations of a peptide/precursor
  group out together gives median peptide-count differences versus the matched
  metagenome search of **−11.9% in PREDICT, −10.7% in CapScan stool and +19.5% in
  intestinal fluid** (9–50 residues). This sensitivity analysis supplements the
  original native-scan comparison. [Results](source_data/v42_updates/group_folds_summary.tsv)
- **Additional error controls.** Paired-entrapment estimates are reported separately
  for complete accepted lists and peptides recovered only by FastaLake, across
  peptide-length ranges and confidence filters. [Results](source_data/v42_updates/entrapment_summary.tsv)
- **Corrected functional annotations.** Figure 5 participant-associated variation
  changes from 45.0% to **44.9%**, and sampling-group variation from 14.8% to
  **15.0%**. All **15 participants** still resemble their own stool more closely.
  Related functional summaries and figures have been refreshed.
  [Figure 5 data](source_data/Fig5/)
- **Expanded validation.** Controls now include all 591 matched acquisitions for
  fold balance and numerical sensitivity, 1,000 label permutations across 58
  healthy-stool acquisitions, and all three EDS split seeds.
  [Validation data](source_data/v42_updates/)

We also refreshed the documentation, linked the preprint and corrected notebook
checks affected by stale source values or machine-level rounding.
