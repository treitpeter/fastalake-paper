# Methods statistics

This unit recomputes the numbers that the main Methods (v28d lines 197-307) prints and that
no figure or table unit carries: reservoir sizes, the sizes of the supplied CAMPI
comparator FASTAs, the five reconverted PREDICT mzML files, the cross-fit roster, the
exact-sequence unions of the selected FASTAs, the resampling settings of the
spectrum-split and antibiotic analyses, and the CAMPI compute times. It has no figure
panels, so it has no pixel checks, `ASSEMBLY.md` or `ARTWORK_NUMBERS.md`.

Run from the repository root: `python paper.py build methods`. The antibiotic test
(line 275) is recomputed only when `FASTALAKE_CONTROLLED_DIR` holds
`paired_antibiotic_profiles.tsv`; otherwise the notebook prints the public aggregate table.

`prepare_source_data.py` documents how `source_data/methods/` was cut from the lab tree
(run once on Slurm with `FASTALAKE_LAB_ROOT`); the notebook does not call it.

## Statistics

| statistic | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| reservoir sizes | public, PREDICT, CapScan and CAMPI reservoirs; assembly FASTAs pooled | 215 | 1 reservoir each; sequence | none (count) | Line 215: reservoirs | `source_data/methods/reservoir_sources.tsv` (check: `source_data/Fig1/acquisition_sizes_726.tsv`) | cohort_database_sizes |
| comparator FASTAs | records and contaminants of the CAMPI reference and MG FASTAs | 227 | 29,673 and 19,435 records; sequence | none (count) | Line 227: CAMPI comparator FASTAs | `source_data/methods/campi_comparator_records.tsv` | campi_comparator_fastas |
| reconverted mzML | which PREDICT searches used reconverted files | 233, 303 | 306 acquisitions | none (count) | Lines 233 and 303 | `source_data/methods/predict_reconversion.tsv` | sage_input_completeness |
| spectrum-split bootstrap | 95% intervals of acquisition-median FL/MG gains | 267 | 296 acquisitions from 97 patients; 282 acquisitions from 15 participants | percentile participant-cluster bootstrap | Line 267: spectrum-split bootstrap | `source_data/methods/split_spectra_mg_per_acquisition.tsv`, `split_spectra_mg_summary_original.tsv` (cluster order: `source_data/Fig4/crossfit_per_acquisition.tsv`) | metagenome_split_spectra_comparison |
| cross-fit roster | matched acquisitions, singular pooled fits, interpreter | 267 | 578 acquisitions | none (count) | Line 267: cross-fit roster | `source_data/Fig4/crossfit_per_acquisition.tsv`, `source_data/methods/crossfit_lda_fallback_347.tsv`, `crossfit_interpreter.tsv` | A95_spectrum_crossfit |
| selected-FASTA unions | exact-sequence union sizes | 269 | 58, 306, 331 and 347 FASTAs; protein sequence | none (count) | Line 269: study unions | `source_data/methods/union_contributions.tsv` | cohort_selected_fasta_unions |
| antibiotic comparison | resampling settings and n of the on-minus-off test | 275 | 42 patients | Monte Carlo paired sign-flip test; percentile patient bootstrap | Line 275: antibiotic comparison | `source_data/Fig4/clinical_paired_statistics.tsv`; controlled `paired_antibiotic_profiles.tsv` | predict_clinical_refit |
| compute | reservoir build, extraction, selection, search and inference times; node | 299 | 1 job; 15 tasks per stage | none (medians) | Line 299: compute times and hardware | `source_data/SFig15/compute_task_records.tsv`, `campi_alphanovo_inference_sacct.tsv`, `source_data/methods/compute_nodes.tsv` | common_resource_accounting |

The cross-fit intervals and the antibiotic estimates themselves are printed in the Results
and Figure 4 and belong to `figures/Fig4`; the split-spectrum intervals are printed in
Supplementary Figure 6. This unit carries the Methods statements about how they were made.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `SPLIT_BOOT_DRAWS` | 10,000 | Methods line 267 | none |
| `SPLIT_BOOT_SEED` | 20261003 | seed of the original analysis (not in the text) | none |
| `SPLIT_CI_PERCENTILES` | 2.5, 97.5 | percentile interval of the original analysis | none |
| cluster order | `participant_sort_rank` | reproduces the original draw order without raw patient identifiers | pseudonym order: PREDICT bounds then differ by up to 0.05 points (`split_bootstrap(order_by_rank=False)`) |
| `CLIN_BOOT_DRAWS` | 50,000 | Methods line 275 | none |
| `CLIN_SIGNFLIPS` | 200,000 (100 batches of 2,000) | Methods line 275 | none |
| `CLIN_SEED` | 20260907 | seed of the original analysis (not in the text); sign flips use seed + endpoint index | none |
| `CLIN_TIE_TOL` | 1e-12 | tolerance of the original exceedance count | none |
| `CAMPI_LAKE_STAGE`, `CAMPI_STAGES`, `POOLED_EXTRACTION_STAGE` | lake; 05_persample, 06_parsimony, 07_sage; 02_evidence | stage names of the CAMPI chain run | none |

## Definitions

- Contaminants (line 227): records whose accession occurs in both comparator FASTAs. The MG
  FASTA otherwise holds only assembly proteins (`k99_` accessions); the 116 shared records
  also have identical sequences in both files.
- Union identity (line 269): exact full sequence, upper case, whitespace removed, I and L kept
  distinct, as in the original `build_union` (`prepare_source_data.py:union_contributions`).

## Known limitations

- The selected FASTAs themselves are not in the repository (about 1 GB). `union_contributions.tsv`
  records each FASTA's sha256 and its contribution to the union; the union recount from the
  FASTAs is in `prepare_source_data.py` and needs the deposited FASTAs.
- Node memory and features (`compute_nodes.tsv`) were read with `scontrol` on 2026-10-06, not at
  run time (2026-09-09). "AMD EPYC" rests on the node feature `genoa`; the CPU model string
  was not recorded. "One NVIDIA L40S GPU" per task is not in the accounting records (the nodes
  carry four L40S each).
- The CAMPI comparator FASTAs are the frozen copies used by the comparator arms; the
  supplied files with decoys were not located, so "We removed supplied decoys" is not checked.
- Methods numbers of other analyses (rosters, thresholds, families) are not in this unit;
  see `OPEN_ITEMS.md`.
