# Supplementary Table 2

Full-route shuffled-control acceptance and control exposure (v28d line 477). The table
reports, for 710 acquisitions in four cohorts, how many artificial shuffled and foreign
control proteins survived each selection stage, how many accepted peptides were assigned
only to controls, and how the ordinary target identifications changed when the controls
were present. It also carries, as `fusion_*` and `indep_*` sheets, the fusion-linked
control experiment of 2026-09-08 and the CAMPI independent-selection comparison.

`STab2.ipynb` rebuilds the published 18-sheet workbook from `source_data/STab2/` and
compares it cell by cell with `expected/Supplementary_Table_2.xlsx` (result: identical,
`checks/table_compare.json`). It also recomputes the Methods counts that cite this table
(v28d line 743).

## Sheets

| sheet | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| Read_me | definitions | 477 | NA | none | Table | constant text in `stab2_table.py` | NA |
| Cohort_overview | per-cohort medians at q ≤ 0.01, lengths 5-50 | 477 | 15; 58; 306; 331 acquisitions | none (medians) | Table | per_acquisition.tsv, stage_control_counts.tsv, baseline_comparison.tsv | cohort_full_route_controls |
| Per_acquisition | control-only counts and fractions | 477 | 7,100 rows (710 acquisitions x 5 q x 2 length floors) | none | Table | per_acquisition.tsv | cohort_full_route_controls |
| Stage_controls | control proteins per stage | 477 | 710 acquisitions x 4 stages x 3 classes | none | Table | stage_control_counts.tsv | cohort_full_route_controls |
| Baseline_comparison | target peptides with and without controls | 477 | 710 acquisitions | none | Table | baseline_comparison.tsv | cohort_full_route_controls |
| Quantile_summary | acquisition quartiles of the fractions | 477 | per cohort, acquisition | none (quartiles) | Quantile summary | quantile_summary.tsv (recomputed from per_acquisition.tsv) | cohort_full_route_controls |
| Control_discoveries | every accepted peptide with a control assignment | 477 | 366,745 peptide rows | none | Table | control_discoveries_<cohort>.tsv | cohort_full_route_controls |
| fusion_* (8 sheets) | fusion-linked control diagnostics, two seeds | 477 | 710 acquisitions x 2 seeds | none | Table | fusion_*.tsv, fusion_validation.json | cohort_fusion_controls |
| indep_* (3 sheets) | CAMPI fixed, self and independent selection | 477 | 2 acquisitions | none | Table | indep_*.tsv, indep_validation.json | campi_independent_selection |
| not in workbook | shuffle-only peptides, healthy-stool pooling | 743 | 12 acquisitions per arm; peptide-acquisition observations | none (count ratio) | Methods line 743 | hstool_pooling_per_acquisition.tsv | healthy_stool_pooling_controls |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| PRIMARY_Q | 0.01 | primary peptide q threshold (Read_me "Primary") | 0.001, 0.005, 0.02, 0.05 in Per_acquisition and Quantile_summary |
| PRIMARY_MIN_RESIDUES | 5 | primary length window 5-50 (Read_me "Primary") | 9-50 in the same sheets |
| N_ACQUISITIONS | 710 | full-route roster as run (15/58/306/331) | v28f roster 726 not run (OPEN_ITEMS.md) |
| CONTROL_SEED | 20260906 | seed of the shuffled-control construction | none |
| QUANTILE_METHOD | exclusive | `statistics.quantiles` method used by collect_cohort.py | fusion sheets use pandas linear percentiles |
| FUSION_SEEDS | 20260907, 20260908 | two fusion control realisations | both reported |
| POOL_ARMS | own, pool_other11 | the two arms quoted at line 743 | pool_all12 in hstool_pooling_per_acquisition.tsv |

## Files

`stab2_table.py` holds the rebuild and statistics functions; `prepare_source_data.py` wrote
`source_data/STab2/` once from the original outputs and is not run by the notebook.
There is no schematic, so there is no ARTWORK_NUMBERS.md. `ASSEMBLY.md` records how the
published workbook was put together.

## Known limitations

- The full-route experiment covers 710 acquisitions (CapScan 331). The v28f roster
  (CapScan 347, 726 FASTAs) was not run through it.
- The counts at line 743 are not in the workbook; they come from the healthy-stool pooling
  experiment (`hstool_pooling_per_acquisition.tsv`).
- Fractions are observed control-labelled fractions, not calibrated false-discovery
  estimates; quartiles describe acquisitions and are not confidence intervals.
- The `fusion_*` sheets are the 2026-09-08 fusion run; Supplementary Figure 5 in v28d uses
  a later reconciliation of the fusion diagnostics (unit SFig05).
