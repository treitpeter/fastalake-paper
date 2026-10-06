# Typed numbers in Supplementary Figure 6

Panel a is a schematic drawn by code (`panels.panel_a`) on the same canvas as b-e and is
pixel-checked with them, so there is no separate artwork file. Its boxes and arrows carry
no numbers.

Numbers typed as labels are listed here with the table that supports them; medians
printed beside the data are computed at drawing time and are not listed.

| artwork file | typed text | value | table it comes from (repo path, row, column) | status |
|---|---|---|---|---|
| SFig06 (b, c x labels) | CAMPI 15 runs | 15 | source_data/SFig06/split_half_acquisitions_726.tsv, rows cohort = CAMPI_MOCK and min_length = 9, count | MATCH |
| SFig06 (b, c x labels) | Healthy stool 58 runs | 58 | same table, cohort = HSTOOL, min_length = 9, count | MATCH |
| SFig06 (b, c x labels) | PREDICT 306 runs | 306 | same table, cohort = PREDICT, min_length = 9, count | MATCH |
| SFig06 (b, c x labels) | CapScan 347 runs | 347 | same table, cohort = CAPSCAN, min_length = 9, count (panels.py `N_CAPSCAN`) | MATCH |
| SFig06 (b annotation) | Two zero-yield runs retained | 2 | same table, cohort = CAPSCAN, min_length = 9, rows with disjoint = 0 | MATCH |
| SFig06 (d x label) | PREDICT (n = 296) | 296 | source_data/SFig06/same_spectrum_benchmark_0_predict.tsv, rows min_length = 9, test_arm = FL_public, reference_arm = MG, count | MATCH |
| SFig06 (d x label) | CapScan (n = 282) | 282 | source_data/SFig06/same_spectrum_benchmark_1_capscan.tsv, same filter, count | MATCH |
| SFig06 (d, e y labels) | 9–50 residues | 9, 50 | both benchmark tables, column min_length = 9 (panels.py `MIN_LENGTH`, `MAX_LENGTH`) | MATCH |

The counts are printed by the notebook (sections Panel b and Panel d). The totals are
asserted in `panels.figure` (726) and `panels.load_benchmark` (296 and 282).
