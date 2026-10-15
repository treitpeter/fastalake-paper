# Artwork numbers: Figure 5

Panel a (`artwork/panel_a.pdf`) is the producer's PDF; for v33 the producer was rerun with the matched-acquisition constant 295 (`final/figs_v33_labels`, Figure 5, only the panel-a digits changed). The notebook cell `Panel a`
extracts its text and checks every count below against `source_data/Fig5/design_347.tsv`
(`panels.stat_a_design`); the matched-acquisition count is checked against the CapScan rows of `source_data/Fig4/crossfit_per_acquisition.tsv`.

| artwork file | typed text | value | table it comes from (repo path, row, column) | status |
|---|---|---|---|---|
| artwork/panel_a.pdf | Physiological intestinal sampling · 15 participants | 15 | source_data/Fig5/design_347.tsv, all rows, distinct `participant` | MATCH |
| artwork/panel_a.pdf | 347 metaproteomes | 347 | source_data/Fig5/design_347.tsv, row count | MATCH |
| artwork/panel_a.pdf | 290 capsule · 57 stool | 290 | source_data/Fig5/design_347.tsv, rows with `region` in duodenum, jejunum, ileum, colon | MATCH |
| artwork/panel_a.pdf | 290 capsule · 57 stool | 57 | source_data/Fig5/design_347.tsv, rows with `region` = stool | MATCH |
| artwork/panel_a.pdf | Region 1 · 74 acquisitions | 74 | source_data/Fig5/design_347.tsv, rows with `region` = duodenum | MATCH |
| artwork/panel_a.pdf | Region 2 · 77 acquisitions | 77 | source_data/Fig5/design_347.tsv, rows with `region` = jejunum | MATCH |
| artwork/panel_a.pdf | Region 3 · 74 acquisitions | 74 | source_data/Fig5/design_347.tsv, rows with `region` = ileum | MATCH |
| artwork/panel_a.pdf | Region 4 · 65 acquisitions | 65 | source_data/Fig5/design_347.tsv, rows with `region` = colon | MATCH |
| artwork/panel_a.pdf | Stool · 57 acquisitions | 57 | source_data/Fig5/design_347.tsv, rows with `region` = stool | MATCH |
| artwork/panel_a.pdf | 295 matched acquisitions | 295 | source_data/Fig4/crossfit_per_acquisition.tsv, rows with `cohort` = CAPSCAN | MATCH |

The producer script typed these counts as constants (`COUNTS`, `N_PARTICIPANTS`,
`N_MATCHED`); it did not read them from a table. The internal region ids map to the printed
names as duodenum = Region 1, jejunum = Region 2, ileum = Region 3, colon = Region 4 (device-defined
sampling regions, not anatomical assignments).
