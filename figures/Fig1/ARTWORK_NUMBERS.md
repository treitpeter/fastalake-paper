# Figure 1, panel a: typed numbers

Every number typed in the schematic, as extracted from the text layer of
`artwork/panel_a.pdf` (words above y = 372 pt). Since v30, `artwork/panel_a.pdf` is the
panel-a crop of `artwork/Figure1_v33.pdf` (the v33 funnel build, Fig1_v3 layout; v30 used `artwork/Figure1_v30.pdf`). Values are recomputed by
`panels.py:artwork_numbers` and printed in the notebook section `Panel a`. The v30 column has
no step badges, no lookup box and no "Pool queries across 58 acquisitions" box; those
rows were removed (the 58 is stated in the legend and in NUMBERS.tsv, line 71).

| artwork file | typed text | value | table it comes from (repo path, row, column) | status |
|---|---|---|---|---|
| artwork/panel_a.pdf | 20,431 reference proteins | 20,431 | source_data/Fig1/human_reference.json, key `entries` | MATCH |
| artwork/panel_a.pdf | Gut sequence reservoir 233 M proteins | 233 M (233,497,702) | source_data/Fig1/stage1_stats.json, key `proteins_searched` | MATCH |
| artwork/panel_a.pdf | table "Proteins in reservoir", Real 233 M | 233 M | source_data/Fig1/stage1_stats.json, key `proteins_searched` | MATCH |
| artwork/panel_a.pdf | table "Proteins in reservoir", Scrambled 233 M | 233 M | source_data/Fig1/stage1_stats.json, key `proteins_searched` (the scrambled arm was looked up in the same reservoir, Methods line 665) | MATCH |
| artwork/panel_a.pdf | table "Study-supported proteins", Real 18 M | 18 M | source_data/Fig1/stage1_stats.json, key `proteins_matched` | MATCH |
| artwork/panel_a.pdf | AlphaNovo peptide predictions, 24,560 unique peptides | 24,560 | source_data/Fig1/selection_grid_acquisitions.tsv, row `arm == tp0.30`, `sample == 20181231_QX7_PhGe_LC95_SA_WorkflowReplicates_2`, column `s3_seeds`; equal to the sum of `distinct_queries` over `arm == REAL` in source_data/Fig1/length_by_arm.tsv | MATCH |
| artwork/panel_a.pdf | Composition-matched scrambled peptides, 36,588 unique peptides | 36,588 | source_data/Fig1/length_by_arm.tsv, sum of `distinct_queries` over rows `arm == SIM_20260901` | MATCH |
| artwork/panel_a.pdf | table "Candidates, this acquisition", Real 6.67 M | 6.67 M (6,671,447) | source_data/Fig1/scramble_end_to_end.tsv, row `arm == authentic`, column `selected_proteins`; equal to `s3_proteins` of the selection-grid row above | MATCH |
| artwork/panel_a.pdf | table "Study-supported proteins", Scrambled 5,841 (new in v30) | 5,841 | source_data/Fig1/scramble_end_to_end.tsv, row `arm == all_position`, `seed == 20260901`, column `selected_proteins` (only this acquisition was scrambled, so its study-supported and candidate sets are the same; legend v30) | MATCH |
| artwork/panel_a.pdf | table "Candidates, this acquisition", Scrambled 5,841 | 5,841 | source_data/Fig1/scramble_end_to_end.tsv, row `arm == all_position`, `seed == 20260901`, column `selected_proteins` | MATCH |
| artwork/panel_a.pdf | table "Sample FASTA after parsimony", Real 12,606 | 12,606 | source_data/Fig1/scramble_end_to_end.tsv, row `arm == authentic`, column `parsimony_proteins` | MATCH |
| artwork/panel_a.pdf | table "Sample FASTA after parsimony", Scrambled 525 | 525 | source_data/Fig1/scramble_end_to_end.tsv, row `arm == all_position`, `seed == 20260901`, column `parsimony_proteins` | MATCH |
| artwork/panel_a.pdf | 29,663 peptides | 29,663 | source_data/Fig1/scramble_end_to_end.tsv, row `arm == authentic`, `gate == q_only`, `q == 0.01`, `min_length == 5`, `max_length == 50`, column `n_peptides` | MATCH |
| artwork/panel_a.pdf | 0 peptides (scrambled) | 0 | source_data/Fig1/scramble_end_to_end.tsv, row `arm == all_position`, `seed == 20260901`, same gate, column `n_peptides` | MATCH |

In v28f the table cell "Study-supported proteins, Scrambled" was typed as "–". In v30 it
reads 5,841, the same number as the candidate cell below it: the scrambled arm was looked
up only for the worked acquisition, so its study-supported and candidate proteins are one
set (v30 legend). The funnel bands carry labels only ("Gut sequence reservoir", "Exact
match · study-supported proteins", "Per-acquisition match", "Parsimony"); all numbers sit
in the table.
