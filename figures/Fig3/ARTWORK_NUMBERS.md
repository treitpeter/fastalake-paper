# Numbers typed in the Figure 3 artwork

Panel a (`artwork/panel_a_design.pdf`, as placed in the v28f composite). Counts per series are
recomputed by `panels.py:stat_a_design_counts` from `source_data/Fig3/transfer_pool_comparisons.tsv`
(rows `kind == transfer`, `arm == own`, `min_length == 5`; distinct `sample` per `group`) and agree
with `source_data/Fig3/design.tsv` (rows per `group`). The number of drawn output symbols per series
equals the typed count.

| artwork file | typed text | value | table it comes from (repo path, row, column) | status |
|---|---|---|---|---|
| artwork/panel_a_design.pdf | 12 injections | 12 | source_data/Fig3/transfer_pool_comparisons.tsv, group TechnicalReplicates, distinct sample | MATCH |
| artwork/panel_a_design.pdf | 15 preparations | 15 | source_data/Fig3/transfer_pool_comparisons.tsv, group WorkflowReplicates, distinct sample | MATCH |
| artwork/panel_a_design.pdf | 10 spatial aliquots | 10 | source_data/Fig3/transfer_pool_comparisons.tsv, group OnePiece, distinct sample | MATCH |
| artwork/panel_a_design.pdf | 11 samples | 11 | source_data/Fig3/transfer_pool_comparisons.tsv, group Intra, distinct sample | MATCH |
| artwork/panel_a_design.pdf | 10 samples | 10 | source_data/Fig3/transfer_pool_comparisons.tsv, group Interindividual, distinct sample | MATCH |
| artwork/panel_a_design.pdf | ×10 (Inter-individual source) | 10 | source_data/Fig3/design.tsv, group Interindividual, distinct subject_id (Ind1 to Ind10) | MATCH |
| artwork/panel_a_design.pdf | Ten donors · one sample each | 10 | source_data/Fig3/design.tsv, group Interindividual, distinct subject_id, one row each | MATCH |
| artwork/panel_a_design.pdf | ×1 (Technical replicates source: one digest) | 1 | no table records the digest; subject_id is empty for this series | OPEN |
| artwork/panel_a_design.pdf | ×1 (Workflow replicates source: one stool aliquot) | 1 | no table records the aliquot | OPEN |
| artwork/panel_a_design.pdf | ×1 (Spatial sampling source: one stool specimen) | 1 | no table records the specimen | OPEN |
| artwork/panel_a_design.pdf | ×1 (Intra-individual source: one donor) | 1 | no table records the donor; subject_id is empty for this series | OPEN |
| artwork/panel_a_design.pdf | 58 acquisitions total (footer; removed in the v28f composite, ASSEMBLY step 12) | 58 | source_data/Fig3/transfer_pool_comparisons.tsv, distinct sample over all five groups | MATCH |
