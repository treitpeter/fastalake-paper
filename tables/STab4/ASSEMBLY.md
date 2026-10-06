# Assembly: Supplementary Table 4

The table is a workbook, not a composite figure. Steps from analysis output to the
published file (paths as in PROVENANCE.md):

1. Input: IT1/eds_refit/input_data/{heatmap_values.npz, heatmap_columns.tsv,
   heatmap_features.tsv, heatmap_85_values.npz, heatmap_85_patient_order.tsv} and
   IT2/shadow/verified50/eds/EDS_approved/source_data/{focused_raw_measurements.tsv,
   focused_features.tsv, selected_pair_summary.tsv, separate_adjustments_all_models.tsv}.
   Tool: IT2/scripts/tables_cr003.py, section "Supplementary Table 3", variant
   verified50 (Slurm job 6048665). What changed: sheets written with pandas
   `to_excel` (20-feature sheets) and openpyxl `append` (Focused_* sheets, freeze panes
   A2), plus a Current_EDS_update note sheet and a Read_me sheet. Output:
   IT2/tables/verified50/Supplementary_Table_3.xlsx, sha256
   46ea54c41fb867a6810e672a391bffd0804aeb190f62ce4e05c48d7b554ccd42.
2. Input: that workbook. Tool: CR/final/figure_sd/source_data/build/build_package.py
   with specials.py function `st3` (selection.py maps it to
   Supplementary_Tables/Supplementary_Table_4.xlsx). What changed: Current_EDS_update and
   Read_me removed, one new Read_me sheet inserted first (text in specials.py
   ST3_README); 426 participant codes replaced by pseudonyms (PREDICT_Pnnn); 1 path
   edit (Read_me last line, "Biology_refit/eds_refit" to "eds_refit"); document
   properties neutralised, re-saved with openpyxl. Output:
   PKG/Supplementary_Tables/Supplementary_Table_4.xlsx, sha256
   4fe89a55c7607b0e5932377cbd488c031bc1f0c7f99c255f6c9eb7d77747e403. This is the v28d
   table; v28f uses it unchanged.
3. Repository: `prepare_source_data.py` (run once) writes source_data/STab4 and
   `expected/Supplementary_Table_4_redacted.xlsx` (sha256
   01314a8f757ae7a996caddfe6c4a6fb20339b94296ead88c583c6a4897b975a0); the notebook
   writes `build/STab4/Supplementary_Table_4.xlsx` (and, with the controlled input,
   `build/STab4/Supplementary_Table_4_full.xlsx`). xlsx files carry save timestamps,
   so equality is checked on cell values, not on file hashes.

No manual step.
