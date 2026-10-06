# Assembly: Supplementary Table 5

The table is a workbook, not a composite figure. It was assembled in two scripted steps,
both reproduced by `panels.build_table` in the notebook section `Table`. No manual step.

1. Producer workbook. Input: `crossomics_same_function_correlations.tsv`,
   `same_function_summary.tsv`, `eligibility.tsv` (revised allocation) and the
   previous-allocation tests. Tool: `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/cr003_followup/scripts/tables_cr003.py`
   (Supplementary Table 4 block, variant verified50; the table was numbered 4 before the
   v21 renumbering). Writes sheets `All_30581_tests`, `Nine_family_summary`,
   `All_eligibility` with pandas, adds a `Read_me` note sheet in front and an
   `Eligibility_common` sheet computed from the outer join of revised and previous tests.
   Output: `.../iteration2/cr003_followup/tables/verified50/Supplementary_Table_4.xlsx`,
   sha256 c7dd4195c2404825cb89e2378e00959885cc467e67e87cc48b479c2899bae65e.
2. Source Data package edit. Tool: `LAB_ROOT/.../final/figure_sd/source_data/build/build_package.py`
   with `specials.py` function `st4`. Removes the producer's `Read_me`, inserts the
   package `Read_me` (ten lines, `panels.PACKAGE_README`), renames four
   `Eligibility_common` columns (`n_new` to `n_newly_eligible`, `median_refit_common` to
   `median_rho_revised_common`, `median_baseline_common` to `median_rho_previous_common`,
   `median_refit_new` to `median_rho_revised_newly_eligible`), and sets neutral document
   properties. Output: `FastaLake_Source_Data/Supplementary_Tables/Supplementary_Table_5.xlsx`
   (copied as `expected/Supplementary_Table_5.xlsx`), sha256
   cef0a5eb4c7774120fc05e021660210a25b4220101c7b4fc72595ca4b63c5641.

The same build writes the two support files from the join of step 1 with the package
functions `elig` and `elig_sum` (labels `both`/`left_only`/`right_only` to
`common`/`newly_eligible`/`lost`; refit/baseline columns to revised/previous allocation):
`expected/eligibility_matched_correlations.tsv` (sha256
d1ad7b4296d71f20fc160a568eaec576631021f87ce4afca0b98f6351fa719aa) and
`expected/common_eligibility_summary.tsv` (sha256
cbbfa654762645f458be66a90768d6af2371caab597fff3f4a7d77e42050db00).

The rebuilt workbook equals `expected/Supplementary_Table_5.xlsx` in every cell
(`checks/table.json`); its bytes differ only by openpyxl timestamps. The v28f table is the
v28d table; nothing changes for v28f.
