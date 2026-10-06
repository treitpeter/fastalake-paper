# Assembly: Supplementary Table 2

How the published workbook was put together. All steps were scripts; no manual edit was
found. The rebuild in `STab2.ipynb` reproduces the cell values of the final file in one
step (`stab2_table.write_workbook`); it does not reproduce the bytes (zip timestamps,
styles, frozen panes and filters are not values).

1. Inputs: `LAB_ROOT/canonical/publication/analyses/full_route_entrapment_2026-09-06/tables/<cohort>/{per_acquisition,stage_control_counts,baseline_comparison,summary,control_discoveries}.tsv`.
   Tool: `LAB_ROOT/../BiorXiv/99_archive/2026-09-11_code_exports_before_edit6.1/08_code/analyses/full_route_entrapment_2026-09-06/build_package.py`, `main` (Slurm job 5877268, `publication_build`).
   Change: concatenated the four cohorts, computed Cohort_overview, wrote seven sheets (Read_me .. Control_discoveries).
   Output: `SupplementaryTable2_full_route_entrapment.xlsx`, sha256 `6b437b2ebcbb7d77f28ed578f98e405381d15c00a1a9fb37e04eccc347a4b42e`.
2. Inputs: step 1 output; `SourceData_fusion_completed.xlsx` (sha256 `de61df5f138e9d6f87ed6de15162823eb864ff9cb8751a1c0bc5100ee5671765`, written by `package_streamlining_sources.py`, `main`); `SourceData_CAMPI_independent_selection.xlsx` (sha256 `835be231c956ed032d9946f2ea273e4034f4ff4e6e439c66e2a4f60c58e163fd`, written by `package_v10_sources.py`, `main`).
   Tool: `LAB_ROOT/../BiorXiv/08_code/edit6.1/figures/upstream_supplement_layout/package_edit3.py` (module-level code, table 2 branch).
   Change: appended every sheet of the two workbooks with prefixes `fusion_` and `indep_` (names cut to 31 characters; `fusion_acquisition_outliers_q001` became `fusion_acquisition_outliers_q00`); set frozen panes and filters on the appended sheets.
   Output: `SupplementaryTable2_full_route_entrapment.xlsx`, sha256 `2d33ce310b38b0a04b1201edaf40cf1bcc193f22fb40f5dcd321321ff024c843` (delivered as `LAB_ROOT/00_PACKAGE_2026-09-17/.../supplements/revision_2026-09-24/delivery/source_data/tables/CONTROL_EXPOSURE/SupplementaryTable2_full_route_entrapment.xlsx`, registry key CONTROL_EXPOSURE).
3. Input: step 2 output. Tool: `LAB_ROOT/00_PACKAGE_2026-09-17/.../v15_MM_coauthor_2026-09-25/package_v15.py` (`scrub`, table_map 2 = CONTROL_EXPOSURE).
   Change: renamed to `Supplementary_Table_2.xlsx`; replaced absolute paths inside the XML parts (none occur in cell values of this workbook; the rebuilt values equal the final ones).
   Output: `coauthor_package/supplementary_tables/Supplementary_Table_2.xlsx`, sha256 `5eb86111a0eacf5983c3f750d302281fba83f6c5c44700ce10b82efa09b1027d`.
4. Input: the v16 Source Data archive `LAB_ROOT/review/FastaLake_v16_SOURCE_DATA.zip`, member `current/Supplementary_Tables/Supplementary_Table_2.xlsx` (14,395,034 bytes). Tool: `LAB_ROOT/review/v15_integration/publish_ready_v16.py` and `.../final/figure_sd/source_data/build/selection.py` (rule Z16, `Supplementary_Tables/` include).
   Change: copied unchanged into the Source Data package.
   Output: `FastaLake_Source_Data/Supplementary_Tables/Supplementary_Table_2.xlsx`, sha256 `4768d39fff4b3fa3605a5cc5af504b084356a98db232f52a2cb8f246dcedbc35` = `expected/Supplementary_Table_2.xlsx`.

The byte change between step 3 (`5eb86111...`) and step 4 (`4768d39f...`) was not traced to a
script (see OPEN_ITEMS.md); cell values of step 4 equal the rebuild.

Final: `expected/Supplementary_Table_2.xlsx` sha256 `4768d39fff4b3fa3605a5cc5af504b084356a98db232f52a2cb8f246dcedbc35`.
There is no composite image or PDF for this table.
