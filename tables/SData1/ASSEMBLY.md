# Assembly: Supplementary Data 1

Supplementary Data 1 is a workbook, not a composite image. It is assembled entirely by
code; there is no manual step.

1. Inputs: `source_data/SData1/v28d_*.tsv` (32 data sheets of the published workbook),
   `CAPSCAN_participant15_join.tsv`, `CAPSCAN_eligibility_audit_347.tsv`,
   `CAPSCAN_raw_header_audit.tsv`, `CapScan_specimen_number_discrepancies.tsv`,
   `source_data/Fig5/design_347.tsv`, and the controlled file
   `SData1_PREDICT_withheld_cells.tsv` when `FASTALAKE_CONTROLLED_DIR` is set.
   Tool: `tables/SData1/panels.py:build_v28d`. Change: none (published values). Output:
   in memory.
2. Tool: `panels.py:cohorts_sheet`, `sheet_directory`. Change: the Cohorts and Sheet
   directory sheets are computed from the roster sheets (counts) and the text constants of
   the original producer. Output: in memory.
3. Tool: `panels.py:write_workbook` (port of the original `sheet()` writer). Output:
   `build/SData1/Sample_rosters_v28d.xlsx` and, via `write_analysis_index`,
   `build/SData1/ANALYSIS_INDEX_v28d.tsv` (sha256
   0e10a9bc4d1dd724592a37bac503df5e156dfb67ecf8d1fe79d05b5c8b7be5cd, byte identical to
   the published file). Checked by `check_v28d` against `expected/` (receipt
   `checks/SData1_v28d.json`: IDENTICAL).
4. Tool: `panels.py:reproduce_v28d_capscan`. Change: none; regenerates the four CapScan
   analyses re-indexed in step 5 on the published 331 roster and requires identical rows
   (positive control, receipt `checks/SData1_v28f.json`).
5. Tool: `panels.py:build_v28f_all` (calls `capscan_acquisitions_347`,
   `capscan_populations`, `membership_rows`, `nonmember_rows`, `regional_pairs`). Change:
   CapScan on 347 acquisitions, as listed in README.md, "What changes from v28d to v28f".
   Output: `build/SData1/Sample_rosters_v28f.xlsx`,
   `build/SData1/ANALYSIS_INDEX_v28f.tsv`.
6. Tool: `panels.py:discrepancy_census`, `discrepancies`. Change: none; the census is
   recomputed from the RAW header audit and recounted against the 347 roster (16
   acquisitions, 14 in the roster). Output:
   `build/SData1/CapScan_specimen_number_discrepancies_v28f.tsv` (sha256
   2f4a14503dd71c1d5ce7440949f2406eefdc36efb1ce0910841de22ba83087e6, byte identical to
   the published v28d table, because none of the sixteen added acquisitions is discordant).

The xlsx files are not byte-stable between runs (openpyxl stores the save time). Their
content is fixed by the per-sheet digests in `checks/SData1_v28f.json` (complete build,
controlled cells present), for example CAPSCAN acquisitions 347 rows, full digest
817f890a337869ab...; CAPSCAN membership 7,753 rows, d130ffa3ef660bcf...; Analyses 30 rows,
ae554b8dd75c9b8b.... Build of 2026-10-06: `Sample_rosters_v28f.xlsx` sha256
bb37b7ab0c29a861e267798589baaeab4ad14f71fa55d25ccf22bf1dcb081362,
`ANALYSIS_INDEX_v28f.tsv` sha256
cd70a0c35202f7231392e6a01fb022d2fad905386435116c66e18e98da34957f. There is no PNG or PDF
rendering of this unit.
