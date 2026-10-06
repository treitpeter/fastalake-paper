# Provenance: Supplementary Table 6

Short names used below:
`CR` = `LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/johannes_main_text_review_2026-09-24/codex_review_2026-09-25`;
`PN` = `LAB_ROOT/MM_review_analysis/A04_4_2_expansion/alphanovo_bags_2026-09-23/paired_nulls_2026-09-23`;
`SD` = `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/figure_sd/source_data`.

## Panels

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| per-acquisition receipts (all sheets) | `PN/runtime/tools/run_bag_nulls.py` (sha256 8677f2b54976df6c9bd09fa82612856bc81db909416b930e64df9c4a81f9b871), function `report`, called by `PN/guard.py report` (a66555bca8c0f23a42289f55abf514f7d7d742626e3044372fbe14b5302430ea) under manifest `PN/manifest.json` (6f02c2677bb6b0b91b58cf41542fc154174877b452de7d2cd928382618a5a0c8) | `PN/outputs/<cohort>/<acquisition>/NULL_COMPARISON.json` (60 used) | per-receipt Sage result hashes in `source_data/STab6/null_receipt_results.tsv` | pilot report 5993697, report array 5993701 (from `PN/SUBMITTED.json`) | CPython 3.13.11 (from `PN/jobs/report.sbatch`); Sage binary sha256 d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625 |
| frozen snapshot | `CR/audit_remaining.py` lines 45-50 (eefb21149c09152b6448a5c7f72c8086063645583e5b91962c197dc88c8f455e): globs every `NULL_COMPARISON.json` present under `PN/outputs` at run time and stores path and content | `CR/NULL_RECEIPTS_SNAPSHOT.json` (60 receipts, written 25 September 2026 15:23) | c2fb0267676b7ecc90c291bb503ca416c75a16d1d93a252974c6dc2b4eb581d3 | not recorded | not recorded |
| Summary and Acquisition_seed_rows TSVs | `CR/package_release.py` lines 20-26 (17ff4cbbe7384650d2300f26b69b5d9798fa1010d5b5e20722815785592d0500) | `CR/release_candidate/supplementary_additions/null_control_summary.tsv`, `.../null_control_by_acquisition_seed.tsv` | e794df7b7cbf3a57004123ee7b09405089c1f8a4b2d4cd942e1731bcdf720e2c; 3efbcb3f3dd8a19c141c836727bae8a2ee6676a8124b72b3fd580008c173e896 | not recorded | pandas (version not recorded) |
| workbook (then Supplementary Table 5) | `CR/v15_MM_coauthor_2026-09-25/package_v15.py` lines 54-58 (0573a006ad4c627a483853fda0052772af26a1dc4fe94007cd94569a143a33fd) | `CR/v15_MM_coauthor_2026-09-25/coauthor_package/supplementary_tables/Supplementary_Table_5.xlsx` | 03d4b914c380c950ce33b1222d2af07f8582f9c11bb186c0566f2ce3ba14ef66 | not recorded | openpyxl (version not recorded) |
| workbook as published (renumbered 5 to 6) | `SD/build/selection.py` (`RENUMBER`, 20e4fabb92433db332e6f426cbb5245f3246157776c0b8d565181f115bea796b) and `SD/build/build_package.py` | `SD/package/FastaLake_Source_Data/Supplementary_Tables/Supplementary_Table_6.xlsx` = `expected/Supplementary_Table_6.xlsx` | 11f7be117046d8eef2361ea05c2c710bed8728fa4587a96328bd6e9d74529ceb | not recorded | package `MANIFEST.tsv`: byte-identical copy of member `current/Supplementary_Tables/Supplementary_Table_5.xlsx` of the Source Data archive of 28 September 2026 |

The coauthor-package workbook (03d4b914...) and the published workbook (11f7be11...) differ
in bytes but not in any cell value; the published one is the copy carried through the
28 September archive.

Reproduction in this repository: the notebook rebuilds both TSVs byte-identical to the
published copies and the workbook identical in every cell (`checks/table.json`, result
`IDENTICAL`; Slurm job 6069116, environment `<work>/env`, Python 3.12.4, pandas 2.3.3,
openpyxl 3.1.5). Earlier independent reproduction:
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/work/gap_unchanged_figures_supp_tables/check_st2_st5.log`.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| Supplementary Table 6 | none | not covered (the census covers the 23 embedded images; this table has no image) | as above | NA | Codex statistics register: line 485 classed as descriptive count with no inferential P or interval (`STATISTICAL_TEXT_COVERAGE.tsv`); no register row for line 561's null-control sentence. Codex number recheck (`out_recheck/REPORT.md`): no row for this table. Nothing to reconcile. |

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| `null_receipt_rows.tsv` | `CR/NULL_RECEIPTS_SNAPSHOT.json`, all 80 rows of each of the 60 entries | c2fb0267676b7ecc90c291bb503ca416c75a16d1d93a252974c6dc2b4eb581d3 | `<work>/tables/STab6/extract.py`: adds `receipt` and `cohort` from the receipt path, drops the path, keeps every row and field (Slurm job 6069108). Package comparison: filtered to q = 0.01 and canonical_peptide it reproduces package `Supplementary_Tables/Table6_support/null_control_by_acquisition_seed.tsv` and sheet `Acquisition_seed_rows` exactly (values agree). |
| `null_receipts.tsv` | same snapshot, receipt-level fields | as above | same script. Not in the package. |
| `null_receipt_results.tsv` | same snapshot, `results` field (Sage result bytes and sha256 per arm) | as above | same script; absolute paths dropped. Not in the package. |

`SOURCE_DATA_MAP_reanchored.tsv` rows ST5, ST5b and ST5-60 (legend line 485, values 13518,
561331 and 60, marked VERIFIED) point to the same package files; the values agree with
this unit's recomputation.
