# Assembly: Supplementary Table 6

Short names `CR` and `SD` as in `PROVENANCE.md`. No manual step.

1. Snapshot. Input: every `NULL_COMPARISON.json` under the paired-null outputs on 25 September 2026 (60 files). Tool: `CR/audit_remaining.py` lines 45-50. Change: stores path and content of each receipt. Output: `CR/NULL_RECEIPTS_SNAPSHOT.json`, sha256 c2fb0267676b7ecc90c291bb503ca416c75a16d1d93a252974c6dc2b4eb581d3.
2. Support tables. Input: step 1. Tool: `CR/package_release.py` lines 20-26 (ported as `panels.py:build_acquisition_seed_rows`, `build_summary`, `write_table`). Change: keeps q = 0.01 and canonical_peptide rows, sums accepted and selected_proteins by route. Outputs: `null_control_by_acquisition_seed.tsv` 3efbcb3f3dd8a19c141c836727bae8a2ee6676a8124b72b3fd580008c173e896; `null_control_summary.tsv` e794df7b7cbf3a57004123ee7b09405089c1f8a4b2d4cd942e1731bcdf720e2c.
3. Workbook. Input: step 2. Tool: `CR/v15_MM_coauthor_2026-09-25/package_v15.py` lines 54-58 (ported as `panels.py:write_table`). Change: sheets Read me (two text rows), Summary, Acquisition_seed_rows, cells written as strings. Output: `Supplementary_Table_5.xlsx`, sha256 03d4b914c380c950ce33b1222d2af07f8582f9c11bb186c0566f2ce3ba14ef66.
4. Archive and renumbering. Input: step 3 as member `current/Supplementary_Tables/Supplementary_Table_5.xlsx` of the Source Data archive of 28 September 2026. Tool: `SD/build/selection.py` (`RENUMBER`: 5 to 6) and `SD/build/build_package.py` (byte-identical copy). Change: file name only. Output: `SD/package/FastaLake_Source_Data/Supplementary_Tables/Supplementary_Table_6.xlsx` = `expected/Supplementary_Table_6.xlsx`, sha256 11f7be117046d8eef2361ea05c2c710bed8728fa4587a96328bd6e9d74529ceb.

v28f uses this workbook unchanged (inventory status "in v28d as-is"). There is no rasterised composite for a table.
