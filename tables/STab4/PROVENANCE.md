# Provenance: Supplementary Table 4

In this file `CR` abbreviates `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03`,
`IT2` abbreviates `CR/iteration2/cr003_followup`, `IT1` abbreviates
`CR/iteration1/fix_eggnog_native/downstream/verified50/B` and `PKG` abbreviates
`CR/final/figure_sd/source_data/package/FastaLake_Source_Data`.

## Panels
| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| workbook, all data sheets | IT2/scripts/tables_cr003.py (sha256 7453d63a10a4cb6f53f52d850c14fd51bcb009f4f263d5316b57b6eb43439f94), variant verified50, section "Supplementary Table 3" (renumbered 4 in v28d) | IT2/tables/verified50/Supplementary_Table_3.xlsx | 46ea54c41fb867a6810e672a391bffd0804aeb190f62ce4e05c48d7b554ccd42 | 6048665 (cr3f_tables_verified50, 2026-10-03 14:26-14:27) | CPython 3.12.4 (miniconda), numpy 2.2.6, pandas 2.3.3, openpyxl 3.1.5 |
| inputs of the Focused_* sheets | IT2/scripts/build_approved_eds_cr003.py (sha256 0917829d274e9cc06b7f2715198ebabcc0f69a2e2cae0a5415b7219d49b871fc), run as IT2/shadow/verified50/eds/build_approved_eds.py | IT2/shadow/verified50/eds/EDS_approved/source_data/ (focused_raw_measurements.tsv afb8db0d..., focused_features.tsv 7945c40a..., selected_pair_summary.tsv 42e89ff3..., separate_adjustments_all_models.tsv fe72b674...) | per file | 6048652 | CPython 3.12.4, numpy 2.2.6, pandas 2.3.3, scipy 1.16.3 |
| inputs of the 20-feature sheets | IT1/eds_biochemistry.py with IT1/eds_refit/analyze.py | IT1/eds_refit/input_data/heatmap_values.npz (b937bb50...), heatmap_85_values.npz (11c17120...), heatmap_columns.tsv (75101454...), heatmap_features.tsv (b883e165...), heatmap_85_patient_order.tsv (43dc72c4...) | per file | 6048399 | CPython 3.12.4, scipy 1.16.3 |
| published workbook | CR/final/figure_sd/source_data/build/build_package.py (sha256 730ded5c10f3fc6964d554e249478711d45586470696b77bdf053166d0e74a04) with specials.py (4b5d700bd04089af4682f851a310f77b882a5d1ed166be19703ecfa9191a36d4, function st3) and selection.py (20e4fabb92433db332e6f426cbb5245f3246157776c0b8d565181f115bea796b) | PKG/Supplementary_Tables/Supplementary_Table_4.xlsx | 4fe89a55c7607b0e5932377cbd488c031bc1f0c7f99c255f6c9eb7d77747e403 | MANIFEST_build.tsv (CR/final/figure_sd/source_data/work/full) was written by job 6055130 (fsd_sd); the xlsx file time (2026-10-04 05:55) predates it, so an earlier build wrote the same bytes; that job id is not recorded | CPython 3.12.4, openpyxl 3.1.5 |
| five-pair decomposition (line 799; not in the workbook) | IT1/eds_refit/separate_adjustments/run_analysis.py (sha256 e27b62da1d938229a6cea7688d7ecbf335417fe0c3b62c16b533bb08832597a3), argument full | IT1/eds_refit/separate_adjustments/all_models.tsv = PKG/Biology_refit/eds_refit/separate_adjustments/all_models.tsv | 1f6400e71348f3d87d5248f407f7b9f8f650f4334e9a5eddf21938d315ea1028 | 6048399 (COMPLETE.json job_id) | CPython 3.12.4, numpy 2.2.6, scipy 1.16.3, statsmodels |

The 20-feature inputs read by tables_cr003.py (IT1/eds_refit/input_data/heatmap_values.npz
and heatmap_85_values.npz) are byte-identical to the copies the Supplementary Fig. 10
unit derived its tables from (EDS_approved/source_data/legacy_259_heatmap_values.npz and
legacy_85_values.npz), so the two units share `source_data/SFig10/`.

Checks recorded by this unit:

- Published workbook against producer workbook, cell for cell
  (`prepare_source_data.py:compare_package_producer`): every data sheet equal except
  patient_id (256 of 259 cells in Observation_metadata, 85 in 85_patient_order, 85 in
  Focused_raw; the 3 equal cells are empty in both), which the package builder
  pseudonymised. The builder's REDACTION_LOG.tsv records 426 participant codes and 1
  path edit for this file, which agrees (256 + 85 + 85 = 426; the path edit is the last
  Read_me line, "Biology_refit/eds_refit" shortened to "eds_refit"). The producer's
  Current_EDS_update and Read_me sheets were replaced by one Read_me sheet.
- Rebuilt workbook against the redacted published workbook: 0 differing cells in all 11
  sheets (`checks/workbook_vs_published_redacted.tsv`).
- Rebuilt complete layout (controlled input) against the published workbook: 0 differing
  cells in all 11 sheets (`checks/workbook_vs_published_full.tsv`), run 2026-10-06 on
  Slurm with the pseudonymised controlled copies listed in OPEN_ITEMS.md.
- The earlier replay control of the producer (IT2/old_vs_new/per_variant/TABLES_replay.json)
  showed that the replay variant equals the frozen pre-fix tables; the verified50 variant
  used here differs from them by design (native annotation check, CR-003).
- Statistics recomputed from source_data: unadjusted pair rho and n, family q values
  (largest difference to the saved screen q over all 11,496 screen rows: 0.0), EDS
  associations of the six features, all 192 adjustment estimates (largest difference
  0.0), joint-adjusted rho, Fixed_feature_key rho and n, z scores (largest difference
  1.3e-15).

The repository notebook ran with the same interpreter and package versions
(`environment/requirements.lock.txt`).

Float parsing: the producer read its own TSV files with the pandas default float parser
before writing them into the workbook; that parser can differ from the exact value in the
last binary digit (3 cells in Fixed_feature_key, 1 in Focused_raw, 169 in
Focused_adjustments). `panels.py` reads those tables the same way for the workbook copy
(`_read(..., producer_parse=True)`) and uses the exact values for all statistics.

## Codex census reconciliation
| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| Supplementary Table 4 (workbook) | none | the census covers the 23 images embedded in the v28d docx; it has no row for supplementary tables | producer IT2/scripts/tables_cr003.py, packaged by build_package.py | NA | REBUILT_CANDIDATES.tsv lists no alternative for this table. The Codex number recheck (out_recheck/REPORT.md) has no row for this table. The Codex statistics register rows S063-S075 cover the statistics; their gaps are filled in STATS.tsv and listed in OPEN_ITEMS.md |

## Source data
| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| source_data/STab4/observations_259.tsv | PKG/Supplementary_Tables/Supplementary_Table_4.xlsx, sheet Observation_metadata | 4fe89a55c7607b0e5932377cbd488c031bc1f0c7f99c255f6c9eb7d77747e403 | tables/STab4/prepare_source_data.py: sample_id, microbiome_id, clinical_group, antibiotic_state removed; observation_index added; values agree with the package sheet exactly |
| source_data/STab4/patients_85_observation_index.tsv | same workbook, sheet 85_patient_order, column index | 4fe89a55c7607b0e5932377cbd488c031bc1f0c7f99c255f6c9eb7d77747e403 | index + 1; checked that each index points to the same acquisition in Observation_metadata |
| source_data/STab4/read_me.tsv | same workbook, sheet Read_me | 4fe89a55c7607b0e5932377cbd488c031bc1f0c7f99c255f6c9eb7d77747e403 | one row per line |
| source_data/STab4/five_pair_all_models.tsv | PKG/Biology_refit/eds_refit/separate_adjustments/all_models.tsv (identical to IT1/eds_refit/separate_adjustments/all_models.tsv) | 1f6400e71348f3d87d5248f407f7b9f8f650f4334e9a5eddf21938d315ea1028 | copied |
| source_data/SFig10/* (shared) | see figures/SFig10/PROVENANCE.md | see there | the rebuilt sheets equal the package sheets (0 differing cells) |

`SOURCE_DATA_MAP_reanchored.tsv` row ST3 points the legend of Supplementary Table 4 to
`Supplementary_Tables/Supplementary_Table_4.xlsx`, sheet All20_features_259_observations,
status LOCATED; this unit confirms the location and the values.

Expected file: `expected/Supplementary_Table_4_redacted.xlsx` (sha256
01314a8f757ae7a996caddfe6c4a6fb20339b94296ead88c583c6a4897b975a0) is the published
workbook with the four withheld columns removed and sample identifiers replaced, written
by `prepare_source_data.py:redact` directly from the published file (not from
source_data).
