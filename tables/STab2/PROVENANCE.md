# Provenance: Supplementary Table 2

## Panels

A table unit has no panels. One row per published artefact or number.

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| workbook, sheets Read_me .. Control_discoveries | `../BiorXiv/99_archive/2026-09-11_code_exports_before_edit6.1/08_code/analyses/full_route_entrapment_2026-09-06/build_package.py` (sha256 3394e0541805e71f519a9e1525802bdfe5ed18de851f56c1f297646280373bff) | `../BiorXiv/99_archive/2026-09-11_before_weekend_organisation/03_tables/supplementary/full_route_entrapment_2026-09-06/SupplementaryTable2_full_route_entrapment.xlsx` | 6b437b2ebcbb7d77f28ed578f98e405381d15c00a1a9fb37e04eccc347a4b42e | 5877268 (publication_build) | OPEN (not recorded); rebuild: Python 3.12.4, pandas 2.3.3, openpyxl 3.1.5 |
| workbook with fusion_ and indep_ sheets | `../BiorXiv/08_code/edit6.1/figures/upstream_supplement_layout/package_edit3.py` (sha256 8b2029987461d3a22b9d998f4ebcc825772b44a809e7f4a26d2eca9e5ec47bd2) | `00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/supplements/revision_2026-09-24/delivery/source_data/tables/CONTROL_EXPOSURE/SupplementaryTable2_full_route_entrapment.xlsx` | 2d33ce310b38b0a04b1201edaf40cf1bcc193f22fb40f5dcd321321ff024c843 | none (login-node script, no job id recorded) | OPEN |
| workbook renamed and path-scrubbed | `00_PACKAGE_2026-09-17/.../codex_review_2026-09-25/v15_MM_coauthor_2026-09-25/package_v15.py` (sha256 0573a006ad4c627a483853fda0052772af26a1dc4fe94007cd94569a143a33fd) | `.../v15_MM_coauthor_2026-09-25/coauthor_package/supplementary_tables/Supplementary_Table_2.xlsx` | 5eb86111a0eacf5983c3f750d302281fba83f6c5c44700ce10b82efa09b1027d | none recorded | OPEN |
| workbook as published (expected/) | `review/v15_integration/publish_ready_v16.py` (sha256 7ba976fa58a6635948eba83e8b0b449ad2e68a3a16028aa5e18046f868d36f8f); `.../final/figure_sd/source_data/build/selection.py` (sha256 20e4fabb92433db332e6f426cbb5245f3246157776c0b8d565181f115bea796b) | `review/FastaLake_v16_SOURCE_DATA.zip` member `current/Supplementary_Tables/Supplementary_Table_2.xlsx`; `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Tables/Supplementary_Table_2.xlsx` | 4768d39fff4b3fa3605a5cc5af504b084356a98db232f52a2cb8f246dcedbc35 | none recorded | OPEN |
| v28d line 743 counts | `FIGURATOR/Figure3/review_2026-09-07/entrapment_pooling/reduce.py` (sha256 8adff074b99f9c7dc6e8652214097487748d3bdf51779e1adb9d6ff1b9286b14) | `FIGURATOR/Figure3/review_2026-09-07/entrapment_pooling/per_acquisition.tsv` | 474de2e78bbfc37fcf132930b688069df06909699b6f7cb9a0c12bf0c42222f8 | 5879093 (search array), 5879104 (reduce) | user miniconda python, version not recorded: OPEN |

Earlier independent recount (not reused): `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/work/gap_unchanged_figures_supp_tables/check_st2_pooling.py` printed 2/376,412 (own) and 12/411,354 (pool_other11) from the same table and found neither denominator in the workbook (`check_st2_st5.log`). This unit agrees with both results.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| Supplementary Table 2 workbook | none | not in census (the census covers the 23 images embedded in the docx) | `build_package.py` -> `package_edit3.py` -> `package_v15.py` -> v16 archive copy | NA | No Codex statistics-register row or number-recheck row refers to v28d line 477 or 743. |

## Source data

Every origin hash of the four full-route cohorts equals the hash recorded in that cohort's `COMPLETE.json` (24 of 24 files checked).


| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| per_acquisition.tsv | `canonical/publication/analyses/full_route_entrapment_2026-09-06/tables/<cohort>/per_acquisition.tsv`; package sheet Supplementary_Table_2.xlsx!Per_acquisition (values identical) | CAMPI_MOCK 7fbbd87610c1a5298f80a7e3a004470e2cda563da260d4f6d8e0e2a8916f761d; HSTOOL 8bd413c252992de908822dceae425229eb18488dac01e4b89939e524cd4541d5; PREDICT 378e2e585174abb1de0d0de4f41679c7276ab6bb47e6090ca5ba90f6677a1af6; CAPSCAN c182e2b1cce6748a8479762f7737fc10eb52f850533d91f27ae6d46d20e434c9 | line concatenation in cohort order, one header (`tables/STab2/prepare_source_data.py`) |
| stage_control_counts.tsv | same directories, `stage_control_counts.tsv`; sheet Stage_controls (identical) | CAMPI_MOCK fbce9906dd3d2d6570b9f29adf3a602bdd7f11eb43637b7c4a7aba622b60f925; HSTOOL 61453e81dfee1380eaf1294494524dbaf804fbe328fd92af2e71295895922a7f; PREDICT 49ab806664aa642aac80ca1df5bb18a9d374a696c17396dc599158b07072b203; CAPSCAN 893029d1a66256745d68af3288f9bcc9cd7ceda0cebe78fa1f970f6d9d0e2922 | line concatenation |
| baseline_comparison.tsv | same directories, `baseline_comparison.tsv`; sheet Baseline_comparison (identical) | CAMPI_MOCK f963ed20745e7564117349001b7e5990f3e9b282218f1e9dfdad6fc87ac855f7; HSTOOL 73a0d370fd8a55e2547a39f2b15d6c5506fec934f87d2466c5e7ff47032a0d91; PREDICT ad00f14eabdfbe7f7a9b9b1d57bc44747f092ec90b946632b66b2330a7f4a14d; CAPSCAN eacc2fa2950f5671898588699f56cbf06d002fd1da19cacdefe828016238dbc8 | line concatenation |
| quantile_summary.tsv | same directories, `summary.tsv`; sheet Quantile_summary (identical) | CAMPI_MOCK 233f8509033ec03ecebfc29b3bbcd77aeafb98274637395c1747cb4b5a11cc56; HSTOOL 90aa4037a641b85bc578de640fb21033fa1adb04424a5f1fbbef9a703310b2d8; PREDICT fcfe131ec3c57d5a8d08bf3e10b59969e88ee79ef944a40c4806b3518203a92e; CAPSCAN a27905b90a74f0619b91d0a8395247c87b8caa6508bfd84bcb95106217d74ef3 | line concatenation |
| control_discoveries_<cohort>.tsv | same directory, `control_discoveries.tsv`; sheet Control_discoveries (identical) | equal to the copies (sha256 in source_data/STab2/README.md) | byte copy |
| fusion_*.tsv, fusion_validation.json | `../BiorXiv/99_archive/2026-09-11_before_weekend_organisation/review_2026-09-08_v10_AP/source_data/completed_fusion/` (archived inputs of `SourceData_fusion_completed.xlsx`); sheets fusion_* (identical) | equal to the copies (byte copy) | byte copy |
| indep_*.tsv, indep_validation.json | `.../review_2026-09-08_v10_AP/source_data/campi_independent_selection/`; sheets indep_* (identical) | f0aa932e... and 6c912a08..., equal to the hashes recorded in VALIDATION.json | byte copy |
| hstool_pooling_per_acquisition.tsv | `FIGURATOR/Figure3/review_2026-09-07/entrapment_pooling/per_acquisition.tsv`; package file `Figure3/per_acquisition.tsv` (byte-identical); SOURCE_DATA_MAP_reanchored.tsv row RES-53 | 474de2e78bbfc37fcf132930b688069df06909699b6f7cb9a0c12bf0c42222f8 | byte copy |

SOURCE_DATA_MAP_reanchored.tsv rows for this unit: ST2 (legend, LOCATED, workbook) and RES-53
(line 743, source `Figure3/per_acquisition.tsv`, VERIFIED, flag "anchor not in text"). Both agree
with this unit. The map does not record that the line-743 values are absent from the cited
workbook.
