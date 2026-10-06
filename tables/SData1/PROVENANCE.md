# Provenance: Supplementary Data 1

`M01` below abbreviates `LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M01_sample_rosters`.
`CR` abbreviates `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03`.
`C2` abbreviates `CR/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15`.

## Panels

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| workbook (first build, 36 sheets) | `M01/build_workbook.py` (sha256 a54c800a217f1fdec7a103642b463928e7e995f6d3d50ece42f6c84a9873d18c) | `M01/delivery/Sample_rosters.xlsx`; check `M01/WORKBOOK_VALIDATION.json` (PASS, 710 acquisitions, 30 analyses, 94 sources rehashed) | 037fa3b26fb3e2966060549669c636c41154c9ba4e482afa4fc1633eebf32c36 | not recorded (run interactively) | not recorded; openpyxl |
| workbook as archived 28 September 2026 | copy into the Source Data archive | `LAB_ROOT/review/FastaLake_v16_SOURCE_DATA.zip` member `current/Supplementary_Data_1_Sample_rosters.xlsx` | 99d4ca952a151cd667e8c7287bc866fdaf32706f98a0a27b519de4737471a66f (member) | none | none |
| workbook as published with v28d (34 sheets) | `CR/final/figure_sd/source_data/build/build_package.py` (sha256 730ded5c10f3fc6964d554e249478711d45586470696b77bdf053166d0e74a04): PREDICT participant and tube codes to random pseudonyms, rows sorted by code re-sorted by pseudonym (5 PREDICT sheets, ORDER_LOG.tsv), paths reduced to file names, 17,504 cell edits; the Sources and Author confirmations sheets were not in the archived member | `CR/final/figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Data_1/Sample_rosters.xlsx` | 5ee76e28adb77d2edd9d424cf450538f72a8a6dc4f66bd07c091f524bafcff72 | not recorded in this unit | openpyxl 3.1.x |
| ANALYSIS_INDEX.tsv as published | `M01/build_workbook.py` (write_tsv), then `build_package.py` (14 path edits) | `CR/final/figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Data_1/ANALYSIS_INDEX.tsv` (origin `M01/delivery/ANALYSIS_INDEX.tsv`, sha256 dd9325dcb5b69c072fc5a3c327a5bcff9e6b5a62dd49e070c285b55f5f96ab27) | 0e10a9bc4d1dd724592a37bac503df5e156dfb67ecf8d1fe79d05b5c8b7be5cd | not recorded | Python csv |
| specimen-number census as published | `CR/iteration3/source_data/build/derive_verified.py` (sha256 14bc7a66a7d178593e288e2e3af6c6d5fdd31b6f9ff2c666e02507461a060ddb), lines 221-231 | `CR/final/figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Data_1/CapScan_specimen_number_discrepancies.tsv` | 2f4a14503dd71c1d5ce7440949f2406eefdc36efb1ce0910841de22ba83087e6 | not recorded | pandas |
| CapScan 347 roster | `C2/scripts/build_fixed_roster.py` (sha256 7157d95b077ddcdc533144290d91338d3d5db8c301b34f45a6b5b691afd625bf), `C2/scripts/make_design347.py` (sha256 de427cad5581da3f19dc7e87110851b79ef100b94997d3bfc5284259845a1d2f) | `C2/roster/participant15_added.tsv`, `C2/roster/design_347.tsv`, `C2/roster/fixed_roster_347.tsv`, `C2/roster/ROSTER_FIX.json` | 688dea12752c41fb95015bb6350156e5d1d9522511cd201a09546ea3718e559c; 30f70233143a370755d1bc4d1dd8f236d94a2fe4b9b27cc10c64e94602461000; 473ff6afcf7348163f1fbaa8d624a8078eed9348a15019952363bc86e70f7aa5; 3e5cde6236a7b831ea47a7fa9408db610c2f05ab49c8bf88f69d743f940bfff6 | not recorded (C2/JOBS.txt lists no roster job) | pandas |
| CapScan 347 nonempty flags | `C2/fig5g/crossomics_capscan.py` (sha256 0ab9a2465b4308b145dcc75d81423601a5fc50ebe4ef16f5aa9e02b79f9591d2) | `C2/fig5g/crossomics_347/eligibility_audit.tsv` | 3c58ea360554a7631bdef00c21fad5b394ebd0b4fbfab72b6ff537e585ceb36b | 6063849 | pandas |
| CapScan 347 paired-region run (cited) | `C2/scripts/capscan_fig5bcd.py` (sha256 901e97146c2491775574f8eff9bb74daa8a78e1ff6f5e63f2a23d072ac6106b2) | `C2/fig5/fig5bcd_347/own_stool_participants.tsv`; `C2/fig5/fig5bcd_347/COMPLETE.json` (346 acquisitions, 15 participants, 75 cells) | 67f22f573d418bdee6e7b2a7ed0bc1ccc2fc140ee5793ab10db2455f6002ff0e; e62cda8f6d68b1a0ee7e86350d410e3f107cda92c3cbd10368133943007eb0ed | 6063597 (profiles), 6063598 | pandas |
| RAW header audit | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/cohort_reconciliation/raw_identity_audit/audit_headers.py` (sha256 3012a0fc5e072d4ebaa1f2300035d651cf6892afbcd2682507a9fdb92207be22) | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/cohort_reconciliation/raw_identity_audit/all_RAW_header_names.tsv` (368 RAW files) | 14bc80a345186d7bb31336cd8cc3f574ea7b98d12bd8d1d526e02a06ddacd97f | not recorded in this unit | not recorded |

This unit's build: `tables/SData1/SData1.ipynb`, executed with `python paper.py build SData1`
under `srun -p p.hpcl92,p.hpcl8,p.hpcl93 --nice=10000 -c4 --mem=4G` (wall 12 s, peak RSS of
the runner 171 MB), CPython 3.12.4, openpyxl 3.1.5, pandas 2.3.3 (environment
`environment/requirements.lock.txt`), with the controlled PREDICT cells present.

Result (`checks/SData1_v28d.json`): all 34 rebuilt sheets IDENTICAL to the published
workbook (full-cell digests including the withheld cells), sheet order identical, the
Cohorts and Sheet directory sheets identical cell by cell, `ANALYSIS_INDEX.tsv` byte
identical. The census table recomputed from the RAW header audit is identical to the
published one. Positive controls (`checks/SData1_v28f.json`, `positive_controls_on_331`):
the four CapScan analyses re-indexed for v28f and the regional pairs, regenerated on the
published 331 roster, equal the published rows.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| Supplementary Data 1 workbook | none | not in scope (the census covers the 23 images embedded in the manuscript; Supplementary Data 1 is a workbook, not an image) | as above | not applicable | Codex number recheck (`out_recheck/REPORT.md`, rows N0134, N0139) independently moves the cohort total from 710 to 726 acquisitions for the split-spectrum experiment; this unit computes 726 from the v28f rosters, in agreement. The Codex statistics register has no row for this unit (no tests). |

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| `v28d_<sheet>.tsv` (32 files) | Source Data package `Supplementary_Data_1/Sample_rosters.xlsx`, sheet `<sheet>` | 5ee76e28adb77d2edd9d424cf450538f72a8a6dc4f66bd07c091f524bafcff72 | `tables/SData1/prepare_source_data.py`: cell values as read by openpyxl, unchanged; PREDICT columns in `panels.WITHHELD_COLUMNS` and rows in `panels.WITHHELD_ROWS` moved to the controlled file. Values agree with the package: the notebook rebuild reproduces every sheet digest. |
| `CapScan_specimen_number_discrepancies.tsv` | package `Supplementary_Data_1/CapScan_specimen_number_discrepancies.tsv` | 2f4a14503dd71c1d5ce7440949f2406eefdc36efb1ce0910841de22ba83087e6 | byte copy; recomputed identically by `panels.discrepancy_census` |
| `CAPSCAN_raw_header_audit.tsv` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/cohort_reconciliation/raw_identity_audit/all_RAW_header_names.tsv` | 14bc80a345186d7bb31336cd8cc3f574ea7b98d12bd8d1d526e02a06ddacd97f | six of 16 columns kept (`prepare_source_data.RAW_AUDIT_COLUMNS`); local file paths and checksums dropped. Not in the Source Data package. |
| `CAPSCAN_participant15_join.tsv` | `C2/roster/participant15_added.tsv` | 688dea12752c41fb95015bb6350156e5d1d9522511cd201a09546ea3718e559c | byte copy. Not in the Source Data package (the package predates the 347 roster). |
| `CAPSCAN_eligibility_audit_347.tsv` | `C2/fig5g/crossomics_347/eligibility_audit.tsv` | 3c58ea360554a7631bdef00c21fad5b394ebd0b4fbfab72b6ff537e585ceb36b | byte copy; agrees with `og_profile_positive` of `source_data/Fig5/design_347.tsv` for all 347 acquisitions (asserted by `panels.build_v28f`). Not in the Source Data package. |
| controlled `SData1_PREDICT_withheld_cells.tsv` (not in the repository) | package `Supplementary_Data_1/Sample_rosters.xlsx`, PREDICT sheets | 5ee76e28adb77d2edd9d424cf450538f72a8a6dc4f66bd07c091f524bafcff72 | `prepare_source_data.py`; 5,254 cells; file sha256 1d52de25c42605a8d8ae90135050921e838131a8ac5cd312f5e2302bc869cfff |

`SOURCE_DATA_MAP_reanchored.tsv` rows SD1 and M-ROSTER point to the workbook and
`ANALYSIS_INDEX.tsv` (LOCATED); rows RES-CS16 and RES-CS14 give 16 and 14 from the census
table (VERIFIED). Both agree with this unit. The map still describes the 710-acquisition
(331 CapScan) workbook; the v28f workbook built here replaces it.

## v30 build (6 October 2026)

| source_data file | derived from (LAB_ROOT/...) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| `v30_capscan_crossomics_specimens.tsv`, `v30_capscan_crossomics_participants.tsv`, `v30_capscan_metabolite_specimens.tsv` | 347 cross-omics run of capscan_crossomics (Slurm 6063849): `capscan_specimen_roster.tsv`, `COMPLETE.json`; extended metabolite table | `source_data/SData1/v30_inputs_receipt.tsv` | `tables/SData1/prepare_v30_inputs.py:main` (Slurm 6075320) |
| `v30_capscan_h_design.tsv`, `v30_capscan_h_finite_n.tsv`, `v30_capscan_h_participants.tsv` | fixed-347 A28 histidine run (Slurm 6064031): `inputs/manifest_revision4.json`, `complete_statistical_families.tsv`, `participant_region_enrichment.tsv`, `participant_reciprocal_controls.tsv` | `source_data/SData1/v30_inputs_receipt.tsv` | `tables/SData1/prepare_v30_inputs.py:summarise_histidine` |

The v30 public workbook (`build/SData1/Sample_rosters_v30.xlsx`) replaced `Supplementary_Data_1/Sample_rosters.xlsx` and
`ANALYSIS_INDEX.tsv` in the Source Data package on 6 October 2026 (package MANIFEST.tsv and SHA256SUMS.txt updated,
`sha256sum -c` clean; Slurm 6075392). Receipt: `checks/SData1_v30.json`.
