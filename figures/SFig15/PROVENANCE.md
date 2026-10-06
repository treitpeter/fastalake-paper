# Provenance: Supplementary Figure 15

`LAB_ROOT/.../` abbreviates
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/`.
`LAB_ROOT/REV/` abbreviates
`LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/supplements/revision_2026-09-24/`.

## Panels

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a-d, before title removal (v28f, 347) | `LAB_ROOT/.../final/suppfigs_n15/SF15/build_sf15.py` (sha256 0b9a2609...; figure block of `LAB_ROOT/REV/build_revision.py`, sha256 e636bd85..., lines 13, 34-36, 51, 59, 229-237, x tick labels derived from the table), run by `LAB_ROOT/.../final/suppfigs_n15/SF15/run.sh` (sha256 a3fc1397...) | `LAB_ROOT/.../final/suppfigs_n15/SF15/new/COMPUTE_RESOURCE_MAP.pdf` | 82177b1587cf86b1d90b70205d3b447bdf03a1faa677a69f7b624830457a61dd | not recorded (run.log has no job id) | CPython 3.12.12, matplotlib 3.10.8, numpy 2.1.1, pandas 2.2.2 |
| composite, title removed (v28f, 347) | `LAB_ROOT/.../final/suppfigs_n15/SF15/redact_sf15.py` (sha256 368f9188...; logic of `LAB_ROOT/00_PACKAGE_2026-09-17/.../reappraisal_2026-09-24/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/prepare_supplement_figures.py`, copy sha256 5c52b4b8...) | `LAB_ROOT/.../final/suppfigs_n15/SF15/new/Supplementary_Figure_15.pdf` and `.png` | pdf 02d9c5ee79cf5e33a87912ec3e6b6612ffbfe01a65dbd67b523b5a91c6a11fdb; png b343e653aa38c78d90f6753f25225eec8ad4077f5ff2b2b6d85d31dfc286fc56 | not recorded | CPython 3.12.4, PyMuPDF 1.27.1 |
| v28d figure (331, for reference) | `LAB_ROOT/REV/build_revision.py` then `prepare_supplement_figures.py` | v28d `word/media/image23.png` (delivered pre-redaction PDF `LAB_ROOT/REV/delivery/figures/COMPUTE_RESOURCE_MAP.pdf`, sha256 01f0dc45...) | da91ed373dd1ae314d899799750a40e288fbd34153840cdb82bdebfca142de16 | not recorded | as above |
| v30 composite, house style (all panels) | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig15/build_sf15_style.py` (sha256 3d7476c4...) with `sf_common.py` (sha256 9a1c9929...), run by `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/run_sf456.sh` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig15/SupplementaryFigure15.pdf` and `.png` | pdf aec64c675e39d1561ab19e13196ad65e35520c9d5945786fa293f013fdf7661b; png fcf6a785b93e5d024844c1a5aec8bc2b142b5390d46e1dd1576f8a258000e469 | not recorded (build.log has no job id) | the repository environment `fastalake-paper_v28_work/env` (run_sf456.sh): CPython 3.12.4, matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3, PyMuPDF 1.27.1 |

v30 result: `panels.py:figure_v30` in the repository environment gives a PNG byte-identical to the v30 image (`checks/SFig15_v30.json`), with every plotted value identical to the producer drawing. The v30 build also wrote `values_check.json` and `pdf_audit.json` (fonts Söhne Buch and Halbfett only, text ink #4e4f51 or white, clip check passed); they are not copied here.

Repository result: with the repository environment (CPython 3.12.4, matplotlib 3.10.7, numpy
2.2.6, pandas 2.3.3, PyMuPDF 1.27.1) `panels.py:draw_figure` gives a PDF pixel-identical at 300
dpi to the producer's `COMPUTE_RESOURCE_MAP.pdf` (0 differing pixels), `panels.py:redact_title`
gives a PDF pixel-identical at 300 dpi to `Supplementary_Figure_15.pdf` and a PNG byte-identical
to `Supplementary_Figure_15.png` (receipts in `checks/`). The producer itself had shown that
`build_sf15.py` with the 331 table reproduces the v28d image byte for byte (sha256 da91ed37...,
`LAB_ROOT/.../final/suppfigs_n15/SF15/run.log`). The notebook also recomputes the stage summary
from the task records with the producer's summary loop and gets a byte-identical table. Notebook
run on Slurm: 4 s wall time, 0.40 GB peak RSS, 4 cores requested.

Changes between the v28d (331) and v28f (347) figure, all in the CapScan column of the
per-acquisition stages (shared 368-acquisition stages and other cohorts unchanged, checked by the
upstream gate `collect.py`): x label CapScan 331 to 347; a, shuffle search 0.83 to 0.85; c,
acquisition lookup 70 to 74, parsimony 57 to 60, Sage 59 to 62, shuffle search 74 to 78 CPU-h;
d, 331/331 to 347/347, 293/331 to 307/347, 248/331 to 262/347, 271/331 to 285/347; b, printed
values unchanged.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| v28d Supplementary Figure 15 | `word/media/image23.png`, sha256 da91ed37... | BYTE_MATCH; producer MATCHING_REPRODUCER `SF15/redact_sf15.py:16` via `run.sh:15-16`; selected origin `SF15/repro/Supplementary_Figure_15.png`; limit: "original historical renderer not independently located" | original renderer is `LAB_ROOT/REV/build_revision.py` (figure COMPUTE_RESOURCE_MAP, lines 229-237) followed by `prepare_supplement_figures.py`; `build_sf15.py` copies that block and reproduces the image byte for byte, and its pre-redaction PDF text equals the delivered `LAB_ROOT/REV/delivery/figures/COMPUTE_RESOURCE_MAP.pdf` | partly | same image and reproducer; this unit names the historical renderer Codex did not locate. Codex also lists byte-identical copies under other figure numbers (release candidate `Supplementary_Figure_11.png`, co-author package `Supplementary_Figure_13.png`), which are renumbered copies of the same image |
| v28f candidate PNG | REBUILT_CANDIDATES.tsv, sha256 b343e653... | REBUILT_CANDIDATE, write route `SF15/redact_sf15.py:16`, WRITE_ROUTE_INSPECTED_NOT_RERUN | `redact_sf15.py` via `run.sh`; reproduced byte-identically by this repository | yes | it is the v28f replacement |
| v28f candidate PDF | REBUILT_CANDIDATES.tsv, sha256 02d9c5ee... | REBUILT_CANDIDATE, write route `SF15/redact_sf15.py:16` | as above; reproduced pixel-identically | yes | |

Codex number recheck (`out_recheck/REPORT.md`, "Resource 331 -> 347: supported"): 347 distinct
completed CapScan acquisition-lookup keys in the task records. This unit agrees (347 tasks in each
CapScan per-acquisition stage). The Codex statistics register has no row for v28d lines 193,
299, 469 or 903.

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| `source_data/SFig15/compute_stage_summary.tsv` | `LAB_ROOT/.../final/reruns_347/A7_other/resource_accounting/fixed347/source_data/compute_stage_summary.tsv` | 37c1d94d821d5ba9defb5ff5c8495905ce93209a136952dec8b600b74816498f | copied unchanged (`figures/SFig15/prepare_source_data.py`) |
| `source_data/SFig15/compute_task_records.tsv` | `LAB_ROOT/.../final/reruns_347/A7_other/resource_accounting/fixed347/source_data/compute_task_records.tsv` | ea24e72fe85b7e4c0755a849d5822a5347d9d3685e73f765a641072c2d7551b8 | column `sample` (raw acquisition file names) replaced by pseudonymous `acquisition` (`<cohort>_A<rank>`); all other fields byte-identical (`prepare_source_data.py:pseudonymise_records`) |
| `source_data/SFig15/campi_alphanovo_inference_sacct.tsv` | `LAB_ROOT/REV/source_data/compute/INFERENCE_SACCT.tsv` (copy of `LAB_ROOT/Hartout_edits/H10_compute_provenance/results_v2/INFERENCE_SACCT.tsv`, written by `frozen/h10_compute_provenance_v2.py` line 47, sha256 f8e9967f...) | 74f88146f2052441da1e227b4b4521ed79bc95bb16eb7ff67c2cd8be98cc1813 | header from the `--format` of that sacct call added, `|` to tab; values unchanged |
| `source_data/Fig1/acquisition_sizes_726.tsv` (owned by Fig1, read only) | see `figures/Fig1/PROVENANCE.md` | fd9ade3bb4ebd2ee6c129f69819b9697cf20abd7c10c77a12644bef25396e386 | none |

Comparison with the Source Data package (`FastaLake_Source_Data/`): the package folder
`Supplementary_Figure_15/` holds only `compute_stage_summary.tsv` with sha256 c7289c38..., which
is the 331 version (byte-identical to `reruns_347/A7_other/resource_accounting/baseline331/`).
Its non-CapScan rows and its shared CapScan rows equal the 347 table; the CapScan per-acquisition
rows differ as listed above. The package has no task-record table, so the summary cannot be
recomputed from it. `SOURCE_DATA_MAP_reanchored.tsv` rows SF15, SF15-scope and SF15-n point to
that 331 file (SF15-n quotes "331 CapScan"). The map has no row for the numbers on v28d lines 193
and 299.

Upstream of the task records (not in this repository): `LAB_ROOT/REV/source_data/compute/original_compute_jobs.tsv`
(sha256 db00c281...; sacct collection by `original_collect_sacct.py`, sha256 096eff50...),
`current_CAMPI_task_records.tsv` (342e3bfc...), the array rosters `capscan_array_roster.tsv`
(584d49cd...), `predict_array_roster.tsv` (e89aaeee...), `hstool_sample_list.tsv` (20790396...),
`hstool_original_design.tsv` (0244ceed...), and the CapScan 347 roster
`LAB_ROOT/.../coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/roster/fixed_roster_347.tsv`
(473ff6af...). The Computing block (`build_revision.py` lines 201-228, with two substitutions for
the 347 roster, `fixed347/block_executed.py`, sha256 38875803...) was run by
`LAB_ROOT/.../final/reruns_347/A7_other/resource_accounting/scripts/compute_block.py` (4c03a3d4...)
as Slurm job 6064882 (CPython 3.12.4, pandas 2.3.3, numpy 2.2.6) and gated by `collect.py`
(08d7a41d...) as job 6064883, with a negative control that failed as intended.
