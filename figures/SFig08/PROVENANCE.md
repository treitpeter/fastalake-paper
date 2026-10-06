# Provenance: Supplementary Figure 8

`P = LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24`
below. Supplementary Fig. 8 of v28d was Supplementary Fig. 6 in the v15 numbering, which
the producer file names keep.

## Panels

v30 (current). `FIN` is `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`.

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a-d, v30 house style | `FIN/v30/figures_style/SFig08/build_sfig08_style.py` (9597973fcac36d3e0290fce6a84a2299c8c8bf12d352607eeeb58f1edeb1e66d) | `FIN/v30/figures_style/SFig08/SupplementaryFigure8.png`; `.pdf` | 4315bca570ac1be5cc223784b8ddd60199f7533be9e59635f6acdbe26f82cd12; ac4e41d6f885100e8cbfdd01ae349d872c26608a7e8b3261dd3109663b2d4ce8 | srun via `FIN/v30/figures_style/_work_sf07_09/run.sh` (job id not recorded) | `fastalake-paper_v28_work/env` python; fl_style with Söhne |
| a at 9-50 (v29) | `FIN/len950/SFig08/build/build.py` (f39f49d5165cb55bef77a6527280856ddfb477779443727e70aedebcaa5e33dc); recount `FIN/len950/SFig08/build/recount_curation.py` (985ab4eb1f5904f8c30e2dea7f33657190f32421b297da7157282090f61753c2) | `FIN/len950/SFig08/SupplementaryFigure8.png`; `FIN/len950/SFig08/curation_counts_5_9.tsv` | 86c8e81db5a4c5f6b71d32d55ac2c04807696a9e31b84994b539ba162425ef6d; faf0e097d3845f94f02caa8ac46a0f01107b2abbd8d6778d27904944f51cb49f | logs `FIN/len950/SFig08/logs/` (job ids not recorded) | as above |

The repository render is IDENTICAL to the v30 PNG and PDF at 300 dpi. v28d producer (superseded, code kept as `panels.py:figure` and used inside `render_v30`):


| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a-d | `P/supplements/revision_2026-09-24/build_revision.py` (sha256 e636bd85164896f53ec1d13ef671754ff800597d2d955a3eb01f5151dced1bbe), block "Clinical comparator scope and fixed-dictionary loss" | `P/supplements/revision_2026-09-24/delivery/figures/CLINICAL_COMPARATOR_SCOPE.pdf` (+ .png) | 10093e473bca0adeb999d59c8bfa2bf644535fba3f627992cdc3b7f861dbe7af (png a862b529688b8f542bd40d195b257ecbe5fb815cb45696f5f236c777fbb15266) | none recorded (run interactively, 2026-09-24 14:33) | matplotlib 3.10.8 (PDF metadata); this repository renders with matplotlib 3.10.7, Python 3.12.4, pandas 2.3.3, numpy 2.2.6, and the result is pixel-identical |
| composite | `P/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/prepare_supplement_figures.py` (sha256 5c52b4b89037916a4503e3b77e8454522186a83c0fdff48cced5f1776141ce4d) with `SUPPLEMENT_FIGURE_SOURCES.json` entry 6 (sha256 69943ed802e30b3182777be1f2a8347afeb7b306fdc82ffff18320a5431c5f15) | `P/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/release_candidate/supplementary_figures/Supplementary_Figure_06.pdf` (+ .png) | a0a98cee7172ae967bbf818a3239cda0352e2d6cf95002842f0a46761cf42113 (png eec3a1f385f9be8572690cdbb824a936a00b443adcb86370b46113289378f2d2) | none recorded | PyMuPDF; this repository uses PyMuPDF 1.27.1 and reproduces the PNG byte for byte |

Notebook run for this unit: Slurm (`srun`, 4 cores, 4 GB), 4.6 s wall time, 0.40 GB
peak RSS.

Upstream jobs (not re-run): panel a counts, collector job 5942965 (repaired searches
array 5942676, pilot 5942675); panels b-d, PREDICT ablation jobs 5940770 (prepare),
5940771 (pilot), 5940772 (full); CapScan ablation jobs 5928242 (pilot), 5928243 (full).

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| Supplementary Figure 8, a-d composite | Supplementary Figure 8, `word/media/image16.png`, sha256 eec3a1f3... | BYTE_MATCH; producer_status COPY_ONLY; cited script `.../v15_MM_coauthor_2026-09-25/build_figures_v15.py:36` only copies the image; "original plotting/rendering producer unresolved" | drawn by `build_revision.py` (Clinical comparator scope block), title removed and rasterised by `prepare_supplement_figures.py` | image identity agrees; producer attribution extends Codex | Codex traced the copy step only. This unit names the drawing and rasterising steps and confirms them by regeneration: the notebook output is byte-identical to image16.png. Codex lists no rebuilt candidate for this figure (REBUILT_CANDIDATES.tsv), consistent with v28f using the v28d image. |

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| `PREDICT_curation_counts_by_length.tsv` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/len950/SFig08/curation_counts_5_9.tsv` (no package sheet yet) | faf0e097d3845f94f02caa8ac46a0f01107b2abbd8d6778d27904944f51cb49f | columns sample, arm, peptides_5_50, peptides_9_50 kept; the `source` column (absolute result paths) dropped; values unchanged. Repository file sha256 01b7f7c288fdb0dc799ab00c299a7d80b70886f37e7b47fc191ad5e0c82ecc3a |
| `PREDICT_curation_counts.tsv` | `P/supplements/revision_2026-09-24/source_data/original_S13/33502dc22447_curation_counts.tsv` = `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/supplement6_repair/collection_parallel/full/curation_counts.tsv`; package sheet `Supplementary_Figure_08/PREDICT_curation_counts.tsv` | 33502dc224471a7e13c77302e65f3ac610640e19417381f0f06260fba6120661 (package sheet identical) | none (byte copy). Cross-check: the `FL_public` and `MG` counts equal `n_test` and `n_reference` of package sheet `Figure4/benchmark_0_per_acquisition.tsv` (PREDICT, min_length 5, FL_public vs MG) in all 296 acquisitions. Map rows RES-CUR1, RES-CUR2, SF8-a, SF8-cur agree (8.19005, 286, 296). |
| `clinical_fixed_dictionary_ablation.tsv` | `P/supplements/revision_2026-09-24/source_data/original_S18/e8f214c89a66_plot_values.tsv` = `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/figures/H05_H09_evidence_contribution/source_data/plot_values.tsv`, written by `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/code/plot_corrected_ablation.py` (sha256 bc46d7c3c575ece8a41b7e906564820d1abf0323c102af7c27f866cb48f5b3dd) from `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/evidence_ablation/full_PREDICT/acquisition_summary.tsv` (6c1b2d94461bf734ba22e03a0994712bddb6860004b48f1ebb1c70a4ad056b55) and `LAB_ROOT/Hartout_edits/H09_evidence_ablation/full_CAPSCAN/acquisition_summary.tsv` (fb76d1c1012d5c6ac1003380971649ef28700f38156c879bf2ffadbe5889c224) | e8f214c89a6670df364906ea34cf653932388d33e35a95342f2c521cccf70893 | none (byte copy). The package sheet `Supplementary_Figure_08/clinical_fixed_dictionary_ablation.tsv` (sha256 e28934baaa76e1571c9ff0e0fd6a5d968760259aac099b462130ad58b8346db2) is the producer's pandas re-export of this file: same rows and columns, values agree to 1.9e-16 relative (last binary digit in 29 rows). Map row SF8-bd agrees (578 rows). |
