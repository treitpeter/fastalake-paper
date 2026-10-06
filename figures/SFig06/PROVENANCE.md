# Provenance: Supplementary Figure 6 (v30)

`SPLIT` below abbreviates `LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/split_half_repair`
and `FINAL` abbreviates `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`.

## Panels

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a-e (v30, one canvas) | FINAL/v30/figures_style/SFig06/build_sfig06_style.py (sha256 66537820...), helper FINAL/v30/figures_style/_sf456_common.py (sha256 c2a03a2c...), run by FINAL/v30/figures_style/run_sf456.sh | FINAL/v30/figures_style/SFig06/SupplementaryFigure6.pdf and .png | pdf ff9407bf6d687a351db16c0108dc38552f8b09626303c5dfc73b5dce56e7450a; png 39b4c460d06b303564740a5c94f72abdc017e3337271ce12d07c8bea30d415e1 | not recorded in the build folder (OPEN_ITEMS.md) | CPython 3.12.4, matplotlib 3.10.7, NumPy 2.2.6, pandas 2.3.3, PyMuPDF 1.27.1 (the repository environment, `fastalake-paper_v28_work/env`) |
| a-e (v29, values) | FINAL/v29/figures/SFig06/build/panels_len950.py (sha256 f478668c...; `_render_abc`, `_render_de`, LEN950_MIN = 9, roster 726), from FINAL/len950/SFig06/build/panels_len950.py | FINAL/v29/figures/SFig06/build/render_9_726/panel_de_plotted_values.tsv | equal to `expected/panel_de_plotted_values.tsv` (sha256 ac8fcc41...) | not recorded | as above |

The v30 build compared the fingerprint of every axes (scatter offsets, IQR and median
lines, labels, limits, ticks) with the v29 producer run unchanged and found them
identical, except the panel-a y range (schematic, no data) and the removed header and
footnote texts (FINAL/v30/figures_style/SFig06/CHANGES.md, build_report.json).

Pixel results (receipt in `checks/`): the repository PNG is IDENTICAL to
`expected/SFig06_v30.png` (0 differing pixels); the repository PDF rasterised at 300 dpi
is IDENTICAL to `expected/SFig06_v30.pdf`. The d,e plotted values equal the v30 table.

## Legend edits

The v30 legend changes (b,c footnote moved into the legend; "Numbers beside the median
bars in d") are drafted in FINAL/v30/figures_style/SFig06/LEGEND_CHANGES.md. The v30 text
(`FINAL/v30/visible_v30f.txt`, line 231) does not yet carry them (OPEN_ITEMS.md).

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| v28d embedded image (old a-d) | Supplementary Figure 6, word/media/image14.png, sha256 b4ec6cf5... | BYTE_MATCH to FINAL/figure_sd/figure/out/SupplementaryFigure6_final_400dpi.png | replaced in v30 by FINAL/v30/figures_style/SFig06/SupplementaryFigure6.png | not applicable | Codex census predates v30 |
| v30 figure | no census row | none | build_sfig06_style.py, regenerated here pixel-identically | not applicable | |

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| split_half_acquisitions_726.tsv | FINAL/public347/spectrum_split347/collect/ACQUISITION_REVIEW.tsv (collect347.py, sha256 146be0c8...) | e0f257046686fb7e551bf2455539ab2f917eec88e5d4deec9339401bb0340e79 | figures/SFig06/prepare_source_data.py: drop `participant`, `specimen` |
| same_spectrum_benchmark_0_predict.tsv | Source Data package sheet Figure4/benchmark_0_per_acquisition.tsv (origin legacy_v12/Figure4/5f5b68c2947a_0_per_acquisition.tsv) | fdd5383ae384b7ec59202820a3e501ab2d80a32fd8e98318d5747d96ca667808 | same script: drop `specimen_id` and `group` (PREDICT visit code) |
| same_spectrum_benchmark_1_capscan.tsv | Source Data package sheet Figure4/benchmark_1_per_acquisition.tsv (origin legacy_v12/Figure4/93b3a849f202_1_per_acquisition.tsv) | 65ac87033c11d58c20246a6a70ac747109ce901c10dc952974262c8051ff8e31 | same script: drop `specimen_id` |
| split_half_structure.tsv | the 726 table's `participant` and `specimen` columns, min_length 9 | e0f257046686fb7e551bf2455539ab2f917eec88e5d4deec9339401bb0340e79 | same script, function `structure`: counts per cohort only (rerun on Slurm, job 6075284) |

The v30 build read the same three data tables from this repository's `source_data/SFig06/`
(the rerun of `prepare_source_data.py` reproduced them byte for byte). The cohort medians
at 9-50 residues on 726 acquisitions equal FINAL/len950/SFig06/NEW_VALUES.tsv
(`build/stats_9_726.json`). The d,e values at 9-50 residues are not yet in the Source Data
package (OPEN_ITEMS.md).
