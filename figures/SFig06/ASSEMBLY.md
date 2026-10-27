# Assembly: Supplementary Figure 6 (v30)

The v30 figure is one matplotlib canvas. No manual editing took place.

1. Draw. Inputs: `source_data/SFig06/split_half_acquisitions_726.tsv` (sha256
   `21a0dea9ce7045e3835ada95627e4d693d843fdd7fa308438391c7cff0cf1007`, rows with
   min_length 9), `source_data/SFig06/same_spectrum_benchmark_0_predict.tsv` (sha256
   `bf59b19aa057f36621559267b74baf70d8923ac8739dce996ae722c40da56455`) and
   `same_spectrum_benchmark_1_capscan.tsv` (sha256
   `37b1cd553132ab1338efdf3dd8b93f5514078fe8d057dd4064e3666c8cc5355c`). Tool:
   `figures/SFig06/panels.py:figure` (`panel_a`, `panel_bc`, `panel_de`, `align_letters`),
   moved from `LAB_ROOT/.../final/v30/figures_style/SFig06/build_sfig06_style.py` (sha256
   `6653782043d04c61810d2e650e583d958414f4aa64940314cc4c53c4dc4b7488`), which restyled the
   v29 producer `LAB_ROOT/.../final/v29/figures/SFig06/build/panels_len950.py`
   (`_render_abc`, `_render_de`; minimum length 9, roster 726) with every plotted value
   unchanged. What changed relative to v28f: 9-50 residues instead of 5-50, 726
   acquisitions instead of 710 (CapScan 347, two disjoint zeros), labels All-prediction /
   Same-half / Disjoint database, one canvas at 180 x 186 mm instead of two producer
   canvases composed with PyMuPDF, house style (Söhne, bold letters, no subpanel titles,
   colours by meaning). The PDF is written by `style/fl_style.save`.
   Output: `build/SFig06/SFig06.pdf`.
2. Rasterise. Tool: `panels.figure` (PyMuPDF, 300 dpi). Output: `build/SFig06/SFig06.png`,
   2126 x 2197 px.

Composite PNG `expected/SFig06_v30.png` (the v30 build's `SupplementaryFigure6.png`):
sha256 `39b4c460d06b303564740a5c94f72abdc017e3337271ce12d07c8bea30d415e1`. The repository
render is pixel-identical to it (`checks/SFig06_v30.json`, 0 differing pixels).

Composite PDF it was rasterised from: `expected/SFig06_v30.pdf` (the v30 build's
`SupplementaryFigure6.pdf`), sha256
`ff9407bf6d687a351db16c0108dc38552f8b09626303c5dfc73b5dce56e7450a`. The repository PDF
differs from it only in PDF metadata.

`LAB_ROOT/.../` abbreviates
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/`.

## v36 label change

Panel d labels print the median paired gain 100 x (FL/MG - 1) per cohort (PREDICT +9.4%, CapScan +40.1%), the statistic the text reports, instead of the per-arm medians (`panels.py:panel_de`). Plotted-values table unchanged. Expected render replaced by the rebuild embedded in v36 (`LAB_ROOT/.../final/figs_v34_fig4b/sf6/SFig06.{png,pdf}`): `expected/SFig06_v36.png` sha256 e9c32339e5f0d46138d496205934487acc53fe7e41e75478f60a63e222b1aae4, `expected/SFig06_v36.pdf` sha256 3fae003aebaa46532671d7a809327b03d45b43c9cf389e7dea1fa24e736df659.

## v38 label change

Panels b and c tick labels "<n> runs" -> "n = <n>" (`LABELS_BC`; "<n> acquisitions" was tried and overlaps neighbouring labels at the placed size). Panel d keeps the v36 median paired change labels. Plotted values unchanged. Expected render replaced: `expected/SFig06_v38.png` sha256 467f7ca18fda57d5e580c619fa18888594c3641fb7bc2b5661c20f9e174096d9, `expected/SFig06_v38.pdf` sha256 1acefc7cd4c9c86aa1d5a74d2f534b9244036d9fd34603b217f54c972ed0aa84. A control build of the unchanged code reproduced the previous expected PNG (= `word/media/image14.png` of `FastaLake_v37.docx`) byte for byte; against it the v38 render differs in 9,296 pixels, all in the second tick-label line of b and c.
