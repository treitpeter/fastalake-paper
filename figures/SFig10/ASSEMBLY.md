# Assembly: Supplementary Figure 10 (v30)

The v30 figure is a single matplotlib canvas written by one script. No panel was edited,
cropped or placed by hand.

1. Render. Inputs: `source_data/SFig10/` (patients.tsv, focused_values.tsv,
   focused_features.tsv, producer_plotted_ab.tsv) plus the controlled
   `patient_order_and_covariates.tsv` (antibiotic track). Tool: LAB_ROOT/.../final/v30/figures/SFig10/build_sfig10.py
   (sha256 86bb5e3018da4553707483c34b97e77c35b7b0394c2b67095e27c56e1a0ab197), here
   `figures/SFig10/panels.py:render`. It calls the earlier producer's drawing functions
   (`tracks`, `matrix`, `colorbar`) on a 180 x 62 mm canvas in fl_style at 6 pt, and saves PDF
   and a 300 dpi PNG with `fl.save`. What changed relative to the v28f figure: panel b (the
   19-feature overview, identical to Fig. 4h) removed; the remaining panel has no letter and now
   carries the patient-position ticks and x label; Söhne via fl_style, lines at least 1 pt;
   values unchanged (the 510 heatmap cells equal the producer's panel-a plotted values).
   Outputs: SupplementaryFigure10.png (2125 x 732 px), SupplementaryFigure10.pdf.
2. Placement in the manuscript. Input: the PNG of step 1, inserted unchanged.

Composite PNG `SFig10_v30.png`: sha256
9b85a672a7c39d87371ec10e5d5a8c5e7361590e150a33303d8c13bebb5a6fd8. Composite PDF
`SFig10_v30.pdf`: sha256 7093d6897d946f253380df3a0f7431dc92bde1e85c3237783a1da427433d909b
(drawn from the same figure object; the PNG was not rasterised from it). The repository render
is pixel-identical to both (`checks/SFig10_v30.json`, `checks/SFig10_v30_pdf.json`). Both
files show antibiotic state per patient and are kept in the controlled directory
(`FASTALAKE_CONTROLLED_DIR`), not in `expected/`, pending the author ruling in OPEN_ITEMS.md.

`LAB_ROOT/.../final/` abbreviates
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/`.
