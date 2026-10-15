# Assembly: Supplementary Figure 15 (v28f)

The figure is one matplotlib page holding panels a-d. It was made in two scripted steps. No
manual editing took place.

`LAB_ROOT/.../` abbreviates
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/`.

1. Draw. Input: `source_data/SFig15/compute_stage_summary.tsv` (sha256
   `37c1d94d821d5ba9defb5ff5c8495905ce93209a136952dec8b600b74816498f`, the 347-acquisition
   table). Tool: `figures/SFig15/panels.py:draw_figure`, which is the figure block of the
   producer `LAB_ROOT/00_PACKAGE_2026-09-17/.../supplements/revision_2026-09-24/build_revision.py`
   (lines 13, 34-36, 51, 59, 229-237) as copied into `LAB_ROOT/.../final/suppfigs_n15/SF15/build_sf15.py`.
   What changed relative to the v28d figure: the input table is the 347-acquisition CapScan
   version (CapScan cells of the per-acquisition stages change, see `PROVENANCE.md`), and the x
   tick labels, hard-coded in `build_revision.py` line 237 as `CapScan\n331`, are derived from
   the table, so they read `CapScan\n347`. The page carries the in-image title
   "Computational resources across manuscript cohorts". Output: `panels_abcd.pdf` (producer
   name `COMPUTE_RESOURCE_MAP.pdf`); producer output sha256
   `82177b1587cf86b1d90b70205d3b447bdf03a1faa677a69f7b624830457a61dd` (`expected/panels_abcd.pdf`).
   The repository render differs from it only in PDF metadata; rasterised at 300 dpi the two are
   pixel-identical (`checks/panels_abcd.json`).
2. Remove the title and rasterise. Input: the PDF from step 1. Tool:
   `figures/SFig15/panels.py:redact_title`, the logic of
   `LAB_ROOT/00_PACKAGE_2026-09-17/.../reappraisal_2026-09-24/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/prepare_supplement_figures.py`
   as copied into `LAB_ROOT/.../final/suppfigs_n15/SF15/redact_sf15.py` (the title box is found
   by its text instead of a stored box; the box is the same, 93.67, 8.08, 417.51, 20.89 pt). What
   changed: the title is covered by a white PyMuPDF redaction; the multiset of numeric strings on
   the page is asserted unchanged; the PDF is saved (`garbage=4, deflate=True`) and rendered at
   300 dpi. Outputs: `SFig15_v28f.pdf` (511.2 x 496.8 pt) and `SFig15_v28f.png` (2130 x 2070
   px, 338 dpi effective at the placed size of 160.00 x 155.49 mm).

Composite PNG `expected/SFig15_v28f.png`: sha256
`b343e653aa38c78d90f6753f25225eec8ad4077f5ff2b2b6d85d31dfc286fc56` (producer
`LAB_ROOT/.../final/suppfigs_n15/SF15/new/Supplementary_Figure_15.png`). The repository render
is byte-identical to it (`checks/SFig15_v28f_png.json`).

Composite PDF it was rasterised from: `expected/SFig15_v28f.pdf`, sha256
`02d9c5ee79cf5e33a87912ec3e6b6612ffbfe01a65dbd67b523b5a91c6a11fdb` (producer
`LAB_ROOT/.../final/suppfigs_n15/SF15/new/Supplementary_Figure_15.pdf`). The repository render is
pixel-identical at 300 dpi (`checks/SFig15_v28f.json`) and differs only in PDF metadata.

## v30 composite (house style)

From manuscript v30 on, Supplementary Figure 15 is the house-style redraw of the same values. No
manual step.

1. Redraw. Inputs: the same `source_data/SFig15/` tables as above. Tool: `figures/SFig15/panels.py:figure_v30`,
   port of `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig15/build_sf15_style.py` (sha256 `3d7476c4e4cccfb328c5d2cc0d1a0433c0348c962e59204d37aa670046e0ffb4`) and its helpers `sf_common.py` (sha256
   `9a1c99295ca64f5f7dec751e889c064a5ce00858f60058fa27bac01fea7fef3e`), with `fl_style` from
   `style/`. What changed relative to the v28f image (style only): drawn at the placed size
   160 x 155.5 mm, 1890 x 1837 px at 300 dpi; Söhne type at 7 pt, ink #4e4f51, heat-map frame removed; bold panel letters a-d; panel titles removed (their content is the legend's a-d sentences); shading white to Okabe-Ito blue #0072B2 on the same log10(1 + value) with the same limits; cell text white or ink by higher contrast on the new cell colour; the in-image title is not drawn. The function first runs the producer drawing
   (`panels.py:draw_figure`) with saving intercepted and compares every plotted array with the redraw exactly
   (4 panels, 112 heat-map values and 112 printed values; identical, as in the v30 `values_check.json`).
2. Save and rasterise. Tool: `panels.py:_v30_finish` (`sf_common.finish`): `fl.enforce_nbt`,
   `fl.normalise_type`, panel letters 7 pt, `fl.clip_check`, PDF save, then the PDF rasterised at
   300 dpi with PyMuPDF. Output: `build/SFig15/SFig15_v30.pdf` and `.png`.

Composite PNG `expected/SFig15_v30.png` (copy of `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig15/SupplementaryFigure15.png`): sha256
`fcf6a785b93e5d024844c1a5aec8bc2b142b5390d46e1dd1576f8a258000e469`. The repository render is byte-identical to it (`checks/SFig15_v30.json`, IDENTICAL, 0
differing pixels; needs `FASTALAKE_FONT_DIR`).

Composite PDF it was rasterised from: `expected/SFig15_v30.pdf` (copy of
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig15/SupplementaryFigure15.pdf`), sha256 `aec64c675e39d1561ab19e13196ad65e35520c9d5945786fa293f013fdf7661b`. The repository PDF differs from it only
in metadata bytes; its 300 dpi rasterisation is the PNG above.

## v35 label change

Column label "Stool" -> "Healthy stool" in a-d, and cell values that `.2g` prints in scientific notation drawn as plain numbers (`_v35_plain`; b, CAMPI shared lookup "1.4e+02" -> "140"). v30 drawing only; the producer tick strings and cell texts are mapped to the new forms for the exact comparison. Expected render replaced: `expected/SFig15_v35.png` sha256 58aefdcd5eee4746dd9c671157ff3ac32a30d8caf89715f103c606ebeeccd7ec, `expected/SFig15_v35.pdf` sha256 1c34d2a4fa2c89210473721bef079b2179c96788c45281c7fd6101383b9d7ce5. Control builds of the previous code reproduced each previous expected PNG byte for byte.
