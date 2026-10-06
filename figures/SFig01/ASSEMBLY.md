# Assembly: Supplementary Figure 1 (v28f)

The figure is a single matplotlib panel. No manual editing took place.

1. Draw. Input: `source_data/Fig1/acquisition_sizes_726.tsv` (sha256
   `fd9ade3bb4ebd2ee6c129f69819b9697cf20abd7c10c77a12644bef25396e386`). Tool:
   `figures/SFig01/panels.py:panel_a` (the SF1 block of the producer
   `LAB_ROOT/review/v15_integration/build_evidence.py`, lines 8-27, as run for the 726-row input by
   `LAB_ROOT/.../final/suppfigs_n15/SF01/build_sf1.py`). What changed relative to the v28d figure:
   input grows from 710 to 726 rows (16 participant-15 CapScan acquisitions added), CapScan tick
   label n=331 becomes n=347. Output: `panel_a.pdf`; producer output sha256
   `95c8a147dcff61366ded17d9b434e65451eb776e58d41d388b51209210bbb0be`
   (`expected/panel_a.pdf`). The repository render differs from it only in PDF metadata
   (creation date); rasterised at 300 dpi the two are pixel-identical (`checks/panel_a.json`).
2. Rasterise. Input: the PDF from step 1. Tool: `figures/SFig01/panels.py:render_composite` (the
   raster step of `LAB_ROOT/.../iteration3/source_data/figures/reexport.py`, as copied into
   `LAB_ROOT/.../final/suppfigs_n15/SF01/render.py`): PyMuPDF render at 2840 px width, which is
   400 dpi effective at the placed width of 180.34 mm (placed height 105.20 mm). Output:
   `SFig01_v28f.png`, 2840 x 1657 px.

Composite PNG `expected/SFig01_v28f.png`: sha256
`29930b6da0b412a612fc66d06d4f8afdd22a70801e0fe5ac315b9c64650963c8`. The repository render is
byte-identical to it.

Composite PDF it was rasterised from: `expected/panel_a.pdf`, sha256
`95c8a147dcff61366ded17d9b434e65451eb776e58d41d388b51209210bbb0be` (the same file as the
producer's `SupplementaryFigure1_n726.pdf`).

`LAB_ROOT/.../` abbreviates
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/`.

## v30 (house style)

1. Input `source_data/Fig1/acquisition_sizes_726.tsv` (same sha256). Tool:
   `LAB_ROOT/.../final/v30/figures_style/SFig01/build_sfig01_style.py` (repository:
   `panels.py:values_unchanged_v30`, `panel_a_v30`). Runs `panel_a` with saving intercepted, then
   redraws at the placed size 180 x 105.2 mm with `fl_style`: Söhne type, ink, selected FASTAs
   blue #0072B2, reservoirs orange #E69F00, reference line grey, value labels beside the diamonds.
   1,480 plotted values, limits and log scale asserted identical; no label string lost or added.
2. Save PDF with `fl.save` and rasterise at 300 dpi (`panels.py:render_png`). Output
   `SupplementaryFigure1.pdf`, sha256 `70152bd0628621386aeb1c2fdb44f6c89a4fdcb560da9455f93d4539ff9dbc8c`
   (= `expected/SFig01_v30.pdf`), and `SupplementaryFigure1.png`, 2126 x 1243 px, sha256
   `8978bee36972f810a7fd8416b268c768ad1b643058597867e328a980c6e7168d` (= `expected/SFig01_v30.png`).
3. Placement at 180.3 x 105.2 mm. The v30 PNG is not among the media of
   `LAB_ROOT/.../final/v30/v30f_TRACKED.docx`; placement in the next docx build is pending.

Composite v30: `expected/SFig01_v30.png` sha256 `8978bee36972f810a7fd8416b268c768ad1b643058597867e328a980c6e7168d`,
rasterised from `expected/SFig01_v30.pdf` sha256 `70152bd0628621386aeb1c2fdb44f6c89a4fdcb560da9455f93d4539ff9dbc8c`.
The notebook render is byte-identical (`checks/SFig01_v30.json`: IDENTICAL).
