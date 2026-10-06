# Assembly: Supplementary Figure 13 (v28f)

The figure is one matplotlib figure with eight panels. No manual editing took place.

`SF13/` abbreviates
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/suppfigs_n15/SF13/`.

1. Select the displayed motif cells. Input: `source_data/Fig6/capscan_histidine_statistical_families.tsv`
   (sha256 `3c24a42eec725e0b43bdeecca3928f971925c25395dde7478d92a1db9e1c28d7`, A28 run on 347
   CapScan acquisitions). Tool: `figures/SFig13/panels.py:derive_display_table` (filters of
   `review_histidine_predict.py` l.16 and `plot_histidine_review.py` l.40, as copied into
   `SF13/build_sf13.py` step 0). Output: the 840-cell display table, byte-identical to
   `source_data/SFig13/capscan_histidine_motif_display.tsv`, sha256
   `269fa990cd52029a6dd4961e7bf21e347436c374c89bee4fa5c9631dd87364ea`.
2. Draw. Inputs: the display table, `source_data/SFig13/predict_histidine_acquisition_enrichment.tsv.gz`
   (patient means, `panels.py:predict_patient_means`) and
   `source_data/SFig13/predict_histidine_statistical_families.tsv` (means and intervals). Tool:
   `figures/SFig13/panels.py:draw_figure` with `panel_a`..`panel_h` (the HISTIDINE_CONTEXT_REPLICATION
   block of `build_revision.py` l.254-275, as copied into `SF13/build_sf13.py` step 1). Figure size
   7.1 x 10.8 in, typeface DejaVu Sans, with an in-image title. What changed relative to the v28d
   figure: only the CapScan input (347 instead of 331 acquisitions; participants in proximal and
   distal fluid 15 instead of 14). The colour limit moves from 13.527 to 13.427 percentage points;
   dots per panel a-f are 109, 107, 104, 107, 86, 85 (v28d 107, 110, 103, 106, 86, 85). Panels g,h
   are pixel-identical to v28d (`SF13/CHANGES_SF13.md`). Output: `HISTIDINE_CONTEXT_REPLICATION.pdf`;
   producer output sha256 `5e1983adc9aa79c778d14111a28634db8993eeb6e67e448b06f8a54bbedf1972`
   (`expected/HISTIDINE_CONTEXT_REPLICATION.pdf`); the repository render is pixel-identical at
   300 dpi and differs only in PDF metadata.
3. Redact the in-image title and rasterise. Input: the PDF from step 2. Tool:
   `figures/SFig13/panels.py:redact_and_rasterise` (the logic of `prepare_supplement_figures.py`, as
   copied into `SF13/redact_sf13.py`): the title "Histidine context and replication across cohorts"
   (box 106.537, 13.692, 404.671, 26.507 pt) is covered by a white redaction, numeric labels are
   asserted unchanged, the page is saved and rendered at 300 dpi. Outputs: `SFig13_v28f.pdf`
   (producer sha256 `25974f9e4c8758f225e94d779126110214512959a3dbfe225f260593edf9fa8d`) and
   `SFig13_v28f.png`, 2130 x 3240 px.

Composite PNG `expected/SFig13_v28f.png`: sha256
`4c26ae6bb7b1e216dbb547dc689a2817ad493fc48b807098d4c9af14e4b0dd0e`. The repository render is
byte-identical to it.

Composite PDF it was rasterised from: `expected/SFig13_v28f.pdf`, sha256
`25974f9e4c8758f225e94d779126110214512959a3dbfe225f260593edf9fa8d`.

The v28d image (`word/media/image21.png`, sha256 `d0ac4534bb5b2d02a22ff4d36fbdb122fae8e3ace08d9abc4f82c3e50730b072`)
is the same three steps on the 331-acquisition input; the producer reproduced it byte for byte
(`SF13/repro/`). For v28f the PNG above replaces it unchanged in size (7.1 x 10.8 in at 300 dpi).

## v30 composite (house style)

From manuscript v30 on, Supplementary Figure 13 is the house-style redraw of the same values. No
manual step.

1. Redraw. Inputs: the same `source_data/SFig13/` tables as above. Tool: `figures/SFig13/panels.py:figure_v30`,
   port of `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig13/build_sf13_style.py` (sha256 `3ce1554dafca059bf4ecb49bc828db4099549cee8aee66d22be4f6e900b9e99b`) and its helpers `sf_common.py` (sha256
   `9a1c99295ca64f5f7dec751e889c064a5ce00858f60058fa27bac01fea7fef3e`), with `fl_style` from
   `style/`. What changed relative to the v28f image (style only): drawn at the placed size
   160 x 243.4 mm, 1890 x 2875 px at 300 dpi; Söhne type at 6-7 pt, ink #4e4f51, lines at least 1 pt; bold panel letters a-h from `fl.panel`; the a-f panel titles moved verbatim into the y-axis labels with the v29 wording "Regions 1–2" and "Regions 3–4" (the v28f image says "Proximal" and "Distal"); the g,h titles moved into the legend entries; heat maps in `fl.div_cmap` at the same symmetric limits; g,h dots in the main Figure 6 route colours (de novo top 30% blue #0072B2, non-specific Sage green #009E73); the in-image title is not drawn. The function first runs the producer drawing
   (`panels.py:draw_figure`) with saving intercepted and compares every plotted array with the redraw exactly
   (8 panels, 1,842 plotted values; identical, as in the v30 `values_check.json`).
2. Save and rasterise. Tool: `panels.py:_v30_finish` (`sf_common.finish`): `fl.enforce_nbt`,
   `fl.normalise_type`, panel letters 7 pt, `fl.clip_check`, PDF save, then the PDF rasterised at
   300 dpi with PyMuPDF. Output: `build/SFig13/SFig13_v30.pdf` and `.png`.

Composite PNG `expected/SFig13_v30.png` (copy of `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig13/SupplementaryFigure13.png`): sha256
`890c191b2649592435d7f9e28c1d191d44cc0d4d1b19e182cf0f4884a848d41e`. The repository render is byte-identical to it (`checks/SFig13_v30.json`, IDENTICAL, 0
differing pixels; needs `FASTALAKE_FONT_DIR`).

Composite PDF it was rasterised from: `expected/SFig13_v30.pdf` (copy of
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig13/SupplementaryFigure13.pdf`), sha256 `6b650c4281252cb75b1c234c70917f49af0d9d90c775b41248f107152ad5daa6`. The repository PDF differs from it only
in metadata bytes; its 300 dpi rasterisation is the PNG above.
