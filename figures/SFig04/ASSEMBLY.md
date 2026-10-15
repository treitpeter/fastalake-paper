# Assembly: Supplementary Figure 4 (v30)

The composite is produced by one script; no manual editing step was applied.

1. Draw. Inputs: `source_data/SFig04/` (the producer inputs listed in PROVENANCE.md). Tool:
   `figures/SFig04/panels.py:_draw` (panels a-g, shared legend, `align_letters`), moved from
   `LAB_ROOT/.../final/v30/figures_style/SFig04/build_sfig04_style.py` (sha256
   `b950b4061c043e52d4c5065a85a2d548d95adfb19612b6a4cbb313753c3a2f71`), which restyled
   the v29 producer `LAB_ROOT/.../final/suppfigs_n15/SF04/build/sf4.py` (sha256
   `550e09f2be5a77b28583c210ac3f4815712a48c93785a1ba67805e1e931f379f`, SF4_MODE=new,
   n347_panel_c) with every plotted value unchanged. What changed relative to v28f/v29:
   house style (`style/fl_style`: Söhne, ink, line floor), bold panel letters at one x
   per column, subpanel titles removed, colours by meaning, legend labels "7-50" / "9-50"
   in place of ">=7" / ">=9", figure 160 x 208 mm. Written by `fl_style.save` as
   `build/SFig04/SFig04_raw.pdf`.
2. Slant taxon labels. Input: `SFig04_raw.pdf`. Tool: `panels.slant` (PyMuPDF): "L.
   plantarum", "C. butyricum", "B. longum" in d are re-set in Söhne Buch with a 0.20
   shear; the text content is asserted unchanged. Output: `build/SFig04/SFig04.pdf`.
3. Rasterise. Tool: `panels.render_png` (PyMuPDF, 300 dpi). Output:
   `build/SFig04/SFig04.png`, 1890 x 2457 px.

Composite PNG `expected/SFig04_v30.png` (the v30 build's `SupplementaryFigure4.png`):
sha256 `fb3641ee8ed2f9f5fc152696d528686b25e76344832764498d7b3fbdf03de684`. The repository
render is byte-identical to it (`checks/SFig04_v30.json`).

Composite PDF it was rasterised from: `expected/SFig04_v30.pdf` (the v30 build's
`SupplementaryFigure4.pdf`), sha256
`93e4d0a7a183e58d5256ede1e5204273d45221c1093f9f23742520fd01f36506`. The repository PDF
differs only in metadata; rasterised at 300 dpi the two are pixel-identical.

The plotted values (`expected/SFig04_plotted_values.tsv`, sha256
`e1a216dda1d258a84ae3e457e9b0d8fe67bb542c975039fae30d3ac0726a6c2c`) are the same in
v29 and v30 and are reproduced exactly.

`LAB_ROOT/.../` abbreviates
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/`.

## v35 label change

Drawn tick labels in the manuscript cohort name (`TICK_RENAME`): c "Healthy FL" -> "Healthy stool FL" (two lines), g "Healthy (58)" -> "Healthy stool (58)" (two lines). Plotted-values table unchanged (items keep the v33 labels). Expected render replaced: `expected/SFig04_v35.png` sha256 3ba4aaa20fd39addbe20cb6b55bc2d1ed52794ae7e6b68d9bc9f410b5da9963c, `expected/SFig04_v35.pdf` sha256 6b9923ef0f04bfb3d84e53ebecbc0c2adbb0c4f70f0c1315f0f36ca31bdf186a. A control build of the unchanged code reproduced the previous expected PNG byte for byte.

## v38 label change

Legend key "5–50 (published)" -> "5–50 (production)" (`LAB`). Plotted values unchanged. Expected render replaced: `expected/SFig04_v38.png` sha256 7e18a9303c2f1cc60734519294ab6e348c99cff0b69d40bce6ceea4f46b49058, `expected/SFig04_v38.pdf` sha256 e41abedf999f9b89c1297571a2d8601d9bb5376c67925abcccfec6cccf33dd17. A control build of the unchanged code reproduced the previous expected PNG (= `word/media/image11.png` of `FastaLake_v37.docx`) byte for byte; against it the v38 render differs in 14,958 pixels, all in the legend row (the centred five-entry legend shifts with the longer first key).
