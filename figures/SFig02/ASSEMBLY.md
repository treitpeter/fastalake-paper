# Assembly: Supplementary Figure 2

Supplementary Figure 2 is two images in the manuscript: panels a-d (v28d `word/media/image8.png`,
placed 160.0 x 135.2 mm) and the continuation with panels e,f (v28d `word/media/image9.png`, placed
154.94 x 117.07 mm). No panel was edited by hand. v28f uses both images unchanged; no rebuilt
Supplementary Figure 2 is queued for v28f. Paths shortened with `...` are given in full in
`PROVENANCE.md`.

## Panels a-d

1. Input `source_data/SFig02/selection_length_lookup.tsv` and the workbook sheets now in
   `selection_grid_summary.tsv`. Tool: `LAB_ROOT/00_PACKAGE_2026-09-17/.../supplements/revision_2026-09-24/build_revision.py`,
   selection block (lines 61-87) and its `save()` (lines 34-36, 51). Draws the 2 x 2 canvas,
   adds the in-image title "Evidence selection and exact-lookup specificity" and saves the vector
   PDF. Output `SELECTION_SPECIFICITY.pdf`, sha256 `f3b4d4a1c620eebbbc979ce2ec46fce2ff2728aa6db79accaa0a45c06496c7bc`
   (= `expected/panel_abcd.pdf`). Repository equivalent: `panels.py:panel_abcd`.
2. Input the PDF of step 1. Tool: `LAB_ROOT/.../codex_review_2026-09-25/prepare_supplement_figures.py`
   (script). Removes the in-image title by a white redaction of its text box (PDF points
   108.65, 6.78, 402.57, 19.60), checks that every numeric label is unchanged, saves the PDF and
   renders it at 300 dpi. Outputs `Supplementary_Figure_01.pdf`, sha256
   `b4148cb2513c2972c650bf2eb3eea59948144aa261727b4107920536212a3d6e` (= `expected/SFig02_v28f.pdf`),
   and `Supplementary_Figure_01.png` (2130 x 1800 px), sha256
   `9dacf0331290f5719b51bc50aff870730fa8318c797d236040d3bb1195be1154` (= `expected/SFig02_v28f.png`).
   Repository equivalent: `panels.py:compose_abcd`.
3. Input the PNG of step 2. Inserted into the manuscript unchanged (v28d `image8.png` has the same
   sha256). The insertion tool is not recorded; the bytes did not change.

## Panels e,f (continuation)

1. Input the 48 completion receipts (now `source_data/SFig02/bags_bounded48_per_acquisition.tsv`)
   and the fragment-audit table (now `bags_fragment_audit_summary.tsv`). Tool:
   `LAB_ROOT/review/.../bags_prominence_review_2026-09-30/integrate_v18.py`, which executes the
   plotting block of `audit_and_plot.py` (lines 46-51) after replacing the panel letters a,b by
   e,f. Output `Supplementary_Figure2ef_bags.pdf`, sha256
   `1f84736efe2a2efa85d270246dd72763ba8a33a5cb8a760e011b2cc5d4a88a00` (= `expected/panel_ef.pdf`).
   Repository equivalent: `panels.py:panel_ef`.
2. Input the PDF of step 1. Tool: `LAB_ROOT/review/.../iteration3/source_data/figures/reexport.py`
   (script, Slurm job 6048949). Renders the PDF with MuPDF at 2440 px width, 400 dpi at the placed
   width of 6.1 in. Output `SupplementaryFigure2_continued_v21_reexport.png` (2440 x 1844 px), sha256
   `f4231f5f090a014bf19e764d156067577c597df499d6efdbb92f31b3705ce806`
   (= `expected/SFig02_continued_v28f.png`). Repository equivalent: `panels.py:compose_ef`.
3. Input the PNG of step 2. Inserted into the manuscript unchanged (v28d `image9.png` has the same
   sha256).

## Composites

| file | sha256 | rasterised from | sha256 of source |
|---|---|---|---|
| `expected/SFig02_v28f.png` (a-d) | `9dacf0331290f5719b51bc50aff870730fa8318c797d236040d3bb1195be1154` | `expected/SFig02_v28f.pdf` | `b4148cb2513c2972c650bf2eb3eea59948144aa261727b4107920536212a3d6e` |
| `expected/SFig02_continued_v28f.png` (e,f) | `f4231f5f090a014bf19e764d156067577c597df499d6efdbb92f31b3705ce806` | `expected/panel_ef.pdf` | `1f84736efe2a2efa85d270246dd72763ba8a33a5cb8a760e011b2cc5d4a88a00` |

The notebook regenerates both composites; the regenerated PNGs have the same sha256 as above.

## v30 (house style)

`FIN` = `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`.

1. Input `source_data/SFig02/` (unchanged). Tool: `FIN/v30/figures_style/SFig02/build_sfig02_style.py`
   (repository: `panels.py:build_v30`, `panel_abcd_v30`, `panel_ef_v30`). Runs `panel_abcd` and
   `panel_ef` with saving intercepted, then redraws both canvases with `fl_style` at the placed
   widths (a-d 160 x 135.2 mm, e,f 155 x 117.1 mm): Söhne type, ink, Okabe-Ito colours by meaning,
   subpanel titles removed, net values of e as right-hand tick labels, ≤/≥ labels composed from
   Söhne glyphs, one c/d key in d, inset raised. Every plotted array (a 504, inset 420, b 84, c 60,
   d 100, e 52, f 8 values), limit and scale is asserted identical to the producers.
2. Save PDF with `fl.save` and rasterise at 300 dpi (`panels.py:render_png`).
   Outputs `FIN/v30/figures_style/SFig02/SupplementaryFigure2.pdf` (sha256
   `d0c213812c86b27018cc958fba357045e6cd8699b2563acdf24feb292908c51d`, = `expected/SFig02_v30.pdf`),
   `SupplementaryFigure2.png` (sha256 `7ac52bafc4f63cab447d6ebb1053e5ecbd79db80065fd972473adc996d463c1a`,
   = `expected/SFig02_v30.png`), `SupplementaryFigure2_continued.pdf` (sha256
   `9e574d2382333a6821e9e6b7d6e762e9bbc8bf4bd86db70ade5ee2e16a9e6df1`, = `expected/SFig02_continued_v30.pdf`)
   and `SupplementaryFigure2_continued.png` (sha256
   `d079324364dcfbd7177a30638e7cea9defc5beb65baa5cc000904a61fb71671e`, = `expected/SFig02_continued_v30.png`).
3. Placement at 160.0 x 135.2 mm and 154.9 x 117.1 mm. The v30 PNGs are not among the media of
   `FIN/v30/v30f_TRACKED.docx`; placement in the next docx build is pending.

The notebook repeats steps 1-2 from `source_data/` and reproduces both v30 PNGs pixel for pixel
(`checks/SFig02_v30.json`, `checks/SFig02_continued_v30.json`: IDENTICAL).
