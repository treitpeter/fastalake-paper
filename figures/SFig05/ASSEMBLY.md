# Assembly: Supplementary Figure 5 (v30)

Supplementary Figure 5 in v30 is the v28f/v29 figure restyled in house style, with every
plotted value unchanged: two pages, panels a-f and the continuation page. Neither was
edited by hand. `FINAL` is defined in `PROVENANCE.md`.

## Panels a-f

1. Draw. Inputs: `source_data/SFig05/fusion_pooled_current.tsv`,
   `donor_entrapment_seed_means.tsv`, `simulation_mixture_outcomes.tsv`,
   `simulation_uncertainty.tsv`. Tool: `figures/SFig05/panels.py:figure_a_to_f`
   (`panel_ab`, `panel_c`, `panel_d`, `panel_ef`, `align_letters`), moved from
   `FINAL/v30/figures_style/SFig05/build_sfig05_style.py` (sha256
   `4d03dcc2b0b7508eb83009e495d41f07f187636dcf84be4cbbb307cee4a663c5`). What changed
   relative to v28f: Söhne instead of DejaVu Sans, ink colour, bold panel letters at one x
   per column, the six subpanel titles and the in-image title removed (their content moved
   to the y-axis labels; d's title to the legend text), colours by meaning (paired =
   FastaLake blue, unresolved = light grey, all dots FastaLake blue), drawn at the placed
   size 160 x 189.3 mm. Written by `style/fl_style.save` as `build/SFig05/SFig05.pdf`.
2. Rasterise. Tool: `panels._rasterise` (PyMuPDF, 300 dpi). Output:
   `build/SFig05/SFig05.png`, 1890 x 2236 px.

## Continuation page

1. Draw. Input: `source_data/SFig05/fusion_decomposition_current375.tsv`. Tool:
   `panels.figure_continued`, from the same v30 build. What changed: house style,
   blue = `fl.FL`, grey = light grey, 180 x 90 mm. Output:
   `build/SFig05/SFig05_continued.pdf`.
2. Rasterise at 300 dpi: `build/SFig05/SFig05_continued.png`, 2126 x 1063 px.

## Composite

- `expected/SFig05_v30.png` (the v30 build's `SupplementaryFigure5.png`): sha256
  `cc3d154ddd9ad10cdd041cbe88052740737ffe0774a7189328fa788e50f4bcf4`, rasterised from
  `expected/SFig05_v30.pdf`, sha256
  `ffad3ec1135bcceb414a06d5ab23f92557b5385c56c65fb39290d0e887a352d5`.
- `expected/SFig05_continued_v30.png`: sha256
  `5cfa2801bed4b38229081f391fbe447c8f255ea6eeb84a4fd107d7d24fcc89ff`, rasterised from
  `expected/SFig05_continued_v30.pdf`, sha256
  `1189ae87bbf80921c45884ed4894ed7a566dd369c33f2e8a5473edf2ce13f078`.

The repository renders are pixel-identical to both PNGs (`checks/SFig05_v30.json`,
`checks/SFig05_continued_v30.json`); the repository PDFs differ from the v30 PDFs only in
metadata and rasterise to identical pixels at 300 dpi.

## v35 label change

Drawn cohort codes replaced by manuscript names in a-c (`COHORT_LABEL`: CAPSCAN -> CapScan, HSTOOL -> Healthy stool; CAMPI, PREDICT unchanged); values and layout unchanged; continuation page unchanged. Expected render replaced: `expected/SFig05_v35.png` sha256 a39777df3a828e0cea6ec6665bb021177239114b4f38f481bc557f6be2fe61c1, `expected/SFig05_v35.pdf` sha256 5650d2bd2bc315c9460c91baa44bccaf47cfbd472bcc640f6f0e5cab7932151b. A control build of the unchanged code reproduced the previous expected PNG byte for byte.

Continuation page (v35): tick labels print the length window and seed as integers ("5–8 aa", "20260907"; were "5.0–8.0 aa", "20260907.0" from float table columns). Expected render replaced: `expected/SFig05_continued_v35.png` sha256 55f53f5565029de0ca4785b512987c4fc4ee5d189ad8651ec5f500210e632d50, `expected/SFig05_continued_v35.pdf` sha256 7fd6836444455e3e6473b96c8d1e2ef8e1031f3a294395fff4839ce829cd8a1a. A control build of the previous code reproduced the previous expected PNG byte for byte.
