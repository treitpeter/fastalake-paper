# Assembly: Supplementary Figure 11 (v30)

The v30 figure is one matplotlib canvas written by one script; no manual step and no schematic
artwork (no ARTWORK_NUMBERS.md). `FINAL/` stands for
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/`.

1. Draw. Inputs: the producer tables `functional_correlations.tsv`, `native_correlations.tsv` and
   `controlled_176_comparisons.tsv` of `FINAL/suppfigs_n15/SF11/inputs_n15/` (here:
   `source_data/SFig11/*_PREDICT.tsv` joined with the CapScan rows in `source_data/Fig5/`, and
   `source_data/SFig11/controlled_176_comparisons.tsv`; all three byte-identical to the build
   inputs). Tool: `FINAL/v30/figures/SFig11/build_sfig11.py` (sha256
   b283403ffaa26ea9a6d0f7b5f530e12304f7cfc01418097e0d80e334c08ebd82), repository
   `figures/SFig11/panels.py:render`. What changed relative to the v28f figure: native, OG and KO
   matrices removed from a and b (identical to Fig. 4i and 5g), so a and b are the KEGG pathway
   matrices; c-f unchanged in content and letter; on-canvas cohort titles, matrix column titles and
   footnotes removed (content moved to the legend); one shared colour bar labelled "Distance
   correlation"; fl_style Söhne, Okabe-Ito blue/orange, 6 pt at 180 x 190 mm placed size.
   Outputs: `SupplementaryFigure11.pdf` and `SupplementaryFigure11.png` (300 dpi, 2125 x 2244 px)
   and `plotted_values.tsv` (174 values).
2. Placement in the manuscript: the PNG of step 1, unchanged.

Composite `expected/SFig11_v30.png`: sha256
cb1e56e0c91cf3f2512a519525581f2bce4dedd4e5bb76af3720d438495b1438.
Composite PDF `expected/SFig11_v30.pdf`: sha256
d03ac4bc92f79b3e3e73662eb353d9bbc80aa83222369f9681564d6132015df2 (drawn from the same figure
object; the PNG was not rasterised from it). Plotted values `expected/plotted_values_v30.tsv`:
sha256 2597c033404714843f4266b5f731bda03d2b75adb6739b0ccf318129e3ded999.

The notebook repeats step 1 from `source_data/`; the PNG is byte-identical (same sha256), the PDF
rasterises identically at 300 dpi, and the plotted-values table is byte-identical (`checks/`).
