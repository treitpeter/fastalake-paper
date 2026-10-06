# Assembly: Supplementary Figure 3 (v30)

The v30 figure is drawn in two automated steps from `source_data/SFig03/`; there was no manual
step. `FIN` stands for
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`.
The v28f assembly (c,d at 5-50 residues, 400-dpi image `SFig03_v28f.png`, sha256
`cb1411f832b003a1fa441800964ddf11792b6728b68b4b469f6a6c7765cdd2bd`, rasterised from the PDF with
sha256 `8b7a2b38b9c854b7debe1d309af3772f50a161e752a9119473bd44a76f57c7fd`) is superseded and its
expected files were removed.

1. Producer step, c,d at 9-50 residues (len950).
   Input: `source_data/SFig03/CAMPI_member_support_recount.tsv`, `CAMPI_reference_recovery.tsv`,
   `CAMPI_within_MSFragger.tsv` (unchanged; the 9-50 rows were already present).
   Tool: `FIN/len950/SFig03/build_sfig03_len950.py len950` with `FIN/len950/SFig03/panels_len950.py`
   (repository: `panels.py:build_figure`, `save_figure`, `slant`, with `PLOTTED_MIN_LENGTH_CD = 9`).
   A control run with `PLOTTED_MIN_LENGTH_CD = 5` reproduced the v28f PNG with 0 differing pixels.
   Change: c,d plot the nested 9-50 rows; a,b unchanged.
   Output: `FIN/len950/SFig03/SupplementaryFigure3_len950.pdf`, sha256
   `8976aefa66807bc0073ef973f959e14619f7a07de82105b33b0f2231c7336fd9` (= `expected/SFig03_len950.pdf`);
   plotted values `FIN/len950/SFig03/work/len950/plotted_values.tsv` (= `expected/plotted_values.tsv`).
2. House-style redraw (v30).
   Input: the same tables; the producer of step 1 is run first and fingerprinted.
   Tool: `FIN/v30/figures_style/SFig03/build_sfig03_style.py` (repository: `panels.py:build_v30`,
   `build_figure_v30`), matplotlib 3.10.7, `fl_style`.
   Change: Söhne type, ink, bold panel letters aligned per column; titles moved into the x-axis
   labels; arm colours as in main Figure 2; a,b value labels moved further from the bar ends;
   c widened and d narrowed (axes positions only). Every plotted array (a 32, b 32, c 170, d 68
   values), axis limit and scale, the plotted-values table and the 16 a,b value labels are
   asserted identical to step 1.
3. Slant and rasterise.
   Tool: `panels.py:slant` (16 taxon labels re-set with 0.20 shear in Söhne Buch) and
   `panels.py:render_png` at 300 dpi.
   Output: `FIN/v30/figures_style/SFig03/SupplementaryFigure3.pdf`, sha256
   `215d8d30e9bd0315d91c82112a6eea8cafd4bfb84d60f43566e1c2719ee4a1a3` (= `expected/SFig03_v30.pdf`);
   `SupplementaryFigure3.png`, 1890 x 1654 px, sha256
   `ebc18430635566789888c1ec031579ca0148aa48faf6953bdc8337d67ccaefd5` (= `expected/SFig03_v30.png`).
4. Placement: 160 x 140 mm (v29 docx extent). The v30 PNG is not among the media of
   `FIN/v30/v30f_TRACKED.docx`; placement in the next docx build is pending.

Composite: `expected/SFig03_v30.png` sha256 `ebc18430635566789888c1ec031579ca0148aa48faf6953bdc8337d67ccaefd5`,
rasterised from `expected/SFig03_v30.pdf` sha256 `215d8d30e9bd0315d91c82112a6eea8cafd4bfb84d60f43566e1c2719ee4a1a3`.

The notebook repeats steps 1-3 from `source_data/` (Slurm, paper environment, Söhne fonts) and
reproduces the step-1 PDF at 300 dpi and the v30 PNG pixel for pixel
(`checks/SFig03_len950.json`, `checks/SFig03_v30.json`: IDENTICAL, 0 differing pixels); the
plotted-values table equals the producer's.
