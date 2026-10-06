# Assembly of the v28f Figure 6 composite

All six panels are drawn by one script on one canvas; there is no manual step. The
original chain is `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/fig6_n15/code/run_all.sh`
(below `F6N15/`); the repository runs the same steps from `panels.py`, and each output is
pixel-identical to the producer's (`checks/`).

1. Inputs for d-f. `F6N15/code/prep_histidine_inputs.py` selected the CapScan histidine
   source-enrichment and matched-reciprocal rows from the A28 fixed-roster family table.
   Repository: `panels.derive_histidine_inputs` (output byte-identical, sha256
   `aa6587cd...845b` and `6b12f31b...4109`). PREDICT rows are the unchanged end-of-day
   review tables. Script.
2. Builder. `F6N15/code/build_figure6.py` draws panels a-f on a 180 x 166 mm canvas and
   writes `panels/Figure6.pdf` and `panels/Figure6.png` (300 dpi). Repository:
   `panels.new_figure`, `panel_a` .. `panel_f`, `save_builder`. Output PDF sha256
   `5506d837465e9654abb9c82386281edac0fc5fe61544aea8909e0daae5d2b3ca`, PNG
   `f466a1e8e2886a352c97a6ef1b604441e06ddac0c2d82194ea606a1bfb4e7a25`
   (`expected/builder_Figure6.*`). Script.
3. Panel letters and trim. `F6N15/code/postprocess_fig6.py`, copied from
   `build_figures_r12.py`: the letters a-f are redacted and re-inserted in Söhne Halbfett
   8 pt; blank horizontal bands taller than 15 px at 72 dpi are removed (two bands, rows
   152-167 and 451-457). Output `post/Figure6_r12.pdf`, sha256
   `e1c961f97ca71f27f0debb5b95da852663c42947f2f000758b81eae4777af718`. Repository:
   `panels.postprocess`. Script.
4. Export width. Same script, copied from `fix_current_figures.py`: page scaled from 180 mm
   to 183 mm width (518.74 x 457.04 pt), minimum font 5.59 pt asserted >= 5 pt. Output
   `post/Figure6.pdf`, sha256 `f814a44b83bc030f486187f722fd8bd89749f97277e43f7d2f5d8fbae5bcc865`.
   Script.
5. Region names. `F6N15/code/relabel_regions.py`: 14 text spans (Duodenum, Jejunum, Ileum,
   Colon and the two panel c legend entries) covered by white boxes and rewritten as
   Region 1-4, Fluid regions 1–2 and Fluid regions 3–4 in Söhne Buch; the original text
   layer is kept and asserted unchanged. Output `region_figures/Figure6.pdf`, sha256
   `5ca34fc246bac700367017d1d5f4d8b65acf9176bf119f6d3080f115a3370576`, and its 300 dpi
   PNG. Repository: `panels.relabel_regions`. Script.
6. Final save. `run_all.sh`: PyMuPDF re-save (garbage=4, deflate) of the PDF and copy of
   the PNG to `F6N15/build_n15/Figure6.{pdf,png}`, copied to `F6N15/Figure6.{pdf,png}`.
   Repository: `panels.final_resave`. Script.

7. v30 key relabel and flatten. `CR/final/v29/figures/_tools/flatten.py --relabel` on the
   v28f PDF: the panel b key span "Host: capsule > stool in 15/15, q < 0.001" (Söhne Buch
   6.1 pt, #4e4f51, origin 313.04, 133.04 pt) is covered with its white background and
   retyped as "Human-assigned: capsule > stool in 15/15, q < 0.001" at the same origin,
   face, size and colour (changed rectangle 312.74-452.64 x 126.70-134.47 pt); the page is
   then written as one 600 dpi raster (no text, no vectors; asserted) and a 300 dpi PNG is
   rendered from the relabelled vector page. Output `CR/final/v29/figures/Fig6/Figure6.pdf`,
   sha256 `07873927e9f96c65375376cbcdfe400f05eae21c44c661c0a04a56cd2c2ae26d`, and
   `Figure6.png`, sha256 `ca90228f7902a7528a5abf340637a9bba6c47793bb3b767d498076c88f0913d1`.
   Repository: the key word is now written by `panels.panel_b` (`KEY_B_HUMAN`), and
   `panels.relabel_and_flatten` (called with no relabel) only flattens the page. Script.

Composite PNG `expected/Fig6_v28f.png`: sha256
`c822e293b35061f00124aae4bf2114d668361d29990841895f699a2c2add0b70` (2162 x 1905 px,
300 dpi), rasterised from the composite PDF `expected/Fig6_v28f.pdf`: sha256
`362a5cd71c47c25c8abafbff99b5bd9ae7a492ec4d310794c7ebda6978540392`.

v30 composite PNG `expected/Fig6_v30.png`: sha256
`ca90228f7902a7528a5abf340637a9bba6c47793bb3b767d498076c88f0913d1` (2162 x 1905 px,
300 dpi), rendered from the relabelled vector page; v30 composite PDF
`expected/Fig6_v30.pdf` (one 4323 x 3809 px raster): sha256
`07873927e9f96c65375376cbcdfe400f05eae21c44c661c0a04a56cd2c2ae26d`. This is the figure in
the v30 manuscript.

v28d embeds a different image (`word/media/image6.png`, sha256 `357bbd56...a365`, the
n = 14 build); the v28f composite above replaces it.
