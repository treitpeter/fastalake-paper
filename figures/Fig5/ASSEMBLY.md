# Assembly: Figure 5

The v30 Figure 5 differs from v28f only in panel f. In v30 the original protein grouping was
removed from the manuscript; panel f shows direct peptide assignment only (all eligible functions
46.4%, strict subset 57.6%, one dot per participant, n = 15). Panels a-e and g are unchanged.

v30 steps (all scripted, Slurm; the notebook repeats steps 2-4 and reproduces every output pixel for pixel):

1. Panel a. Unchanged from v28f (step 1 of the v28f assembly below). `artwork/panel_a.pdf`,
   sha256 5072d68a9cd1fb4ef6e3579de7706b58f4192a9b416f9ce15b1d8b707f10bc6d.
2. Panels b-g. Input: the source tables in `PROVENANCE.md`, with `participant_leaders_resolved.tsv`
   and `participant_leaders_strict_both_methods.tsv` replaced by the v30 direct-assignment recount.
   Tool: `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures/Fig5_strict/build.py` calling the v30 producer copy of `panels.py`
   (`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures/Fig5_strict/producer/figures/Fig5/panels.py`, sha256
   b40a123faca7f26cc04af27bd919c7c318dafd9de5b812985c0fc880a7de7067; only `panel_f` differs from
   the v28f code), launched by `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures/Fig5_strict/job.sh`. Here: `panels.py` functions
   `panel_b` to `panel_g` and `save_body`. Output: `expected/panels_bg_v30.pdf`, sha256
   bf72dd1007d04ea21ee5ed2fc53f635eadc7cd4647015cadb5d866de65f99328.
3. Composition. Tool: `panels.compose` (PyMuPDF `show_pdf_page`). Output: `expected/Fig5_v30.pdf`,
   sha256 362a81c55c6a6973be90e00d2f557ebd458550fb37dabfb3e88a82e4049c0f36 (the rebuilt PDF differs in
   bytes because of PDF metadata and is identical at 300 dpi).
4. Rasterisation. Tool: `panels.compose` (PyMuPDF `get_pixmap(dpi=300)`). Output:
   `expected/Fig5_v30.png`, 2162 x 2868 px; the notebook's rebuilt PNG is byte-identical.
5. Per-panel references. `expected/panel_b.png` to `expected/panel_g.png` are crops of
   `expected/panels_bg_v30.pdf` made with `panels.crop_panels` (boxes `panels.PANEL_BOXES_PT`),
   the same crop the notebook applies to its own render. `expected/panel_a.png` is unchanged.

Composite `expected/Fig5_v30.png` sha256: b2ab579d475b74dbbd83d73e9cdc0839bd9ab20eb4f8ca93cd236172216e85e0.
Composite PDF it was rasterised from, `expected/Fig5_v30.pdf` sha256: 362a81c55c6a6973be90e00d2f557ebd458550fb37dabfb3e88a82e4049c0f36.

## v28f assembly (superseded for panel f; the v28f expected files were removed from `expected/`)

The v28f Figure 5 is the v28d Figure 5 unchanged. Steps 1 to 4 were scripted; the tool of
step 5 is not identified. The notebook repeats steps 2 to 4 and reproduces every output pixel for pixel.

1. Panel a. Input: icon vector paths from `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/small/fig5c/out/fig5_verified/Figure5_v21_refit_bcd.pdf`
   (sha256 6340bdf3581bcaa3f1c1222e2392c457a5de701e74fb0ebc3bf2bed53e573023) and the Söhne fonts.
   Tool: `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v27/n15/build_fig5a_n15.py` (PyMuPDF drawing; counts typed as constants).
   Change: draws the 518.74 x 140 pt study-design band with the 347-acquisition counts.
   Output: `artwork/panel_a.pdf`, sha256 5072d68a9cd1fb4ef6e3579de7706b58f4192a9b416f9ce15b1d8b707f10bc6d.
2. Panels b-g. Input: the twelve source tables listed in `PROVENANCE.md`. Tool: `panels.py`
   functions `panel_b` to `panel_g` and `save_body` (originally `build_figure5_n15.py` lines 56-292).
   Change: draws all six panels on one 518.74 x 548.28 pt page. Output: `panels_bg_n15.pdf`,
   sha256 2dcd79a1e939190ae791b19a3d40551a0b14ceb526f7014d5e454ca5449e1b60 (the rebuilt PDF differs
   in bytes because of PDF metadata, and is identical at 300 dpi).
3. Composition. Input: outputs of steps 1 and 2. Tool: `panels.compose` (originally
   `build_figure5_n15.py` lines 294-299, PyMuPDF `show_pdf_page`). Change: places panel a in the
   band y = 0-140 pt and panels b-g in y = 140-688.28 pt of a 183 x 242.8 mm page. Output:
   `Fig5_v28f.pdf`, sha256 f180eda7935f2939d5eefa27f2d94394a42b7af2dee60a525995465eb997e9b8.
4. Rasterisation. Input: the composite page. Tool: `panels.compose` (originally lines 300-302;
   PyMuPDF `get_pixmap(dpi=300)`, then Pillow sets the 300-dpi tag). Output: `Fig5_v28f.png`,
   2162 x 2868 px.
5. Insertion. The PNG is `word/media/image5.png` of `FastaLake_v28d.docx`, byte-identical
   (sha256 below). The tool that swapped it into the v27 working docx, and from there into v28d, is
   not identified (`final/v27/n15/REPORT.md` states only that it was "swapped into v27"); see
   `OPEN_ITEMS.md`. The image bytes are unchanged by insertion.

v28f composite `Fig5_v28f.png` sha256: 57c575ada96e16beabc0503f03b641719a7e65ba9f2b454af0c9357e21f48a69.
v28f composite PDF `Fig5_v28f.pdf` sha256: f180eda7935f2939d5eefa27f2d94394a42b7af2dee60a525995465eb997e9b8.
