# Figure 1: assembly of the v28f composite

The v28f Figure 1 is `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/fig1a_perfecter/Figure1.pdf`
(round 8 of `build_fig1a.py`). It is a chain of PDF edits: each step opened the previous
step's PDF and redacted, painted over or overlaid parts of it. Below, `FINAL` stands for
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`
and `FAR` for `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/four_author_integration_2026-10-02/actual_FAR_261002`.

Steps 1 to 5 are taken from `FINAL/fig1a_perfecter/SOURCE.md` and the receipts it names;
they were not re-run here. Steps 6 to 9 were checked against their receipts and output
hashes.

1. Root artwork. Input: Illustrator file `LAB_ROOT/FIGURATOR/Figure1/FastaLake_Figure1_v4_JBMR.ai`
   (Adobe Illustrator 29.8, 7 Sep 2026). Tool: manual, Adobe Illustrator (J. Mueller-Reif).
   Three columns, reservoir, query boxes, five-band funnel, icons, spectra, output boxes.
2. Label substitutions. Tool: `LAB_ROOT/MM_review_analysis/figure_restoration_2026-09-16/Figure1/panel_a_jbmr_v10/build_panel_a_jbmr_v10.py`
   (v10 text edits applied in memory), composed with the v10 bottom row by
   `LAB_ROOT/MM_review_analysis/figure_restoration_2026-09-16/Figure1/v10_compose/compose_figure1_v10.py`
   (Slurm jobs 5966808, 5966809, 18 Sep 2026). The bottom row is
   `LAB_ROOT/MM_review_analysis/figure_restoration_2026-09-16/Figure1/panels_bcd_v10/Figure1_panels_bcd_v10.pdf`
   (sha256 `1dc4f23dde58a06724099b56609909b77e0753daa498763853119d322636c064`), which is
   the source of panel e in every later step. Copied here as `expected/panels_bcde_row_v10.pdf`.
3. Study-level stage, table and typography (20 to 21 Sep 2026). Tools: `edit_figure1a.py`,
   `refresh_figure1_typography.py`, `build_figures.py` in
   `LAB_ROOT/fable_5.1_review_2026-09-16/scoping_review_2026-09-19/`. Inserted the
   study-level band, the four-row Real/Scrambled table, retyped panel-a text in Söhne.
4. Lookup box and arrows (23 Sep 2026). Tool: `build_figures_r12.py`, function
   `edit_figure1()`, in `LAB_ROOT/00_PACKAGE_2026-09-17/.../v11_edits2_multiomics_integration_2026-09-23/`.
5. Width 180 to 183 mm (25 Sep 2026). Tool: `fix_current_figures.py` in
   `LAB_ROOT/00_PACKAGE_2026-09-17/.../codex_review_2026-09-25/`, scale 1.016667. Output:
   v15 co-author package `Figure1.pdf`, sha256 `2cdddb3cb87b4eece22d98b9f623f93a0bc079bac22e4dd35469f9ce09955873`.
6. Panel c replaced (2 Oct 2026). Tool: `FAR/build_FAR12.py` (sha256 `757b492fe3224100eb4319dea63e6032384ea767079e52834cfe66513dd1e8af`).
   Rendered `FAR/FAR12_figure1c_applied/Figure1c.pdf` (sha256
   `ecff68781c968ced8e2ef9439f2226383fb7d114a44e774bfc900a1349289671`, = `expected/panel_c.pdf`),
   painted a white rectangle over (132, 374, 264, 517) pt of the step-5 page and overlaid
   the new panel there with `show_pdf_page`. The old panel-c text was not removed, only
   covered (see `OPEN_ITEMS.md`). Receipt `FAR/FAR12_figure1c_applied/VALIDATION.json`:
   all other pixels unchanged at 300 dpi. Output `Figure1.pdf` sha256
   `e95846727c28bdc40ac82f6898b2bdb66397ebe089cc82b72541ebcf2560b260`.
7. Floating-text removal (5 Oct 2026). Tool: PyMuPDF text-only redaction; the script was
   not saved (method in `FINAL/v26b/fig_clean/REPORT.md`, removed strings in
   `FINAL/v26b/fig_clean/REMOVED.tsv`). Removed five floating labels from panel a; line
   art untouched. Output `FINAL/v26b/fig_clean/Figure1.pdf`, sha256
   `302cef6ccd2584857152e68c492173c99ea2c6a1a045e7f49bae7aa85cde3fe2`.
8. Panel letter c (5 Oct 2026). Tool: `FINAL/fig_v27/build_fig1_3_v27.py` (sha256
   `c2c8e023671c6aa6bacf6d464954a0f7c3b2501c75b132c337024d3b612747c8`). Replaced the
   regular 8.0 pt "c" by Sohne-Halbfett 8.13 pt to match b, d, e. Output
   `FINAL/fig_v27/Figure1.pdf`, sha256 `432699dea4b0f12142de3ef76c96e9adffa8bb26b3ac7b7debaf336a66e6140d`;
   its 300 dpi PNG (sha256 `26ec1aaf189e5a09f364766e5bf28243dd13e3188f71206875bfdce50920b7ba`)
   is the image embedded in v28d.
9. FastaLake column redrawn and panels b, d at n = 347 (5 Oct 2026, rounds 2 to 8). Tool:
   `FINAL/fig1a_perfecter/build_fig1a.py` (sha256
   `26714201117a52ccf6ddc0c0aeec77a45b7df110dbaab61a02c18f305a13ccfa`), run with PyMuPDF
   1.27. Changes are listed in `FINAL/fig1a_perfecter/CHANGES.md` and the receipt
   `FINAL/fig1a_perfecter/_build_receipt.json`:
   - panel a, FastaLake column only (bbox 125.5, 6.0, 513.6, 361.9 pt): step order
     MS/MS, de novo, top-30% queries, study-level lookup, funnel; scrambled path through the
     same reservoir; one arrow helper for all arrows; floating text removed; "Any search
     engine" replaced by "Database search"; middle column redrawn as three lanes; title
     "AlphaNovo peptide predictions" (round 7). Comparison table kept as source objects.
   - panels b and d: the old vectors in (-1.57, 375.37, 127.42, 526.48) and
     (265.31, 375.37, 390.55, 526.48) pt were removed and the producer crops
     `FINAL/fig1bd_n15/panel_b.pdf` (sha256 `1fcc8730a6f64fb6320ad1027f10e08256bfdd6de11a2f22ea0d2e302a3e916c`)
     and `panel_d.pdf` (sha256 `88cfb3382ef98b9d8fcee01d7617c7b143e0a8fc7109819d3591a307753562fe`)
     were placed at scale 1.02244 in (-1.472, 375.573, 127.317, 526.282) and
     (265.413, 375.573, 390.454, 526.282) pt. These are `expected/panel_b.pdf` and
     `expected/panel_d.pdf`.
   - receipt assertions: 0 changed pixels at 300 dpi in panels c and e and in columns 1
     and 2.

   Output `FINAL/fig1a_perfecter/Figure1.pdf`, 518.74 x 525.55 pt (183 mm wide), and
   `Figure1.png`, 2162 x 2190 px at 300 dpi.

10. v29 funnel (6 Oct 2026). Tool: `FINAL/v29/figures/Fig1/build_fig1a_funnel.py`; FastaLake
   column redrawn as one funnel on the step-9 page; band 3 "Per-acquisition match"; one
   orange #d94d00. Output `FINAL/v29/figures/Fig1/Figure1.pdf`, sha256
   `cb28f0d27cfcae47f432e7357aa0ad5f9b5b2fee6daed1774f60efff8d57a968` (kept as
   `artwork/Figure1_v29_intermediate.pdf`).
11. v30 funnel, Fig1_v2 (6 Oct 2026). Tool: `FINAL/v29/figures/Fig1_v2/build_fig1a_funnel.py`
   (sha256 `de4b1216...`), ported unchanged except for paths as `artwork/build_fig1a_funnel.py`.
   Base: the step-9 PDF (`artwork/Figure1_v28f.pdf`); art clips from `FINAL/fig_v27/Figure1.pdf`
   (`artwork/Figure1_v27_art_source.pdf`). Changes, FastaLake column only: AlphaNovo and
   scrambled boxes side by side above the rim with one arrow AlphaNovo to scrambled; MS/MS
   and De novo centred on the column axis; one arrow from each box into the rim; band 2
   "Exact match · study-supported proteins"; table cell Scrambled / Study-supported "–" to
   "5,841"; scrambled result-box spectrum recoloured #e85d0f to #d94d00. Build assertions:
   0 changed pixels at 300 dpi outside the column against step 9 and step 10. Output
   `FINAL/v29/figures/Fig1_v2/Figure1.pdf` (sha256 `0f686c53...`) and `Figure1.png`
   (sha256 `b1481ff6...`). Panels b to e are unchanged from step 9.

## In this repository

- `artwork/Figure1_v28f.pdf` is the step-9 PDF, byte-identical (sha256
  `8bc5a81b489f1ab5fb7ec23f695c0cc60e92d26346faff4255973acb952c5128`).
- `artwork/panel_a.pdf` is the same file with its CropBox set to (0, 0, 518.74, 372) pt
  by PyMuPDF `set_cropbox` (2026-10-06); no object was changed or redrawn. sha256
  `92a13d812d9cf35b2b048c88cb3071b88f4e5b8ec3487c6fb47e9389bde4c23a`.
- `expected/Fig1_v28f.png` is the step-9 PNG, byte-identical. The notebook rasterises
  `artwork/Figure1_v28f.pdf` at 300 dpi with PyMuPDF 1.27.1 and finds it pixel-identical
  to this PNG (`checks/Fig1_v28f.json`).

- `artwork/Figure1_v30.pdf` and `expected/Fig1_v30.png` are the step-11 PDF and PNG,
  byte-identical. `artwork/panel_a.pdf` is now the step-11 PDF with CropBox (0, 0, 518.74,
  372) pt (sha256 `26dd41a2...`). The notebook reruns `artwork/build_fig1a_funnel.py`
  (Slurm job 6075277) and the result is pixel-identical to `expected/Fig1_v30.png` at
  300 dpi (`checks/Fig1_v30_rebuild.json`); `artwork/Figure1_v30.pdf` is too
  (`checks/Fig1_v30.json`).

v30 composite: `expected/Fig1_v30.png` sha256
`b1481ff689fb4439824ec5ec4f567b03e2d14860fa2445e6dfbabb94cd4fcb95`, rasterised from
`artwork/Figure1_v30.pdf` sha256
`0f686c53a29d7c360ace01aceae7c2918d8aed450f40311c4016c35e4e6cc23f`.

v28f composite (base of step 11): `expected/Fig1_v28f.png` sha256
`6e92d894478bc9e532bf9fa42e51bc84e76fdc7e4ed6df70d1dda84a40eb3f4a`, rasterised from
`artwork/Figure1_v28f.pdf` sha256
`8bc5a81b489f1ab5fb7ec23f695c0cc60e92d26346faff4255973acb952c5128`.
