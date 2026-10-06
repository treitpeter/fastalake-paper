# Assembly of Figure 3 (v28f, then v30)

No Figure 3 rebuild is queued for v28f, so the v28f image is the v28d image (`word/media/image3.png`
of `FastaLake_v28d.docx`, sha256 `6d72728a...`). It is the end of a chain of PDF edits. Each step
below names its input, tool, what changed and its output. All tools are Python scripts using
PyMuPDF unless stated. Paths are under `LAB_ROOT`. "b-e producer" and "f producer" are the two
quantitative panel builders reproduced by `panels.py`.

1. Panel renders.
   - Panel a artwork: `MM_review_analysis/figure_restoration_2026-09-16/Figure3/panels_61/a_design.pdf`
     (sha256 `1503196b...`; = `artwork/panel_a_design.pdf`), extracted from the v6.1 manuscript figure; no drawing script.
   - Panels b-e: b-e producer `.../figure_restoration_2026-09-16/Figure3/build_figure3_cdef.py` (sha256 `2ba6574e...`)
     wrote `panels_cdef/Figure3_panels_cdef.pdf` (sha256 `cf1fba8b...`; = `expected/panel_bcde.pdf`).
   - Former panel f: `panels_61/b_agreement.pdf` (sha256 `0b7446d0...`), replaced in step 7 and again in step 10.
2. `.../figure_restoration_2026-09-16/Figure3/compose_figure3.py` (sha256 `0dfbe3b2...`): placed a at the top,
   b-e at 33 mm (scaled to 180 mm width) and the former f at 130 mm on a 180 x 195 mm page; typed the letters
   a and f in Söhne. Output `Figure3.pdf` (sha256 `8529b8de...`, 510.24 x 552.76 pt).
3. Untraced: between step 2 and step 4 the page became
   `fable_5.1_review_2026-09-16/scoping_review_2026-09-19/v9_reading_2026-09-21/figures/Figure3.pdf`
   (sha256 `e9f9c858...`, same page size). The script that wrote it was not found (OPEN_ITEMS).
4. `00_PACKAGE_2026-09-17/.../audit_completion_2026-09-23/correct_figure3_label_r2.py` (sha256 `38785696...`):
   added the panel letter f in Helvetica 8 pt at (36.96, 378.0). Output
   `figure3_correction_r2/Figure3_panel_f_corrected.pdf` (sha256 `4088d54a...`).
5. `.../v11_edits2_multiomics_integration_2026-09-23/build_figures_r12.py` (sha256 `220a6b7b...`),
   functions `labels_and_small_fixes` and `trim_blank_bands`: removed the blank band 361-367 pt.
   Output `figures_r12/Figure3.pdf` (sha256 `c3120185...`, 510.24 x 546.76 pt).
6. `.../codex_review_2026-09-25/fix_current_figures.py` (sha256 `df29937a...`): scaled the 180 mm page to
   183 mm (factor 183/180) so the smallest label exceeds 5 pt. Output `release_candidate/figures/Figure3.pdf`
   (sha256 `a9997b8f...`, 518.74 x 555.87 pt), copied unchanged into the v15 co-author package by
   `v15_MM_coauthor_2026-09-25/build_figures_v15.py`.
7. `review/.../top_down_review_2026-09-30/JAG18_okabe_ito/build_figure.py` (sha256 `3cc1aea2...`): covered
   everything below 130 mm with white and placed a new f (`Figure3f_sequential.pdf`, sha256 `313b17dd...`,
   white-to-blue scale) built from the A33 protein-group and earlier KO tables. Output `JAG18_okabe_ito/Figure3.pdf`
   (sha256 `c99c9bff...`).
8. `.../JAG19_transfer/update_legend.py` (sha256 `91349c1c...`): covered the two b legend labels with white and
   typed "Same-category source (T1)" and "Inter-individual-series source (T4)" in Helvetica 5.8 pt. Output
   `JAG19_transfer/Figure3.pdf` (sha256 `5a34fb78...`).
9. `review/.../four_author_integration_2026-10-02/actual_FAR_261002/build_FAR21.py` (sha256 `bfc47ba1...`):
   replaced a by a compact strip (`Figure3a_compact_preview.pdf`, sha256 `b11b8c5f...`), moved b-e up 24 pt and
   re-typeset their 76 text spans in Söhne (letters 8 pt, axis labels 7 pt, ticks and legend 6.5 pt; each old span
   covered with white), and rebuilt f (`FAR21_Figure3_applied/Figure3f.pdf`, sha256 `3d42b7fb...`). Output
   `FAR21_Figure3_applied/Figure3.pdf` (sha256 `f5c022d8...`).
10. `review/.../claude_review_2026-10-03/iteration2/cr003_followup/scripts/build_fig3f_cr003.py verified50`
    (sha256 `65b68dac...`, Slurm job 6048652): the f producer. Drew f with the revised KO allocation
    (`shadow/verified50/fig3/Figure3f_proposed.pdf`, sha256 `de646679...`; = `expected/panel_f.pdf`) and placed it
    below 344 pt, keeping the FAR21 page above. Output `Figure3_proposed.pdf` (sha256 `48829e54...`), copied to
    `figures/corrected_verified50/Figure3.pdf`.
11. `.../final/v26b/fig3_restore/compose_fig3.py` (sha256 `7e480476...`): took the top 91 pt (panel a from step 8)
    of `JAG19_transfer/Figure3.pdf` and the step-10 page from 67 pt down, placed at 91 pt. Output
    `fig3_restore/Figure3.pdf` (sha256 `6f4d94ee...`, 518.74 x 579.87 pt).
12. `final/v26b/fig_clean`: manual, no script on disk. `REPORT.md` (sha256 `bad9e288...`) describes a PyMuPDF text
    redaction and a crop of the bottom band. Removed from Figure 3 (`REMOVED.tsv`, sha256 `6433756c...`): the panel a
    footer "Each small output symbol represents one acquired replicate · 58 acquisitions total", the f texts
    "Class labels:", "above / below" and the KO footnote. Output `fig_clean/Figure3.pdf` (sha256 `e89a7d97...`,
    518.74 x 570.88 pt).
13. `.../final/fig_v27/build_fig1_3_v27.py` (sha256 `c2c8e023...`): typed the row keys "Upper: host" /
    "Lower: microbial" and "Upper: all" / "Lower: non-human" in Söhne Buch 6.5 pt at x 458 pt, and rasterised the
    page at 300 dpi. Output `fig_v27/Figure3.pdf` and `Figure3.png`.

v28f composite PDF: `expected/Fig3_v28f.pdf`, sha256 `c5790787cc48fb1b5d97b895258436152502e2da1d6703060e999bc0441f1458`.
v28f composite PNG (300 dpi, 2162 x 2379 px; byte-identical to `word/media/image3.png` of the v28d docx):
`expected/Fig3_v28f.png`, sha256 `6d72728af6aeecaf99383770a51999a18a746e4124c5835fc26325b2cafe3d55`.

## What the band checks show

The notebook places each regenerated panel where steps 2, 6, 10 and 11 put it and compares that band of the page
with `expected/Fig3_v28f.pdf` at 300 dpi (`checks/composite_band_*.json`, `*.diff.png`). These are visual checks.

- a: 9,010 differing pixels (1.2%): the letter a (step 2), the footer removed in step 12, and anti-aliasing at icon
  edges from the two rescalings.
- b-e: 80,753 differing pixels (3.3%): all text (re-typeset in step 9, legend relabelled in step 8) and the edges of
  points and lines (sub-pixel offset from the rescalings). The plotted marks coincide.
- f: 10,629 differing pixels (0.6%): the texts removed in step 12, the keys added in step 13, and one pixel row at the
  top edge of the band. All heatmap tiles and labels coincide.

## v30 composite (9-50 residues, relabelled legend)

14. `final/fig3_labels/relabel.py` (sha256 `a7993917...`, run 6 October 2026, no Slurm job recorded): input step 13
    `fig_v27/Figure3.pdf` (sha256 `c5790787...` = `expected/Fig3_v28f.pdf`). Covered the b legend box
    (54.0, 186.8, 156.0, 204.7 pt) with white and typed "Same series" / "Different donor" in Söhne Buch 6.5 pt,
    #4e4f51. Output `fig3_labels/Figure3.pdf`, copied to `LAB_ROOT/Figure3_v28f.pdf` (sha256
    `2b9a92d6d9d769282bbfd0974b12635123ff63c3a561a400029fd5e24a96dfb8`; = `expected/Fig3_v28f_relabel_base.pdf`).
15. `final/len950/Fig3/build_fig3_len950.py 9 Figure3_len950` (sha256 `827f4359...`) with `panels_len950.py`
    (sha256 `54e4ee83...`; the repository `panels.py` of v28f with `MIN_LENGTH = 9`, nothing else changed):
    drew b-e at 9-50 residues (`Figure3_len950_panels_bcde.pdf`, sha256 `27247c89...`; = `expected/panel_bcde.pdf`),
    covered the b-e band (300-dpi rows 333-1519) of the step-14 page with white, placed the new canvas at
    `composite_rects()['bcde']`, re-typeset its 72 text spans as step 9 did (the two producer legend spans are
    skipped) and retyped the step-14 legend. Ported here as `panels.assemble_v30`. Output `Figure3_len950.pdf`
    (sha256 `a4f9c36a41187c8fd86e4e593ca2143acda83103fa7703df730fbbc5ab36bf9c`; = `expected/Fig3_v30_vector.pdf`)
    and 300-dpi `Figure3_len950.png`. The producer's receipt: 0 differing pixels outside the band against step 14;
    a control run at `MIN_LENGTH = 5` (`work/control5/`) reproduced the base. No Slurm job id recorded.
16. `final/v29/figures/_tools/flatten.py` (sha256 `d73f26fc...`): `Figure3.png` is the step-15 page rendered at
    300 dpi (0 changed pixels vs `Figure3_len950.png`); `Figure3.pdf` holds only one 600-dpi raster of that page,
    so the hidden covered text is no longer extractable (`flatten_receipt.json`). Outputs `final/v29/figures/Fig3/`,
    used unchanged in v30.

v30 composite PNG (300 dpi, 2162 x 2379 px): `expected/Fig3_v30.png`, sha256
`13589443fb69fa47dfde7a592c3ea3221548f63f4cc2ff07d4287cfc3753ec9d`. It was rasterised from the vector page
`expected/Fig3_v30_vector.pdf` (sha256 `a4f9c36a...`). The flattened v30 PDF: `expected/Fig3_v30.pdf`, sha256
`d6b5c2174750d64ff96067e76bc58c3e2bb0cc593d69d0ccbab0b44c891fa60e`.

The notebook rebuilds step 15 from the regenerated b-e canvas (`panels.assemble_v30`) and compares it with
`expected/Fig3_v30.png`: IDENTICAL, 0 differing pixels (`checks/composite_v30.json`). The band checks now compare
against `expected/Fig3_v30_vector.pdf`: a 9,037 (1.2%), b-e 65,654 (2.7%) and f 10,629 (0.6%) differing pixels,
at re-typeset or removed text and mark edges, as described above.
