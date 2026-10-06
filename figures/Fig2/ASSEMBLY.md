# Assembly of Figure 2 (v30)

The v30 Figure 2 is the len950 rebuild (panels b, d-g at 9-50 residues) with the panel b and e axis
titles changed from "identified" to "accepted" (v29). It is built by one script that plots the
canvas (step 0, `panels.build_canvas`) and repeats the six v28f PDF edit steps below (steps 1-6,
`assembly.assemble`), with two generalisations: the stream patterns of steps 2 and 3 also accept
`Tf 0 0 Td[<`, and step 2 removes however many "<10" frames and n/8 counts the canvas has.
Producer: `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v29/figures/Fig2/build_fig2_len950.py`
(sha256 e45b08327ec98300908234cfd082149e745f451e662112e8d893648c8c03076e) with `panels_len950.py`
(sha256 5dd8ccf9eadb7857d8aef504ec8a1e030bc8a795c7bd432376e21a7f02c64fa6); its BUILD_LOG.json records:

| step | tool (this repository) | what changed | producer output sha256 |
|---|---|---|---|
| 0 | `panels.py:build_canvas` | canvas, b and d-g at 9-50, "Accepted"/"accepted" axis titles | e2a75f5acfdbce0d6d201b384596b8e59f3ddd99b25c3894d54fe17db5a5a1f3 |
| 1 | `assembly.py:assemble` step 1 | letters, panel d legend keys, Söhne letters, blank band 544-553 px removed | b6e8f7818154f1a5c319053ae067449be11b87a763551761be85790d94d4ed85 |
| 2 | step 2 | panel f: 8 frame drawings ("<10" cells plus key), "Members with", "10 or more peptides", "Fewer than 10", 6/8, 7/8, 8/8 x3 removed | 7ab88d21156f80719d94900097c7a78c72474f2c6064561196f15c03b642ce7b |
| 3 | step 3 | panel g axis title, panel e keys, scaled to 183 mm | 7faef6cce9150e077fb92a2125bfe42b617b38bd2187bc0ebfb49243464e2c2d |
| 4 | step 4 | "REF" -> "Reference" overlays | 0b317b28eac327dc35c19510910bfa0edcaed311f9cfcbd683d351bf3f5852a8 |
| 5 | step 5 | 16 species labels italic | bb90e293f3b1c4d31d1acad0424d04f5177e14ecbab8dd7b88385b1ffe84726b |
| 6 | step 6 | four panel a strings redacted; 300 dpi PNG | PDF d74e1be857204c462eabba8fd96442007890c4a18c5fad66631a8259b495a336; PNG ec359e7ff04f0ae2e8434e2da64da108736bb2dc8631539b2ce0d66bbaf2965f |

Composite: `expected/Fig2_v30.png` sha256 ec359e7ff04f0ae2e8434e2da64da108736bb2dc8631539b2ce0d66bbaf2965f
(2,162 x 2,964 px), rasterised from `expected/Fig2_v30.pdf` sha256
d74e1be857204c462eabba8fd96442007890c4a18c5fad66631a8259b495a336 (copies of the producer's
`Figure2.png` and `Figure2.pdf`). The regenerated `build/Fig2/Figure2.png` has the same sha256
(pixel check IDENTICAL, `checks/composite.json`). Intermediate PDFs differ in bytes from the
producer's (PDF object metadata), not in rendering.

---

# Assembly of the v28f Figure 2 (history)

Figure 2 is unchanged between v28d and v28f (LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/fig_v27/CHANGES.md: "Figure 2, 5, 6 - No changes").
The composite was made from one plotted canvas by six PDF edit steps. No step replotted data.
Paths in steps 1 to 7 are original lab locations; `LAB_ROOT` is the lab project directory.
Abbreviation used below: `P = LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23`,
`B = LAB_ROOT/review/v15_integration/jess_integration_2026-09-29`.

| step | input | tool | what changed | output | output sha256 |
|---|---|---|---|---|---|
| 0 | `source_data/Fig2/*`, `source_data/SFig03/CAMPI_membership_*.tsv`, `artwork/bottle_61.pdf` | `figures/Fig2/panels.py:build_canvas` (moved from LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/figure2_redesign/build_campi_figure.py, sha256 dc20ed0b1b9f7edcf5d95d902907d9f3674228a24c0298864441de4939200400) | plots all panels on a 180 x 250 mm canvas with review letters a, b, c, d-f, g, h, i; places the community artwork at 3-50 x 7-40.5 mm | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/figure2_redesign/CAMPI_Figure2_review.pdf (producer; was pixel-identical to the regenerated 5-50 canvas before the v30 update) | b8c042347e03c06f1621bb9d035cb51e869daa04df602490d1b599084ccf3632 |
| 1 | step 0 output | `P/v11_edits2_multiomics_integration_2026-09-23/build_figures_r12.py` (sha256 220a6b7ba6b7a9867703a52d77189353cd10b48b50654de32a7f8d06b1a452a5), figure 2 branch, `labels_and_small_fixes`, `trim_blank_bands` | letters d-f, g, h, i redacted and replaced by d, f, g, h; letter e inserted at (260, 295) pt; all panel letters redrawn in Söhne Halbfett 8 pt; the two panel d legend keys redrawn apart (they overlapped); white bands taller than 15 px at 72 dpi removed (page 510.24 x 699.66 pt) | `P/v11_edits2_multiomics_integration_2026-09-23/figures_r12/Figure2.pdf` | 8330d46d13854401b70ecf0bae5e9a88a54644775201a3ce54f02ea98b06d40b |
| 2 | step 1 output | `P/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/FIG2F_rebuild_2026-09-25/rebuild_fig2f.py` (sha256 69e671397fd865874bbeebf1f550dc1f440604c88b09f5b27ac853c7e66a0f3d) | panel f: deleted the six orange "fewer than 10" frames (content-stream edit) and the text "Members with", "10 or more peptides", the five "n/8" counts and "Fewer than 10"; asserted the 40 cell numbers unchanged | `.../FIG2F_rebuild_2026-09-25/Figure2.pdf` | 2b0a3e4eb49fa871008c2b75b6a129be7fa8096871e05edd0705e459754aea1f |
| 3 | step 2 output | `P/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/fix_current_figures.py` (sha256 df29937ae8378539ac5c00e394b3fdf5737a16e5e57f64c485b1954bd713ce0b), figure 2 branch | panel g y-axis title "FL peptides (thousands)" deleted from the stream and reinserted at (75, 618) pt in Söhne Buch 6.3 pt; panel e three-class key (One member, Several members, No member) redrawn with spacing; numeric labels asserted unchanged; page scaled from 180 to 183 mm width | `.../codex_review_2026-09-25/release_candidate/figures/Figure2.pdf` (byte-identical copy in `.../v15_MM_coauthor_2026-09-25/coauthor_package/figures/Figure2.pdf`) | 15e0d741228e7978808d72c47613392bd833a23d83206fc00918cc13032530e1 |
| 4 | step 3 output | `B/four_author_integration_2026-10-02/actual_FAR_261002/build_FAR15_figures.py` (sha256 59a116a6fec8dd3392857e8148075f413f49fa935160b6d92da74401e3451ab5) | visible "REF" labels covered with white boxes and overwritten as "Reference" (tick labels in b and c wrapped onto two or three lines); "Unrecovered REF" became "Unrecovered reference"; two panel a text lines reworded. The original "REF" text stays in the PDF text layer under the overlays | `B/four_author_integration_2026-10-02/actual_FAR_261002/FAR15_reference_label/Figure2.pdf` | 81b7724f84b926bd164e0891baac3770003b64def02048e8819575d4e4544e4d |
| 5 | step 4 output | `B/four_author_integration_2026-10-02/actual_FAR_261002/build_FAR16_figures.py` (sha256 4a7139b0b4e520bd094d23ccb1d404766fd0a31df5fc399ac67fe5c801ec6437) | 16 species labels in panels f and g set in italic (Söhne Buch with a 0.20 shear); panel a labels were already italic | `B/four_author_integration_2026-10-02/actual_FAR_261002/FAR16_species_italics/Figure2.pdf` | 9167ec03ac927db66c1f95c58dac8dcefa272c70f3861ea0a71c23a4a9b60a1d |
| 6 | step 5 output | manual: PyMuPDF text redaction, no script saved (`B/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v26b/fig_clean/REPORT.md`, `REMOVED.tsv` rows for Figure 2) | removed four panel a text strings: "8-member gut model · 15 acquisitions" (in the artwork), "FL: selection from 375 M gut-catalogue and UniProtKB proteins;", "no community composition, CAMPI reference proteomes or specimen MG supplied.", "Reference proteomes and MG were added to the selected FL database after parsimony." | `.../final/v26b/fig_clean/Figure2.pdf` | d8a474b81d1cea981c5a3e2a1c96b28287ffe4564e722b2c2fc1662ab6466508 |
| 7 | step 6 output | file copy | none (byte-identical) | `.../final/fig_v27/Figure2.pdf` and `Figure2.png`; the PNG is `word/media/image2.png` in the v28d docx | PDF d8a474b8...; PNG deb175c906964d338090151f0760c4fcb492a9435c18b8fa17609163155a97c0 |

Each step's input hash was checked against the previous step's output: step 2 records its source
sha256 8330d46d... (REPORT.json), step 3 records 2b0a3e4e... (LATEST_FIGURE_MANIFEST.json), step 4
records the step 3 PNG c64dccea... and PDF (INTAKE.json), step 5 records 81b7724f...
(FAR16 FIGURE_VALIDATION.json), step 6 names the step 5 PDF as its source (REPORT.md).

The PNGs of steps 4 and 5 were made by pasting re-rendered text regions into the previous PNG.
The final PNG is nevertheless a plain 300 dpi render of the final PDF: rendering
`expected/Fig2_v28f.pdf` with `tools/check_pixels.py` at 300 dpi gives a result IDENTICAL to
`expected/Fig2_v28f.png` (2,162 x 2,964 px, 0 differing pixels).

Composite: `expected/Fig2_v28f.png` sha256 deb175c906964d338090151f0760c4fcb492a9435c18b8fa17609163155a97c0,
rasterised from `expected/Fig2_v28f.pdf` sha256 d8a474b81d1cea981c5a3e2a1c96b28287ffe4564e722b2c2fc1662ab6466508
(518.74 x 711.32 pt, 183 x 250.9 mm).
