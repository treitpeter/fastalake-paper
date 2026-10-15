# Assembly: Supplementary Figure 14

The composite is the producer's single figure with its in-image title removed. There is no
manual step.

1. Render. Inputs: `source_data/SFig14/PREDICT_spectrum_case_review.tsv` and
   `source_data/SFig14/spectra/HIS_00{1..4}/{processed_peaks,matched_fragments}.tsv` (the
   producer read the same bytes from its own copies; the case table there had CRLF line
   ends, same values). Tool:
   `LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/supplements/revision_2026-09-24/build_revision.py`
   lines 278-291 and its `save()` (lines 34-52), ported as `panels.py:figure_abcd`. What
   changed: nothing; this is the first render. It adds the title "Measured PREDICT
   histidine-boundary spectra" and moves overlapping ion labels up with gray leader lines.
   Output: `HISTIDINE_MEASURED_SPECTRA.{pdf,svg,png}` (png 180 dpi, 1278 x 1584 px), pdf
   sha256 `4b027ebdadb153ca2ce3d297d7fdfd5843e4335ded4db21132032f5de31e7a50`, png sha256
   `cf3013f9404667db0d069770d69d658edb7e6d2dcd212ae5533a4e27362507c3`.
2. Title removal and rasterisation. Tool:
   `LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/prepare_supplement_figures.py`
   lines 7-17, with the title span located by `inspect_package_inputs.py` line 16 (same
   folder); ported as `panels.py:remove_title`. What changed: the title text span (bbox
   110.858, 10.812, 400.334, 23.627 pt) was covered by a white redaction; numeric labels were
   checked unchanged; the PDF was rasterised at 300 dpi (2130 x 2640 px). Output:
   `release_candidate/supplementary_figures/Supplementary_Figure_10.{pdf,png}` (numbered 10
   at that time), pdf sha256
   `1239e2c7b56f37bdf0765d5d2aecdb7212ab8db15b93ea5df1af6420a2879ddd`, png sha256
   `4b3e84a7d0791c0841e7866a490c1f8488b17e5ac123f9e1cfe0e0ce90051f37`.
3. Renumbering copy. Tool:
   `LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/v15_MM_coauthor_2026-09-25/build_figures_v15.py`
   line 36 (`shutil.copy2`, old number 10 to new number 12). What changed: file name only.
   Output: `coauthor_package/supplementary_figures/Supplementary_Figure_12.png`, same sha256.
4. Embed. The PNG bytes were carried unchanged through later manuscript versions (the figure
   became Supplementary Figure 14 by renumbering). In `FastaLake_v28d.docx` the image is
   `word/media/image22.png`, sha256
   `4b3e84a7d0791c0841e7866a490c1f8488b17e5ac123f9e1cfe0e0ce90051f37` (checked 2026-10-06).
   The script that wrote this media part into the v28d package was not identified (see
   `OPEN_ITEMS.md`); the bytes are identical, so the step did not change the image.
5. v28f. No change; Supplementary Figure 14 is not in the v28f rebuild queue.

Composite: `expected/SFig14_v28f.png` sha256
`4b3e84a7d0791c0841e7866a490c1f8488b17e5ac123f9e1cfe0e0ce90051f37`, rasterised at 300 dpi
from `expected/SFig14_v28f.pdf` sha256
`1239e2c7b56f37bdf0765d5d2aecdb7212ab8db15b93ea5df1af6420a2879ddd`. The notebook's
regenerated composite PNG has the same sha256; its PDF differs in bytes but rasterises
pixel-identically (`checks/SFig14_v28f_pdf.json`).

## v30 composite (house style)

From manuscript v30 on, Supplementary Figure 14 is the house-style redraw of the same values. No
manual step.

1. Redraw. Inputs: the same `source_data/SFig14/` tables as above. Tool: `figures/SFig14/panels.py:figure_v30`,
   port of `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig14/build_sf14_style.py` (sha256 `66cd5621524cb3ef1a4597aafff0ab1e71f5ee6d9a53be239b3cb15f0dc78c38`) and its helpers `sf_common.py` (sha256
   `9a1c99295ca64f5f7dec751e889c064a5ce00858f60058fa27bac01fea7fef3e`), with `fl_style` from
   `style/`. What changed relative to the v28f image (style only): drawn at the placed size
   160 x 198.3 mm, 1890 x 2343 px at 300 dpi; Söhne type at 5.04-7.2 pt (charge-2 ion labels as mathtext at 7.2 pt), ink #4e4f51, lines 1 pt; bold panel letters a-d; the panel titles (case ID, peptide, origin, boundary) moved into the x-axis labels; b ions blue #0072B2, y ions orange #E69F00, processed peaks grey; the producer label-displacement rule re-applied on the new layout (1 label displaced, with a leader line); the in-image title is not drawn. The function first runs the producer drawing
   (`panels.py:figure_abcd`) with saving intercepted and compares every plotted array with the redraw exactly
   (4 panels, 916 plotted values and 30 ion labels; identical, as in the v30 `values_check.json`).
2. Save and rasterise. Tool: `panels.py:_v30_finish` (`sf_common.finish`): `fl.enforce_nbt`,
   `fl.normalise_type`, panel letters 7 pt, `fl.clip_check`, PDF save, then the PDF rasterised at
   300 dpi with PyMuPDF. Output: `build/SFig14/SFig14_v30.pdf` and `.png`.

Composite PNG `expected/SFig14_v30.png` (copy of `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig14/SupplementaryFigure14.png`): sha256
`0ce9a0c46aa000236b3dbe5fa60e93a7fcc6024915e05d7f4320d5b0c7f14b14`. The repository render is byte-identical to it (`checks/SFig14_v30.json`, IDENTICAL, 0
differing pixels; needs `FASTALAKE_FONT_DIR`).

Composite PDF it was rasterised from: `expected/SFig14_v30.pdf` (copy of
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig14/SupplementaryFigure14.pdf`), sha256 `f5aa8c121daf384713fd195e1546f682f3a0aff27001304a98ced6493feb25e2`. The repository PDF differs from it only
in metadata bytes; its 300 dpi rasterisation is the PNG above.

## v35 label change

x-axis label case IDs HIS_001-HIS_004 drawn as "Case 1"-"Case 4" (panel order a-d; Source Data folders keep the HIS_00k names). v30 drawing only; plotted values unchanged. Expected render replaced: `expected/SFig14_v35.png` sha256 cbc4afbd25d8615c711c145e5c385206524b1f2a1a805e52c2f14d74b2d7c384, `expected/SFig14_v35.pdf` sha256 dc8a9c563bb4675a5e98c4289f3693938a001deced4f43f100f4bb92c299c8ab. A control build of the previous code reproduced the previous expected PNG byte for byte.
