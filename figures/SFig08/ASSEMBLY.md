# Assembly: Supplementary Figure 8 (v30)

The v30 composite is one matplotlib canvas. No panel was combined or edited by hand. `FIN` is
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`.

1. Recount at 9-50 residues. Tool: `FIN/len950/SFig08/build/recount_curation.py` on the four
   PREDICT arms' Sage results (rank 1, target, peptide q <= 0.01, modifications removed, I to L;
   9-50 taken from the accepted 5-50 set). Output `FIN/len950/SFig08/curation_counts_5_9.tsv`,
   sha256 faf0e097d3845f94f02caa8ac46a0f01107b2abbd8d6778d27904944f51cb49f; its 5-50 column equals
   `source_data/SFig08/PREDICT_curation_counts.tsv` in all 1,184 rows. Copied without its
   `source` column as `source_data/SFig08/PREDICT_curation_counts_by_length.tsv`.
2. Producer drawing at 9-50 (v29). Tool: `FIN/len950/SFig08/build/build.py` (sha256
   f39f49d5165cb55bef77a6527280856ddfb477779443727e70aedebcaa5e33dc), which calls the v28d producer
   functions (`panels.py:panel_a` to `panel_d`) with panel a from `peptides_9_50` (here
   `panels.py:read_tables_v30` and `v30_reference`). Output
   `FIN/len950/SFig08/SupplementaryFigure8.png`, sha256
   86c8e81db5a4c5f6b71d32d55ac2c04807696a9e31b84994b539ba162425ef6d.
3. House-style restyle. Tool: `FIN/v30/figures_style/SFig08/build_sfig08_style.py` (sha256
   9597973fcac36d3e0290fce6a84a2299c8c8bf12d352607eeeb58f1edeb1e66d; here `panels.py:render_v30`).
   Runs step 2 unchanged into a hidden figure, redraws the same values at 160 x 134 mm with
   `fl_style` (Söhne, ink, colours by meaning, panel letters, titles removed) and compares every
   plotted array exactly (`values_check.json`). Outputs `SupplementaryFigure8.pdf`, sha256
   ac4e41d6f885100e8cbfdd01ae349d872c26608a7e8b3261dd3109663b2d4ce8, and
   `SupplementaryFigure8.png` (300 dpi), sha256
   4315bca570ac1be5cc223784b8ddd60199f7533be9e59635f6acdbe26f82cd12.
4. Manuscript insertion: the PNG from step 3 at 160 x 134 mm.

`expected/SFig08_v30.png`: sha256 4315bca570ac1be5cc223784b8ddd60199f7533be9e59635f6acdbe26f82cd12.
`expected/SFig08_v30.pdf` (same canvas, same step): sha256
ac4e41d6f885100e8cbfdd01ae349d872c26608a7e8b3261dd3109663b2d4ce8.

The repository render is IDENTICAL to both at 300 dpi (`checks/SFig08_v30.json`,
`checks/SFig08_v30_pdf.json`). The v28d producer files (`CLINICAL_COMPARATOR_SCOPE.pdf/.png`)
and the title-redacted v28d composite are superseded and no longer kept in `expected/`;
`panels.py:figure` and `redact_title` remain as the v28d producer code.

## v35 label change

Panel a arm tick labels replaced by short forms of the legend names (`ARM_TICK`): FL_cohortMG -> "FastaLake + cohort MG", FL_public -> "Public FastaLake", MG -> "Specimen MG", MG_own_curated -> "Curated specimen MG"; drawn order unchanged (alphabetical by arm key). The producer/restyle tick check maps the producer labels through `TICK_RENAME`; all 40 value checks pass. Expected render replaced: `expected/SFig08_v35.png` sha256 041f4960c1bf572e6b1278de0969f3d95dbc333b769c2064e239c63fa1a72fb6, `expected/SFig08_v35.pdf` sha256 832bfaa2a6a44a60be942aa0a8e7ad4fb23f9a2bfd0bae3ea249f55025e395e4. A control build of the unchanged code reproduced the previous expected PNG byte for byte.
