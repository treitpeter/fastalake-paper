# Assembly: SFig07 (v30)

The v30 composite is one matplotlib canvas. No panel files were combined and no manual
editing took place. `FIN` below is
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`.

1. Recount at 9-50 residues. Inputs: the accepted 5-50-residue identity tables of the
   healthy-stool searches and the per-acquisition clinical Sage results (outside this
   repository). Tools: `FIN/len950/SFig07/work/hstool_recount.py` (srun job 6072769) and
   `FIN/len950/SFig07/work/clin_recount.py` (Slurm array 6072755), aggregated by
   `FIN/len950/SFig07/build_sfig07_len950.py`. What changed: panel a and the peptide axis of
   panel b now count peptides of 9-50 residues. Outputs: the tables now in
   `source_data/SFig07/identity_per_acquisition.tsv`,
   `common_search_participants.tsv`, `completeness_by_length.tsv`,
   `identity_summary_by_length.tsv` (hashes in PROVENANCE.md).
2. Producer drawing at 9-50. Tool: `FIN/len950/SFig07/panels_len950.py` (the v28d producer with
   `MIN_LENGTH_PANEL_A = 9`; here `panels.py:render`). Output: the v29 figure
   `FIN/len950/SFig07/Supplementary_Figure7_len950.png`, sha256
   7bb2e76578f216a49c49cb2d44bd0430952a65e3175dd9701371c57a9e7af8c2.
3. House-style restyle. Tool: `FIN/v30/figures_style/SFig07/build_sfig07_style.py` (here
   `panels.py:render_v30`). Runs step 2 unchanged into a hidden figure, redraws the same values
   with `fl_style` at 160 x 144 mm (Söhne, ink #4e4f51, colours by meaning, panel letters a-d,
   subpanel titles moved to x labels or the legend) and compares every plotted array exactly
   (54 checks, `values_check.json`). Outputs: `SupplementaryFigure7.pdf`, sha256
   c364e69ef304c0695f21faa2b11bd414d76fe3c4e997830701d2c5eb1494e4f4, and
   `SupplementaryFigure7.png` (300 dpi), sha256
   ab999c1d525baddfd84e2357c5271faa12180e9fea793c00f4994ca7739c19e3.
4. Manuscript insertion: the PNG from step 3, placed at 160.0 x 144.2 mm.

`expected/SFig07_v30.png`: sha256 ab999c1d525baddfd84e2357c5271faa12180e9fea793c00f4994ca7739c19e3.
`expected/SFig07_v30.pdf` (written from the same canvas in the same step; the PNG is not
rasterised from it): sha256 c364e69ef304c0695f21faa2b11bd414d76fe3c4e997830701d2c5eb1494e4f4.

The repository render (`build/SFig07/SupplementaryFigure7.png` and `.pdf`) is IDENTICAL to
both at 300 dpi (`checks/SFig07_v30.json`, `checks/SFig07_v30_pdf.json`). The v28d/v28f
composite (`Supplementary_Figure7_revised.png`, sha256 b88113c0...) is superseded and no longer
kept in `expected/`.
