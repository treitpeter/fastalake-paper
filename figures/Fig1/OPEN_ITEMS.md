# Open items: Fig1

## v30 sync (6 Oct 2026)

- AUTHOR: v30 Figure 1 is `expected/Fig1_v30.png` (sha256 b1481ff6...) from `artwork/Figure1_v30.pdf`; panel a changed (Fig1_v2 funnel), panels b to e unchanged. Rebuilt in the notebook: IDENTICAL at 300 dpi.
- AUTHOR: NUMBERS.tsv keeps the v28d columns and adds the four `v30_*` columns of tools/check_numbers.py (`v30_para` = paragraph index in visible_v30f.txt: v28d 51 -> 25, 53 -> 26, 55 -> 27, 59 and 61 -> 28, 63 -> 29, 65 -> 30, 71 -> 33) ; `check_numbers.py --manuscript-text visible_v30f.txt`: 44 rows, 0 failures, 43 MATCH, 1 NOT_IN_V30. No value changed.
- AUTHOR: v30 paragraph 25 no longer prints "the human reference contained 20,431 entries"; that row is `NOT_IN_V30` (20,431 remains in paragraph 27 and the legend). "more than 10,000-fold fewer" is now "more than 10,000-fold larger than the human reference proteome".
- UNTRACED: v30 paragraph 30 adds "In five acquisitions ... the corresponding median estimates were 1.03% and 0.51% (Source Data)" (self-selected T0 tier); not recomputed here and not in NUMBERS.tsv.
- UNTRACED: the Slurm job id of the v30 producer run of `build_fig1a_funnel.py` was not recorded; the repository rerun (job 6075277) reproduces its PNG byte-identically.

## Manuscript (v28f)

- AUTHOR (v28f, superseded by v30 above): v28d still embeds the old Figure 1 (`word/media/image1.png`, sha256 26ec1aaf..., old panel a and CapScan n = 331 in panels b and d); v28f must embed `expected/Fig1_v28f.png` (sha256 6e92d894...) from `artwork/Figure1_v28f.pdf`. The v28d legend (line 71) already describes the new figure (726 acquisitions, three-route middle column, scrambled path) and needs no change; no Results number of lines 51 to 65 changes.
- TEXT-GAP: the three donor-selected entrapment seeds (20260910, 20260911, 20260912) are not stated; Methods line 729 says only "Three seeds were averaged" (STATS.tsv rows 1 and 2; Codex S001, S002).
- TEXT-GAP: the software that computed the seed means and medians of Fig. 1e is not stated (A11 `analysis.py`, pandas mean then numpy median); Methods line 727 names only FDRBench v0.0.4 and Sage (STATS.tsv rows 1 and 2).
- TEXT-GAP: Fig. 1e shows medians only; the interquartile ranges (5-50: 0.942 to 1.036%; 9-50: 0.438 to 0.565%) are computed but not stated in text or legend (STATS.tsv rows 1 and 2).
- TEXT-GAP: the software computing the scrambled match fractions of Fig. 1c is not stated (A02 `report.py`, production Rust extractor for the lookup) (STATS.tsv rows 3 to 5; Codex S131 to S133).
- TEXT-GAP: Fig. 1b and 1d medians are over acquisitions, and CapScan and healthy stool contain several acquisitions per participant; the legend says "726 acquisitions" but does not say that points are acquisitions for b and d (it does for e).

## Check tooling (shared requests)

- AUTHOR: `tools/check_numbers.py` fails 11 NUMBERS.tsv rows whose printed values are the recomputed values rounded to two or three significant figures (250,000 for 250,197; 25,000 for 24,560; 12,600 for 12,606 twice; 18,500-fold for 18,522.7; 29,700 for 29,663 three times; 3,100 for 3,057; 12,200 for 12,212.5; 86,000 for 85,835). The rule rounds to decimals only. Request: a significant-figure rule for integers with trailing zeros, or a status for stated approximations (NUMBERS.tsv rows for lines 53, 55, 59, 61).
- AUTHOR: the printed text "more than 10,000-fold" is encoded as printed_value `>10,000` so that the bound is checked; that string is not literally inside the quote, so the `--manuscript-text` check will flag it (NUMBERS.tsv, line 51).
- PIXEL: `style/fl_style.py` is a later revision of the module the Figure 1 producers used (lab copy of fl_style.py, sha256 58b64c34..., PROVENANCE.md): 1 pt instead of 0.6 pt lines and ticks, bold instead of regular panel letters, savefig.bbox None instead of tight. With it alone, panels b, c, d differ in 2,573 to 2,997 pixels (spines, ticks, letter). `panels.py` therefore applies the producer values locally (`PRODUCER_RC`, `_style`, `_panel`) and all four panels are pixel-identical. Request: decide whether `style/` should carry the producer-era module for units whose producers used it.

## Panel a and the composite

- UNTRACED: the script of ASSEMBLY step 7 (floating-text redaction that produced the v26b fig_clean Figure1.pdf) was not saved; only its method note and REMOVED.tsv exist.
- UNTRACED: Slurm job ids and interpreter versions for `build_FAR12.py` (panel c) and `build_fig1a.py` (composite) were not recorded. The panel c render is nevertheless reproduced pixel-identically.
- AUTHOR: the composite PDF still contains the hidden text of the v10 panel c under the white rectangle painted by `build_FAR12.py`; panel-c labels and the letters c and e occur twice in the text layer. Visual output is unaffected; accessibility and text extraction are.
- AUTHOR: panel a is Illustrator artwork changed by about ten PyMuPDF edit passes; there is no single editable source (ASSEMBLY.md). Since v30 the FastaLake column is regenerated by `artwork/build_fig1a_funnel.py` on top of `artwork/Figure1_v28f.pdf`; the other two columns remain artwork.

## Data and coverage

- RESOLVED 6 Oct: the per-peptide table behind the 4,258 is now `source_data/Fig1/worked_acquisition_accepted_peptides.tsv` (also Source Data package Figure1/, map RES-17 VERIFIED); the two raw inputs (71 MB, 58 MB) stay outside the repository.
- AUTHOR: PREDICT acquisitions appear in `acquisition_sizes_726.tsv` and `acquisition_support_summary_726.tsv` by raw-file name (for example `20191024_QX7_PeTr_LC02_SA_MicrobPred_10_191025211701`). These are acquisition names, not patient identifiers, and carry no clinical column; please confirm that the run numbers do not encode patient ids.
- PIPELINE: cohort_database_sizes, database_selection_specificity, cohort_prediction_support, healthy_stool_donor_entrapment and the terminal-scramble, reference-58 and selection-grid runs are not in `pipelines/`; UPSTREAM.json lists their original code with `pipeline_dir` null.
- AUTHOR: the id `cohort_prediction_support` follows the lab folder `A09_4_1_prediction_support`; another unit uses `minlen7_search_sensitivity` for a different analysis. One of the two ids needs renumbering.
- AUTHOR: Results line 57 (reservoir sizes 233.5, 272.9, 264.9 and 375.3 million; Supplementary Fig. 1) and line 67 (bag-aware lookup, Supplementary Fig. 2e,f) sit in the Figure 1 Results section but cite supplementary figures; they are not covered here and should be covered by SFig01 and SFig02. Line 67 values may change with the pending bag reruns.
- MISMATCH: Methods line 675 ("the same 710 selected FASTAs", Supplementary Fig. 1) and line 405 ("all 710 acquisition-specific FASTAs") still use the 710 roster while Fig. 1b and line 673 use 726; outside this unit, noted for SFig01.
