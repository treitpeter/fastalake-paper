# Release check (v30 manuscript state)

Date: 2026-10-06. Checks on the exact staged file set (`git add -A`, then `git ls-files`):

| # | check | result | detail |
|---|---|---|---|
| H1 | `tools/scan_release.py` on the staged files | PASS | 0 findings, 2,075 files |
| H2 | no controlled or clinical data | PASS | no `controlled/` file except `CONTROLLED_INPUTS.md`; no per-patient antibiotic state, raw patient id or clinical covariate (age, sex, MELD, CRP, BMI, visit date, clinical group) column in any staged table; patient keys are pseudonyms (`PREDICT_Pnnn`); no image of the antibiotic-state panels (Fig. 4e,g,h; Supplementary Fig. 10) |
| H3 | file sizes | PASS | 447 MB in 2,075 files; largest `source_data/Fig6/capscan_residue_participant_region.tsv` 19.7 MB; none above 100 MB |

## Changes made for this check

1. `tools/check_numbers.py`: the trailing-zero rounding rule now also lists v30 values as ROUNDED (it already applied to `v28f_change` values; Fig. 2 17,100, 18,000 and 18,200 pass as ROUNDED). New `v30_status` value `REMOVED_V30` (sentence cut deliberately between v28d and v30; `v30_quote` holds the reason; listed, never fails).
2. Figure 6: the panel b key "Human-assigned" is written by the plot code (`panels.panel_b`, `KEY_B_HUMAN`); the cover-and-retype on the composite is gone and the v30 step only flattens. Pixel check against the v30 manuscript figure: 5,884 differing pixels (PNG, 300 dpi), all inside the panel b key (x 1266-1882, y 535-589 px): glyph rasterisation and a ~1 px lower "Microbial" key row. Panels a, c-f IDENTICAL.
3. Figure 4: the controlled package copy of panel h now carries "Enterotype dysbiosis score (0–1)". The published Fig. 4 composites (`Fig4_v28f.*`, `Fig4_v30.*`), the vector panel-a artwork (which carries the whole page) and the v30-vs-v28f diff image show panels e, g, h; they were moved to the controlled package and added to `.gitignore`. Panel a ships as `figures/Fig4/artwork/panel_a.png` (crop of the v30 composite). The notebook reads the composites from `FASTALAKE_CONTROLLED_DIR` when present. The committed Fig4 and Fig6 notebooks were executed without the controlled directory.
4. The 149 number rows not located in v30: 95 relocated and MATCH, 25 relocated and OPEN (v28d row already OPEN, or v30 prints a value the code does not return), 29 REMOVED_V30 with reason, 0 genuinely missing. `check_numbers.py` against `visible_v30f.txt` afterwards: 33 files, 1,354 rows, 154 OPEN, 0 failures, 0 NOT_IN_V30, 38 ROUNDED, 29 REMOVED_V30.

## Gate on a fresh copy

GATE_RESULTS_PENDING
