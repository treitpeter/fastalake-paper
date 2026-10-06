# Supplementary Figure 10

PREDICT metabolite-enzyme pairs for 85 patients ordered by the metagenomic Enterotype
Dysbiosis Score (EDS). The figure (v30) is a single panel: three metabolite-enzyme (KO) pairs
selected post hoc from an exploratory screen of the same cohort, with EDS, human-signal and
antibiotic tracks. It is a descriptive display; the statistics behind it (Spearman
correlations, covariate adjustments and the candidate rule) are recomputed here from the saved
screen tables and the six displayed measurements. The former panel b (the 19-feature overview,
identical to Fig. 4h) was removed in v30; its tables stay in `source_data/SFig10/` because the
Methods, the legend and Supplementary Table 4 still print its counts (20 features, 259
observations).

In this unit the column `v28d_line` of NUMBERS.tsv and STATS.tsv holds the paragraph index of
the v30 text (`visible_v30f.txt`, first column), not a v28d line. Every quote was checked
against that text with `tools/check_numbers.py --manuscript-text` (one paragraph per line).

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a (single panel, printed without a letter) | six-row heatmap (three metabolite-KO pairs) under EDS, human-share and antibiotic tracks | inverse metabolite-enzyme associations, ρ=−0.53, −0.44, −0.35; depth adjustment weakens them (v30 paragraphs 73, 145, 239, 409, 410) | 85 patients (84 for tryptophan), one observation each | Spearman rank correlation; partial rank correlations; Benjamini-Hochberg in the complete screen family | Panel a; Statistics for panel a | patients.tsv, focused_values.tsv, focused_features.tsv, augmented_pair_screen.tsv, all_EDS_associations.tsv, producer_plotted_ab.tsv; controlled patient_order_and_covariates.tsv | predict_eds_screen, predict_eds_focused_pairs |

The notebook renders the figure (`build/SFig10/panel_a.png`, `.pdf`) with the v30 build
layout and pixel-checks it against the v30 figure (`SFig10_v30.png`, `SFig10_v30.pdf`, sha256
in ASSEMBLY.md). Both checks are IDENTICAL (0 differing pixels) with the Söhne fonts and the
controlled input (receipts `checks/SFig10_v30.json`, `checks/SFig10_v30_pdf.json`).

The figure shows antibiotic state per patient. Under decision 2 of the release rules, no
image of it is kept in this repository: the v30 images are read from
`FASTALAKE_CONTROLLED_DIR` like the covariate table, the render and the diff images stay in
the gitignored `build/SFig10/`, and the notebook does not display the figure. Whether the
published figure itself may be shipped here awaits an author ruling (OPEN_ITEMS.md).

The clip at ±2 is read from the producer's plotted values (`producer_plotted_ab.tsv`, rows of
both former panels): the largest absolute drawn value is 2.0 while z scores reach 7.66 before
clipping, and 87 of 1,921 drawn cells (4.5%) are clipped. The 510 cells drawn in v30 equal
that table's panel-a rows exactly (`panels.render`).

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| N_PATIENTS | 85 | earliest matched proteome-metabolome observation of each identified patient (v30 paragraph 409) | 259 observations with repeated visits: Supplementary Table 4, legacy_259_values.tsv |
| CLIP_Z | 2 | display range, values standardized within feature and clipped at ±2 (v30 paragraphs 145, 239) | none |
| MIN_POSITIVE_N | 43 | association tests need positive finite values in at least 43 patients (v30 paragraph 409) | none |
| Q_CUTOFF | 0.05 | candidate pairs pass q<0.05 in each complete unadjusted family (v30 paragraph 409) | none |
| PAIRS | glucuronic acid-K01812, tryptophan-K01667, glycochenodeoxycholic acid-K01442 | chosen post hoc among 21 candidates for interpretable reaction links and distinct topics (v30 paragraph 409) | all 21 candidates: augmented_pair_screen.tsv, Supplementary Table 4 |
| COVARIATES | human share, protein-group depth, antibiotic state, clinical group | adjustment set of the screen (v30 paragraph 409) | all 16 combinations: separate_adjustments_all_models.tsv |
| NORMALIZATIONS | native, annotated_relative | KO signal as measured, or divided by total KO signal (v30 paragraph 409) | both in separate_adjustments_all_models.tsv |
| SCOPES | pair_complete, common_84 | pair-specific complete subsets and the subset shared by all three pairs (v30 paragraph 409) | both in separate_adjustments_all_models.tsv |
| JOINT_MODEL | host_depth_antibiotic_clinical | joint adjustment named in v30 paragraphs 73 and 409 | EDS-adjusted model in augmented_pair_screen.tsv |
| SEPARATORS_A | [1.5, 3.5] | producer layout: lines between the three pairs | none |
| FIG_MM, FONT_PT, PNG_DPI | 180 x 62 mm, 6 pt, 300 dpi | v30 build: authored at the placed size, fl_style | none |
| HEATMAP_LR, TRACK_TOP_HEIGHT, HEATMAP_RECT_Y, COLORBAR_RECT | (0.30, 0.905), (0.975, 0.055), (0.19, 0.60), [0.925, 0.25, 0.010, 0.48] | v30 build layout (figure fractions) | none |
| CONTROLLED_FILE | patient_order_and_covariates.tsv | antibiotic state and clinical group are controlled (CONVENTIONS section 8) | none |

## Known limitations

- The three pairs were selected post hoc after inspection of the complete screen; the
  same cohort was used for selection and display. No held-out patients.
- Depth-only and single-covariate adjustments are descriptive point estimates without P
  values or intervals; the code computes no confidence interval for any correlation in
  this figure.
- The original 20-feature selection screen (predict_eds_feature_selection) is not traced to
  code in this unit; its counts are taken from the fixed feature list.
- Without the controlled input the figure is not rendered (the antibiotic track is part
  of it) and the joint and single-covariate estimates are skipped; all printed numbers in
  NUMBERS.tsv are still recomputed.
- The pixel check needs the v30 images from the controlled directory; without them it
  prints SKIPPED and the committed receipts are left as they are.
