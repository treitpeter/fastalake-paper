# Supplementary Figure 14

The figure shows four measured MS2 spectra from PREDICT stool, one for each combination of
origin (human, microbial) and boundary side (N, C), for peptides identified by the
non-specific Sage search with a histidine at a protein boundary. It lets a reader see the
fragment evidence behind the histidine-boundary peptides discussed with Fig. 6f (v28d line
177). The figure is in v28d unchanged and is not rebuilt for v28f; it is PREDICT only and does
not depend on the CapScan roster.

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | HIS_001, human CTRC peptide SLFVGVDTIHVH, C-terminal boundary histidine | "Four selected PREDICT spectra" (465); "Supplementary Figs. 13,14 show sequence-context and selected-spectrum checks" (177) | 1 PSM, 1 acquisition | none | Panel a | `source_data/SFig14/spectra/HIS_001/`, `PREDICT_spectrum_case_review.tsv` | histidine_spectrum_evidence |
| b | HIS_002, human IPSP peptide SNIQVSEM[+15.9949]VHK, N-terminal boundary histidine inferred from the source protein | same | 1 PSM, 1 acquisition | none | Panel b | `source_data/SFig14/spectra/HIS_002/` | histidine_spectrum_evidence |
| c | HIS_003, microbial peptide VIEDAGYGPYFGH, C-terminal boundary histidine | same | 1 PSM, 1 acquisition | none | Panel c | `source_data/SFig14/spectra/HIS_003/` | histidine_spectrum_evidence |
| d | HIS_004, microbial peptide FHPALAPIK, N-terminal boundary histidine inferred from the source protein | same | 1 PSM, 1 acquisition | none | Panel d | `source_data/SFig14/spectra/HIS_004/` | histidine_spectrum_evidence |

The producer drew all four panels in one matplotlib figure, so the notebook renders them
together (`build/SFig14/panel_abcd.{pdf,svg,png}`), then removes the figure title and
rasterises at 300 dpi as the release step did (`build/SFig14/SFig14_v28f.{pdf,png}`). All four
pixel checks are IDENTICAL, and the regenerated composite PNG is byte-identical to the image in
v28d. The notebook also re-derives each case's processed peaks and matched b/y ions from the
measured peaks, and checks the case-selection rule against the selection tables.

## v30 composite

From manuscript v30 on, the figure is the house-style redraw of the same values (`panels.py:figure_v30`, notebook section `v30 composite (house style)`; ASSEMBLY.md, second part). Panels and lettering are unchanged; only type, colours, size and the placement of panel titles changed. The function compares every plotted array with the producer drawing and raises on any difference. The regenerated PNG is byte-identical to the v30 image (`checks/SFig14_v30.json`).

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `PRODUCER_RC` | DejaVu Sans, 8 pt, savefig 180 dpi | producer's matplotlib settings | none |
| `COL` | #2878a3 (b ions), #df8f44 (y ions) | producer palette | none |
| `PEAK_COLOR`, `PEAK_LW`, `ION_LW` | #b5b5b5, 0.6, 0.9 | producer drawing | none |
| `LABEL_MIN_RELATIVE` | 0.08 | ions above 8% of the maximum processed peak are labelled | none |
| `YLIM` | 0 to 137% | producer axis | none |
| `FIGSIZE` | 7.1 x 8.8 in | producer figure | none |
| `LABEL_STEP`, `LABEL_ATTEMPTS`, `LABEL_CEILING` | 10, 12, 120 | producer's displacement of overlapping ion labels (legend: "thin gray leader lines") | none |
| `COMPOSITE_DPI` | 300 | release rasterisation of the title-removed PDF | none |
| `TITLE_SPAN_MAX_Y`, `TITLE_SPAN_ASSERT_Y` | 35 pt, 20 pt | how the release step found the title span | none |
| `FRAGMENT_MZ` | 200 to 2000 m/z | Sage fragment window of the non-specific search (A25 search_revision2.py lines 72, 80); not stated in Methods for the plotted spectrum (TEXT-GAP in OPEN_ITEMS.md) | none |
| `MAX_PEAKS` | 150 | Sage max_peaks; processed spectrum = 150 most intense peaks in the window (A25 line 71); Methods state 150 only for healthy stool, v28d line 625 (TEXT-GAP) | none |
| `FRAGMENT_PPM` | 20 | Sage fragment tolerance; Methods, v28d lines 625 and 855 | none |
| `FRAGMENT_CHARGE_MAX` | 2 | Sage max_fragment_charge (A25 line 71); a precursor of charge z uses fragment charges 1 to min(z - 1, 2); Methods state it only for healthy stool, v28d line 625 (TEXT-GAP) | none |
| `PPM_ATOL` | 1e-6 | producer's ppm consistency check | none |
| `PROTON_F64` | 1.0072764 | proton mass in the producer's ppm check; `PROTON` is its float32 form used by Sage and the A25 code | none |
| `SAME_SCAN_Q` | 0.01 | peptide and spectrum q gate for `accepted_in_semi` (case audit line 158); Methods, v28d line 855 | none |
| `V30_SIZE_MM`, `V30_H_PAD`, `V30_W_PAD` | 160 x 198.3 mm; constrained-layout pads | v30 house style at the placed size | v28f producer size 7.1 x 8.8 in |
| `V30_PEAK_COLOR`, `V30_ION_PT`, `V30_ION2_PT`, `V30_LEADER_COLOR` | #b5b5b5; 7 pt; 7.2 pt; .5 | v30 drawing: grey peaks, ion labels (charge-2 labels at 7.2 pt so the superscript is at least 5 pt), grey leader lines | none |
| `V30_LETTER_PT`, `V30_DPI` | 7 pt, 300 dpi | house panel letters; PNG = PDF rasterised at 300 dpi | none |
| `CASES` | HIS_001 to HIS_004 | one case per origin x side in PREDICT stool (A25 selection rule) | strata table `spectra/all_strata.tsv` |

## Known limitations

- Selection is illustrative. Four spectra are not four independent estimates of a cleavage
  mechanism, and no test is made.
- The cases were chosen by a fixed rule (median acquisition, then median matched-intensity
  PSM), not for visual quality. The notebook confirms the rule against the selection tables;
  the step before (which peptides were candidates) depends on upstream site-evidence tables
  that are not in this repository.
- The bounded alternatives cover single substitutions, adjacent swaps and Sage C/M
  modification placements present in the local FASTA. They do not exclude every sequence or
  modification, chimeric spectra, or I/L ambiguity.
- For N-terminal cases (b, d) the histidine is in the source protein just before the
  peptide, not in the observed peptide.
- The v28f image keeps the producer's DejaVu Sans font, panel titles and colours. From v30 the
  figure is the house-style redraw (Söhne; panel titles in the x-axis labels); its plotted values
  are checked identical to the producer drawing.
