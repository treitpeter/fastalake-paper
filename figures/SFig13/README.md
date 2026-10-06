# Supplementary Figure 13

Shows the sequence context of histidine-ending peptide boundaries in CapScan and the replication of
C-terminal histidine enrichment in the PREDICT stool cohort. Panels a-f are heat maps of how often
each residue occurs at flank positions -4..+4 of boundaries whose P1 residue is histidine, compared
with the same source proteins, as participant means in percentage points, for proximal fluid, distal
fluid and stool and for N- and C-terminal boundaries; dots mark q < 0.05 in the complete 32,400-test
motif family. Panels g,h show PREDICT histidine enrichment at N- and C-terminal boundaries as 93
patient means, for de novo predictions and for the non-specific search.

This is the v28f version (347 CapScan acquisitions, 15 participants). It replaces the v28d image
(331 acquisitions, 14 participants in proximal and distal fluid); panels g,h are unchanged.

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | CapScan proximal fluid, N boundary, H-conditioned motif deviations | legend: dots q<0.05 in the 32,400-test family, 840 cells displayed (461) | 15 participants per cell | two-sided exact sign-flip, BH in 32,400 tests | Panel a | `source_data/SFig13/capscan_histidine_motif_display.tsv` (from `source_data/Fig6/capscan_histidine_statistical_families.tsv`); family recomputed from `source_data/SFig13/capscan_participant_motif_differences.tsv.gz` | terminal_histidine_inference |
| b | as a, proximal fluid, C boundary | 461 | 15 participants | as a | Panel b | as a | terminal_histidine_inference |
| c | as a, distal fluid, N boundary | 461 | 15 participants | as a | Panel c | as a | terminal_histidine_inference |
| d | as a, distal fluid, C boundary | 461 | 15 participants | as a | Panel d | as a | terminal_histidine_inference |
| e | as a, stool, N boundary | 461 | 15 participants | as a | Panel e | as a | terminal_histidine_inference |
| f | as a, stool, C boundary | 461 | 15 participants | as a | Panel f | as a | terminal_histidine_inference |
| g | PREDICT H enrichment, de novo top 30%, N and C boundaries: patient means, medians, mean and 95% bootstrap interval | all 93 patients positive with both routes (175); 93 patient means (461) | 93 patients per side | two-sided sign-flip (200,000 assignments), BH in 1,560 tests; percentile bootstrap | Panel g | `source_data/SFig13/predict_histidine_acquisition_enrichment.tsv.gz`, `source_data/SFig13/predict_histidine_statistical_families.tsv` | terminal_histidine_inference |
| h | as g, non-specific Sage with peptide and spectrum q <= 0.01 | 175; 461 | 93 patients per side | as g | Panel h | as g | terminal_histidine_inference |

## v30 composite

From manuscript v30 on, the figure is the house-style redraw of the same values (`panels.py:figure_v30`, notebook section `v30 composite (house style)`; ASSEMBLY.md, second part). Panels and lettering are unchanged; only type, colours, size and the placement of panel titles changed. The function compares every plotted array with the producer drawing and raises on any difference. The regenerated PNG is byte-identical to the v30 image (`checks/SFig13_v30.json`).

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `DISPLAY_ROUTE` | `denovo_top30` | displayed motif route, de novo top-30% predictions (Supplementary Methods 9, v28d line 867) | `sage_semi_bothq` and `sage_unspecific_bothq` are in the primary flank table and all six routes in the family; Source Data only |
| `DISPLAY_ORIGIN` | `All` | all peptide origins pooled | `Human`, `Microbial` in the family; Source Data only |
| `DISPLAY_POPULATION` | `H_only` | boundaries whose P1 residue is histidine (legend, v28d line 461) | `all_nonKR` in the family; Source Data only |
| `CONDITIONED_POSITION` | -1 | the fixed histidine position is omitted from display and tests (v28d lines 461, 873) | none |
| `Q_CUTOFF` | 0.05 | BH q threshold for dots (v28d lines 461, 889) | none |
| `MOTIF_COMPARTMENTS` | proximal = duodenum + jejunum, distal = ileum + colon, stool | compartments of panels a-f; a participant counts only with all regions of the compartment | single regions in the participant table; not tested |
| `PREDICT_ROUTES` | `denovo_top30`, `sage_unspecific_bothq` | the two PREDICT boundary routes (v28d line 869) | `sage_unspecific` (peptide q only) and score bins in the 1,560-test family; Source Data only |
| `PREDICT_ORIGIN`, `PREDICT_RESIDUE` | `All`, `H` | histidine, all origins | other residues and origins in the family |
| `SIGNFLIP_MIN_N` | 8 | n < 8 is not estimable, P = 1 (A28 `signflip`) | none |
| `SIGNFLIP_EXACT_MAX_N` | 16 | n <= 16 enumerates all sign assignments (v28d line 869) | none |
| `SIGNFLIP_DRAWS` | 200,000 | random sign assignments for n > 16, plus-one correction (v28d line 869) | none |
| `SIGNFLIP_SEED_BASE` | 20260914 | seed = base + n; not stated in the text | none |
| `BOOTSTRAP_DRAWS` | 20,000 | whole-unit bootstrap resamples (v28d line 869) | none |
| `BOOTSTRAP_SEED_BASE` | 20261014 | seed = base + n; not stated in the text | none |
| `PRODUCER_RC`, `COL`, `FIGSIZE_IN`, `HEIGHT_RATIOS` | as in `panels.py` | producer drawing settings (`build_revision.py`), DejaVu Sans | none |
| `V30_SIZE_MM`, `V30_HEIGHT_RATIOS`, `V30_W_PAD`, `V30_H_PAD`, `V30_HSPACE`, `V30_WSPACE` | 160 x 243.4 mm; 1.1, 1.1, 1.1, 0.95; constrained-layout pads | v30 house style at the placed size (v30 restyle script) | v28f producer size 7.1 x 10.8 in |
| `V30_REGION_LABEL`, `V30_ROUTE_LAB` | "Regions 1–2", "Regions 3–4"; "PREDICT: de novo top 30%", "PREDICT: unspecific Sage; both q gates" | v29 panel-title wording, moved to y-axis labels (a-f) and legend entries (g,h) | none |
| `V30_LETTER_PT`, `V30_DPI` | 7 pt, 300 dpi | house panel letters; PNG = PDF rasterised at 300 dpi | none |
| `RASTER_DPI` | 300 | raster of the redacted PDF (producer `prepare_supplement_figures.py`) | none |
| `PANEL_CROP_MARGIN_IN` | 0.06 | margin of the per-panel crops (display and per-panel checks only) | none |

## Known limitations

- The producer drew all eight panels as one matplotlib figure. The repository does the same; the
  panel sections show crops of the composite, and the per-panel pixel checks compare crops cut with
  the same box from the repository render and from the producer's PNG. The whole composite (both
  PDFs at 300 dpi and the 300 dpi PNG) is also checked and is identical; the PNG is byte-identical.
- The plotted means and intervals of panels g,h are read from the published PREDICT family table,
  as the producer did. The notebook recomputes the whole 1,560-test family from the acquisition table
  and asserts agreement within 1.1 x 10^-16.
- The CapScan motif family is recomputed from participant-region motif differences. The upstream
  per-acquisition motif counts and source-protein backgrounds (A28 collection over the acquisition
  outputs) are not rerun here.
- The PREDICT acquisition table carries pseudonyms. Visit codes are withheld; the analysis averages
  over all specimens of a patient, so no value depends on them.
- The PREDICT panels are labelled "unspecific Sage" in the figure (v30: legend entry of h), while
  the text says "non-specific"; the wording is not changed.
