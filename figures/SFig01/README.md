# Supplementary Figure 1

Shows how small the acquisition-specific databases selected by FastaLake are, compared with the
full sequence reservoirs they were selected from and with the human reference proteome. Each blue
point is the number of FASTA entries in one acquisition-specific database; bars are cohort medians,
orange diamonds the full cohort reservoirs, and the dashed line the 20,431-entry human Swiss-Prot
proteome, on a logarithmic axis. The unit builds the v30 figure: the same values redrawn in the
house style (Söhne, ink, Okabe-Ito colours by meaning, drawn at 180 x 105.2 mm), proven identical
to the original producer and pixel-identical to the v30 image.

## Panels

The figure has one unlettered panel; this repository calls it `a`.

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | FASTA entries per acquisition-specific database by cohort, cohort medians, full reservoir sizes, human reference line | all databases smaller than the human reference (55); reservoirs of 233.5, 272.9, 264.9 and 375.3 million sequences (57) | 726 acquisitions (15 CAMPI, 58 healthy stool, 306 PREDICT, 347 CapScan); one reservoir per cohort | none (descriptive) | Panel a (producer `panel_a`, v30 `panel_a_v30`) | `source_data/Fig1/acquisition_sizes_726.tsv`, `source_data/SFig01/human_swissprot_canonical_2026_03.fasta.gz` | cohort_database_sizes (not in `pipelines/`) |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `HUMAN_REFERENCE_ENTRIES` | 20431 | entry count of the reviewed human Swiss-Prot canonical proteome, release 2026_03, isoforms excluded (Methods, v28d line 673); recounted from the shipped FASTA | none |
| `Y_LIMITS` | (150, 1.2e9), log scale | sizes span 241 to 375 million entries (Methods, v28d line 675) | none |
| `JITTER_SEED` | 20260928 | horizontal point jitter, display only; no value depends on it | none; the seed is not stated in the text |
| `JITTER_HALF_WIDTH` | 0.14 | horizontal spread of points, display only | none |
| `COHORT_ORDER` | CAMPI, HSTOOL, PREDICT, CapScan | producer order; CapScan last, so adding its 16 participant-15 rows leaves the jitter of the other cohorts unchanged | none |
| `FIGSIZE_IN` | (7.2, 4.2) | producer figure size | none |
| `PLACED_WIDTH_MM`, `PLACED_HEIGHT_MM` | 180.34, 105.20 | placed size of the image in the manuscript | none |
| `RASTER_DPI` | 400 | effective resolution of the placed v28f image (producer raster step) | none |
| `V30_PLACED_W_MM`, `V30_PLACED_H_MM` | 180.0, 105.2 | v30 figure drawn at its placed size, no scaling | none |
| `V30_PNG_DPI` | 300 | v30 PNG rasterised from the PDF at 300 dpi | none |

## Known limitations

- Reservoir sizes are read from the `reservoir_entries` column of the size table (one constant
  value per cohort). They were counted by the upstream size audit (cohort_database_sizes) from the
  reservoir statistics files; the reservoirs themselves are not in this repository.
- The 233,497,702-entry public-sequence benchmark reservoir is traced as the healthy-stool
  reservoir value; the text (v28d line 57) states that the public-sequence benchmarks used that
  reservoir.
- The producer drew with matplotlib defaults (typeface DejaVu Sans, which ships with matplotlib),
  not with the house style, so the pixel checks do not depend on the Söhne fonts and run always.
  `panel_a` switches off the house style set by `fl.init()` to reproduce this.
- The v28d legend and Methods say 710 acquisitions and 331 CapScan acquisitions; the v30 text
  (v30f) prints 726 and 347, as the figure plots.
- The v30 pixel check needs the Söhne fonts (`FASTALAKE_FONT_DIR`); the producer checks do not.
