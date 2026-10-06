# Supplementary Figure 4

Peptide-length sensitivity of FastaLake results. Every production search was repeated
with a minimum peptide length of seven residues, with confidence re-estimated in each new
search. The figure compares these searches, and the production results filtered to at
least 7 or 9 residues, with the published 5-50-residue results: the gain over the matched
metagenome (a, b), where short peptides were accepted (c), the weakest CAMPI members (d),
the CAMPI fusion control (e), healthy-stool search strategies (f) and protein grouping (g).
The 347-acquisition CapScan roster is used in c and g. This unit reproduces the v30 figure: the same values as v29, restyled (house style, no subpanel titles, colours by meaning, 7-50 / 9-50 labels in place of >=7 / >=9).

The composite is drawn as one figure by `panels.build_figure`, moved from the v30 build.
Panel images are crops of the composite at boxes computed from the drawn axes
(`panels.panel_boxes`). The regenerated 300 dpi composite is byte-identical to the v30
render (`checks/SFig04_v30.json`).

## Panels

| panel | what it shows | claim (v30 paragraph) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | per-acquisition gain of public FastaLake over the matched metagenome, by length population | 398 ("reduced the median gains to 14.5% (PREDICT) and 45.0% (CapScan)"; "gains of 67% and 96%") | 296 PREDICT, 282 CapScan acquisitions | none (medians) | Panel a | benchmark_matched_per_acquisition.tsv | minlen7_search_sensitivity |
| b | recovery of metagenome-search identities | 398 ("identity recovery remained 85.4% and 84.0%") | 296, 282 acquisitions | none (medians) | Panel b | benchmark_matched_per_acquisition.tsv | minlen7_search_sensitivity |
| c | median change in accepted peptides vs 5-50, filtered or re-searched at 7-50, by database | 225 (legend) | 15 to 347 acquisitions per database (CapScan FL 347) | none (medians) | Panel c | per_acquisition_length_counts.tsv, per_acquisition_length_counts_participant15.tsv | minlen7_search_sensitivity |
| d | pooled member-specific peptides, three least represented CAMPI members | 225 (legend) | peptides pooled over 15 acquisitions | none | Panel d | campi_member_counts.tsv | minlen7_search_sensitivity |
| e | CAMPI fusion-control diagnostic at q <= 0.01, per seed | 225 (legend) | 15 acquisitions per seed, 2 seeds | none (pooled ratio) | Panel e | campi_fusion_pooled_q01.tsv, campi_fusion_per_acquisition_q01.tsv | minlen7_search_sensitivity |
| f | healthy-stool common-FASTA and joint searches vs own database | 225 (legend) | 58 acquisitions | none (medians) | Panel f | healthy_strategies_per_acquisition.tsv | minlen7_search_sensitivity |
| g | change in study-level protein groups, new 7-50 vs 5-50 | 126, 225 | one grouping per cohort and database; CapScan 347 acquisitions | none | Panel g | grouping_study_groups.tsv | minlen7_search_sensitivity |

The legend search count (2,022) is computed from `search_inventory.tsv`
(`panels.stat_search_count`).

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| CAPSCAN_ROSTER | 347 | v28f CapScan roster (15 participants); c CapScan FL and g CapScan use it | 331-acquisition version is the v28d figure (Codex census image11) |
| MATCHED_N | PREDICT 296, CapScan 282 | matched FL_public vs MG pairs; participant 15 has no metagenome | none |
| POPS_AB | baseline5, fresh7, nested9 | populations drawn in a and b | nested7 and fresh7_nested9 are in the source table |
| GAIN_CLIP | 250 | legend: values above 250% drawn at 250% and counted | none |
| LOG_OFFSET | 0.5 | legend: 0.5 added for the logarithmic axis in d | none |
| FUSION_Q | 0.01 | legend: peptide q <= 0.01 in e | other thresholds in the fusion run receipts |
| FUSION_REF_LINE | 1 | legend: dashed line at 1% in e | none |
| D_MEMBERS | L. plantarum, C. butyricum, B. longum | three least represented CAMPI members | all eight members in campi_member_counts.tsv |
| JITTER_SEED | 1 | horizontal jitter of points in a and b (display only) | none |
| JITTER_HALF_WIDTH | 0.25 | jitter width in a and b | none |
| SLANT_SHEAR | 0.20 | synthetic italic for the taxon labels in d (no Söhne italic installed) | none |
| PLACED_W_MM, PLACED_H_MM | 160, 208 | size at which the figure is placed in the manuscript | none |
| PNG_DPI | 300 | v30 raster render of the PDF (1890 x 2457 px) | none |
| C (colours) | 5-50 grey, 7-50 filtered light grey, 7-50 new search sky blue, 9-50 filtered vermilion, 9-50 subset reddish purple | v30 colours by meaning: FL blue and MG orange are kept for databases | v28f colours (blue, orange, purple) in the v28f figure |
| LAB (legend labels) | "7-50, filtered", "7-50, new search", "9-50, filtered", "9-50 subset of new search" | Söhne has no >= glyph; same meaning (maximum 50 residues) | ">=7" / ">=9" in v28f |
| LETTER_X_MM | 1, 82 | one panel-letter x per column | none |

## Known limitations

- Panels a and b stay on the matched 296 PREDICT and 282 CapScan pairs: participant 15 has
  no metagenome and no public-sequence database among the 282.
- All values are descriptive medians or counts; no test or interval is computed. Several
  acquisitions come from one participant, and this is not modelled.
- The 9-50 values of the production searches reuse the original q values (filtered, not
  re-estimated).
- Panel e is a fusion-linked diagnostic, not a certified false-discovery rate.
- Pixel-exact text needs the licensed Söhne fonts (`soehne-buch.otf`, `soehne-halbfett.otf` and the TrueType `soehne-buch.ttf` for the slanted taxon labels).
