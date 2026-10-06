# Figure 6

Figure 6 shows that peptides from intestinal fluid often end at internal, non-tryptic
protein positions, that histidine is strongly over-represented at those C termini in
host and microbial peptides, that this is not explained by the histidine content of the
source proteins, and that histidine boundaries more often carry a complementary product
on the other side. CapScan data are 347 acquisitions from 15 participants; the PREDICT
stool comparison has 93 patients. The notebook regenerates all six panels with the
producer's plot code, rebuilds the v28f composite, recomputes every test family behind
the figure, and checks every panel pixel by pixel against the original render
(`checks/`). The panel b key reads "Human-assigned" (the v30 term) from the plot code;
the last composite step flattens the page, as in the v30 manuscript figure
(`expected/Fig6_v30.png`). Panels a and c-f are identical to the producer renders; the
panel b key differs from the v30 figure only in glyph rasterisation (see Known limitations).

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | share of 9-50-residue peptides with a non-K/R C terminus, regions 1-4 and stool, three routes; dashed line, tryptic search | non-K/R termini 5.0-5.5% (v28f 5.1-5.6%) vs 3.1% stool, higher in all fourteen (v28f fifteen) participants; tryptic 0.33% (v28f 0.34%) (165, 169) | 15 participants per region (participant) | two-sided paired sign-flip test of capsule mean minus stool, exact; BH within `terminal_fraction__length_9` (4,536 tests) | Panel a | `capscan_terminal_participant_region.tsv`, `capscan_terminal_statistical_families.tsv` | terminal_participant_inference |
| b | histidine share of non-K/R C termini, human-assigned vs microbial, non-specific search (key "Human-assigned", written by the plot code) | higher in fluid in all participants, mean paired differences 15.2 and 18.0 pp (v28f 17.6) (169, 173) | 15 participants (participant) | as a | Panel b | as a | terminal_participant_inference |
| c | composition of non-K/R C termini, regions 1-2, regions 3-4, stool | histidine 40.3% and 40.8% (v28f 39.9%, 40.3%) vs 24.4% (173) | 15 participants per group (participant) | interval only: 95% percentile bootstrap, 2,000 participant resamples | Panel c | `capscan_residue_participant_region.tsv` | terminal_participant_inference |
| d | histidine at internal C boundaries, observed minus source expectation | about 40% observed vs about 2.4% expected; every region (175) | 15 participants (CapScan), 93 patients (PREDICT) | interval: 95% bootstrap, 20,000 whole-unit resamples; test (STATS.tsv): sign-flip, BH within `source_P1_all_residues` (9,000 CapScan; 1,560 PREDICT) | Panel d | `capscan_histidine_primary_source_enrichment.tsv`, `predict_histidine_primary_source_enrichment.tsv` | terminal_histidine_inference |
| e | as d for internal N boundaries | C enrichment greater than N (175) | as d | as d | Panel e | as d | terminal_histidine_inference |
| f | complementary-product rate, histidine minus matched non-histidine boundaries | more complementary products in fluid on all three routes (177) | as d | as d; BH within `all_reciprocal_controls` (90 CapScan; 12 PREDICT) | Panel f | `capscan_histidine_primary_reciprocal_tests.tsv`, `predict_histidine_primary_reciprocal_tests.tsv` | terminal_histidine_inference |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `MIN_LENGTH` | 9 (upper bound 50 applied upstream) | terminus analyses use 9-50 residues (lines 285, 855) | lengths 0 and 5 families in `capscan_terminal_statistical_families.tsv` |
| `TERMINUS` | `C_last` | C-terminal residue (lines 169, 855) | `N_first`, `N_cleavage_P1` rows kept in source data |
| `PRIMARY_ORIGINS` | All, Human, Microbial | conservative origin rule (lines 287, 861) | Mixed, Unknown, NotInSelectedFASTA kept in source data |
| `DISPLAY_ROUTES` | `denovo_top30`, `sage_semi_bothq`, `sage_unspecific_bothq` | displayed routes require peptide and spectrum q <= 0.01 for Sage (lines 285, 869) | peptide-only and score-bin routes in the A28 family table |
| `SEED_PLOT` | 6 | builder generator for point jitter (a, b) and panel c bootstrap | not stated in v28d |
| `N_BOOT_PANEL_C` | 2,000 | Fig. 6c intervals (line 869) | none |
| `A24_N_BOOT_TERMINAL`, `A24_N_BOOT_COMPOSITION` | 50,000; 2,000 | A24 interval draws (line 869) | none |
| `A24_SEED_SALT` | "20260914" | A24 seeds are the first 4 bytes of sha256 of the key plus this salt | not stated in v28d |
| `EXACT_MAX_N`, `N_SIGN_DRAWS` | 16; 200,000 | exact sign enumeration for n <= 16, otherwise random signs with plus-one (line 869) | none |
| `MIN_PAIRED_N` | 8 | fewer paired units: test kept with p = 1 | not stated in v28d |
| `A28_N_BOOT` | 20,000 | internal-boundary and complementary-product intervals (line 869) | none |
| `A28_SEED_NULL`, `A28_SEED_BOOT`, `A28_SEED_RATIO` | 20260914 + n, 20261014 + n, 20261114 + n | A28 generator seeds | not stated in v28d |
| `CI_QUANTILES` | (0.025, 0.975) | 95% percentile bootstrap intervals; used by every interval in panels.py (line 169) | none |
| `FINAL_WIDTH_MM` | 183 | double-column export | builder canvas 180 mm |
| `KEY_B_HUMAN` | "Human-assigned" | v30 term for human-origin peptides (v30 paragraph 87), written by `panel_b` in the legend | the v28f producer wrote "Host"; the v29 tool covered and retyped it on the composite (`KEY_B_OLD`, `KEY_B_NEW`, kept only for that optional path of `relabel_and_flatten`) |
| `FLATTEN_DPI` | 600 | v30 PDF holds one raster so the covered text is gone | none |
| `REGION_LABELS` | Duodenum..Colon to Region 1-4 | manuscript region naming (line 169) | none |

## Known limitations

- The panel b key "Human-assigned" is written by the plot code. The v30 manuscript figure
  was made by covering "Host" on the composite and retyping the new word with PyMuPDF.
  Against that figure (300 dpi PNG) 5,884 pixels differ, all inside the panel b key
  (bounding box x 1266-1882, y 535-589 px): matplotlib and PyMuPDF rasterise the glyphs
  differently, and the "Microbial" key row sits about 1 px lower because matplotlib
  measures the new first row. Wording, position, face and size agree; the diff image is
  `checks/Fig6_v30_png.diff.png`. For the same reason panel b, the builder composite and
  the v28f composite now differ from their v28f references (which read "Host") in the
  key only. The v30 PDF is a 600 dpi raster without a text layer.

- The repository style module (`style/fl_style.py`) has moved to 1 pt lines since this
  figure was drawn; `panels.producer_style` restores the producer's settings (0.6 pt
  lines, Söhne Buch only) so that the panels are reproduced, not restyled.
- Region names are replaced in the PDF by white boxes with new text over the original
  vector text, as in the producer. The original words remain in the PDF text layer. In
  panel d, "Region 3" and "Region 4" are wider than "Ileum" and "Colon" and touch the
  axis line.
- The per-panel crops in `expected/` are cut from the producer's builder PNG with the
  boxes in `PANEL_BOXES_MM`; the full builder output and the v28f composite are checked as
  whole files as well.
- The non-specific searches, origin assignment and boundary localisation (A05, A20, A27)
  were not re-run; the notebook starts from participant- and patient-level tables.
- PREDICT complementary-product values (panel f, PREDICT rows) are read from the shipped
  table; patient-level values for them are not in this repository.
