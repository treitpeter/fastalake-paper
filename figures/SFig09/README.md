# Supplementary Figure 9

The figure relates the FastaLake peptide gain in CapScan to the size of the specimen-matched
metagenome (MG) database it is compared with. Gain is 100 x (FastaLake peptides / MG-search
peptides - 1) at 9-50 residues (v30), with FastaLake databases selected from the scored
spectra. Panel a shows specimens, panel b medians over intestinal-fluid acquisitions in thirds
of MG database size. Participant 15 has no MG database, so the 347-acquisition roster adds no
matched acquisition; the figure will change again with the 295-acquisition CapScan MG benchmark.

## Panels

Claims are cited by v28d line; printed values are those of the v30 text (v30f), in which the
peptide counts are 9-50 residues.

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | specimen median gain (9-50 residues) against predicted-protein records in the specimen-matched MG database | legend, 250 matched specimens (445) | 250 specimens (282 acquisitions, 15 participants) | none | Panel a | `source_data/SFig09/capscan_acquisitions_9_50.tsv` (aggregated by `derive_tables`) | capscan_metagenome_yield_join; capscan_public_only_search |
| b | median acquisition-level gain in thirds of intestinal-fluid MG database size | v30: gains 70%, 58%, 42% (129); 70.4%, 58.3%, 42.0% (445) | 226 acquisitions (76, 75, 75; 14 participants) | none (descriptive medians) | Panel b | `source_data/SFig09/capscan_acquisitions_9_50.tsv`, `capsule_assembly_gain_tertiles_9_50.tsv` (expected values) | capscan_metagenome_yield_join; capscan_public_only_search |

Both panels are one canvas (`build/SFig09/SupplementaryFigure9.{png,pdf}`), checked against the
v30 render. The notebook also recomputes the related numbers on v28d line 129 (MG database
sizes, DNA quantities, within-participant comparison, pooled yields at 9-50) from the same
table. The input table is a parameter (`ACQ_9_50_TSV`, environment `SFIG09_ACQ_TSV`): the
figure will be redrawn when the CapScan MG benchmark is extended to 295 acquisitions.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `GAIN_SCALE` | 100 | gain in percent, 100 x (FL/MG - 1) (legend, line 445) | none |
| `SPECIMEN_SUMMARY` | median | repeated acquisitions of a specimen summarised by their median (line 445) | none |
| `N_BINS` | 3 | "increasing thirds" of MG database size (lines 129, 445) | none reported |
| `BIN_LABELS` | Smallest, Middle, Largest | bin labels of panel b | none |
| `STOOL` | stool | stool acquisitions are excluded from panel b; bins are over intestinal fluid only | none |
| `LINTHRESH` | 50 | gain axis linear up to 50% and logarithmic above (line 445) | none |
| `N_MATCHED_ACQUISITIONS` | 282 | matched CapScan acquisitions with an MG database (line 129) | 347-acquisition roster: participant 15 has no MG database |
| `N_SPECIMENS` | 250 | specimens behind the 282 acquisitions (line 445) | none |
| `N_FLUID_ACQUISITIONS` | 226 | intestinal-fluid acquisitions in panel b | none |
| `EXPECTED_BIN_N` | 76, 75, 75 | bin sizes from pandas.qcut on 226 values (line 445) | none |
| `REGION_COLORS`, `REGION_LABELS` | producer colours; Region 1-4, Stool | producer reference render (checked) | none |
| `REGION_COLOUR` | Region 1 sky blue, 2 blue, 3 green, 4 reddish purple, Stool grey (Okabe-Ito) | v30 house style; orange is reserved for MG in the main figures | producer colours |
| `ACQ_9_50_TSV`, `TERTILES_9_50_TSV`, `REFERENCE_PNG` | `capscan_acquisitions_9_50.tsv`, `capsule_assembly_gain_tertiles_9_50.tsv`, `expected/SFig09_v29_reference.png` | inputs of the v30 figure; overridable by `SFIG09_ACQ_TSV`, `SFIG09_TERTILES_TSV`, `SFIG09_BASE_PNG` for the 295-acquisition rerun | none |
| `MARKER_S` | 11 pt^2 | producer marker size; also sets the lowered panel-a y minimum (-1.07%) so the -0.03% specimen shows | none |
| `V30_W_MM`, `V30_H_MM` | 180 x 91.2 mm | placed size | none |
| `V30_TICK_PT`, `V30_LAB_PT`, `V30_LEG_PT`, `V30_VAL_PT`, `V30_N_PT`, `LOG_EXP_PT` | 6.5, 7.0, 6.5, 6.5, 6.0, 5.0 pt | v30 type sizes; log tick exponents drawn at 5 pt | mathtext exponent would be 4.6 pt |
| `FIGSIZE`, `DPI` | 7.2 x 3.65 in, 300 dpi | as saved by the producer | none |
| `FONT_FILE` | soehne-buch.ttf | the font file the producer registered | DejaVu Sans without `FASTALAKE_FONT_DIR` |

## Known limitations

- Panel b medians are descriptive. Acquisitions are repeated within specimens and
  participants, and no interval or test is given.
- MG database size and gain share the MG search as comparator, so the association is not
  evidence of cause (stated in the legend).
- The table joining peptide counts, MG database size, BioSample, region and DNA quantity per
  acquisition was made upstream without saved code. Peptide counts, regions and MG database
  sizes are cross-checked here against independent tables; BioSample assignment, participant
  numbers and DNA quantities are not (see `OPEN_ITEMS.md`).
- The v30 figure is restyled (Söhne, Okabe-Ito region colours, panel titles removed); the
  producer drawing is run alongside and every plotted value is compared exactly.
- At 9-50 residues one specimen has a gain of -0.03%; the panel-a gain axis starts at -1.07%
  so its marker shows in full.
- SF09 will change again when the CapScan MG benchmark is extended to 295 acquisitions; the
  input paths are parameters for that rerun (`OPEN_ITEMS.md`).
