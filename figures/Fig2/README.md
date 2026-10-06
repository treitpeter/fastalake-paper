# Figure 2

Figure 2 tests whether FastaLake selects a database that matches a community it was not told
about. Fifteen CAMPI SIHUMIx acquisitions of a defined eight-member community were searched
against five databases: FastaLake (FL), the supplied member reference proteomes (REF, printed
"Reference" in v28d), the metagenome co-assembly proteins (MG), and FL with REF appended alone or
together with MG. The panels compare database size, peptide and protein-group yield, where
reference peptides were lost, and which members the identified peptides belong to.

This unit follows the v30 figure. The producer draws all eight panels on one 180 x 250 mm canvas;
in v30 panels b and d-g show peptides of 9-50 residues and c and h of 5-50 residues (len950
rebuild), and the panel b and e axis titles read "accepted". `panels.py` holds the v30 producer code
unchanged except for paths; `assembly.py` holds the six PDF edit steps that turn the canvas into the
composite (labels and removed elements; no replotting). The notebook redraws the canvas from
`source_data/`, runs the edits and compares the composite with the v30 figure: IDENTICAL, 0 differing
pixels (`checks/composite.json`). `NUMBERS.tsv` traces every number of the v28d and the v30 text and legend (v30 columns
`v30_para`, `v30_quote`, `v30_printed_value`, `v30_status`).

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | community artwork and the five databases with median protein entries | median 3,100 entries vs 29,700 REF and 19,400 MG (77; v30 para 37) | 15 acquisitions per database | none (medians) | Panel a | `Fig2/metrics_per_acquisition.tsv`; artwork: `ARTWORK_NUMBERS.md` | campi_database_metrics |
| b | accepted 9-50-residue peptides per acquisition, paired across databases | median 17,100 vs 18,000 and 18,200 (77; v30 para 37) | 15 acquisitions per database (repeated) | none (medians) | Panel b | `Fig2/paired_reference_recovery.tsv` (min_length 9) | campi_reference_gap |
| c | quantified protein groups per acquisition (5-50); median % groups with >=2 assigned peptides | median 2,936 vs 3,968 and 4,048; 4,356 and 4,454 after add-back, higher in all 15 (77; v30 para 37) | 15 acquisitions per database (repeated) | none (medians, counts) | Panel c | `Fig2/arm_<cohort>__reaggregated_per_acquisition.tsv`, `Fig2/arm_<cohort>__direct_acquisition_summary.tsv` | campi_simplified_grouping |
| d | stage at which 9-50-residue reference peptide-acquisition observations missed by FL were lost | 94.2% of 44,217 lacked the sequence in the selected database (79; v30 para 38) | 44,217 peptide-acquisition observations, 15 acquisitions | none | Panel d | `Fig2/stage_summary.tsv` (min_length 9) | campi_reference_gap |
| e | pooled 9-50-residue peptides by member assignment; pooled recovery of 63,348 reference identities | 84.6% pooled, 98.6% with REF appended, 98.1% with MG (79; v30 para 38) | distinct peptides pooled over 15 acquisitions | none | Panel e | `Fig2/recovery_summary.tsv` (min_length 9) | campi_reference_gap |
| f | pooled member-specific 9-50-residue peptides per member and database | all eight members for FL, six for MG; 2,538 vs 217 for A. caccae (83; v30 para 40) | distinct pooled peptides; 8 members | none | Panel f | `Fig2/member_summary.tsv` (min_length 9) | campi_reference_gap |
| g | sharing of 9-50-residue FL peptides among members (UpSet, 14 sets) | 84% of 9-50-residue FL peptides map to one member (83; v30 para 40) | 55,100 FL peptides in the drawn sets | none | Panel g | `Fig2/upset_intersections_FL_9_50.tsv`, `SFig03/CAMPI_membership_<sample>.tsv` | campi_reference_gap |
| h | FL peptide length (5-50) by member-assignment class | median 7 versus 13 residues (83; v30 para 40) | 63,640; 3,842; 10,160 FL peptides | none (medians) | Panel h | `Fig2/length_by_member_class_a14.tsv` | campi_reference_gap |

Text statistics of this unit that are not drawn in a panel (all recomputed in the notebook's last
cell): reservoir size and member containment (line 75; cohort_database_sizes, campi_reservoir_containment),
add-back mechanism and depth standardization (line 77; campi_variability), same-spectrum audit
(line 81; campi_near_miss_audit) and AlphaNovo accuracy (line 81; campi_alphanovo_accuracy).

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `MIN_LENGTH` | 5 | 5-50 window for panels c and h and the add-back and subsampling statements (v30 legend, para 42; para 37) | 9-50 (`PLOT_MIN_LENGTH`) |
| `PLOT_MIN_LENGTH` | 9 | v30 window for panels b and d-g (legend, para 42); the min_length 9 rows are the accepted 5-50 lists filtered to length >= 9 | 5 reproduces the v28f canvas (`expected/Fig2_v28f.png`) |
| `MIN_LENGTH_SENSITIVITY` | 9 | 9-50 window of the member statements (v30 para 40) | 5-50 in v28d |
| `MAX_LENGTH` | 50 | upper peptide length (line 87) | none |
| `GROUPING_METHOD` | simplified | protein groups from the simplified grouping used throughout the manuscript; checked equal to the direct export | the `original` rows in the same tables (not shown) |
| `N_LEAST_REPRESENTED` | 3 | "three least represented members" (v30 para 40), the members with fewest REF-search member-specific 9-50-residue peptides | none |
| `N_MOST_REPRESENTED` | 5 | "four of the five most represented members" (v30 para 40), ranked by REF-search member-specific 9-50-residue peptides | none |
| `MG_MISSED_MEMBER` | Anaerostipes caccae | the member "which MG largely missed" (v30 para 40), excluded from the four | none |
| `COMPOSITE_WIDTH_MM` (`assembly.py`) | 183 | composite width (assembly step 3) | none |
| `COMPOSITE_DPI` (`assembly.py`) | 300 | composite raster, as the v30 PNG | none |
| `JITTER_SEED` | 20260922 | producer's seed for horizontal point jitter in b and c; display only | none |
| `JITTER_HALF_WIDTH` | 0.1 | producer's jitter width; display only | none |
| `PAGE_MM` | (180, 250) | producer canvas size | v28d page is 183 x 250.9 mm after scaling (`ASSEMBLY.md` step 3) |
| `BOTTLE_RECT_MM` | (3, 7, 50, 40.5) | producer placement of the community artwork | none |
| `PRODUCER_LINE_PT` | 0.6 | line weight of the house style module the producer imported; the repository style module later raised it to 1 pt | none |
| `RENDER_DPI` | 300 | resolution of the producer's PNG and of the pixel checks | none |

## Known limitations

- The producer canvas carries the review letters (a, b, c, d-f, g, h, i), the panel f orange frames
  and "10 or more" row, and several panel a text lines that the v28d figure no longer shows. These
  were changed by PDF edits after plotting. `assembly.py` re-executes those edits as the v30 producer
  ran them (`ASSEMBLY.md`).
- Panels c and h stay at 5-50 residues while b and d-g are 9-50 (as the v30 legend states).
- Medians in b and c, and all text numbers, are descriptive. No hypothesis test is made in this figure.
- The fifteen acquisitions are not independent biological replicates: S11 is one sample measured
  as four fractions, and other acquisitions repeat the same community under different conditions.
- Pixel identity of text requires the licensed Söhne fonts (`FASTALAKE_FONT_DIR`).
- Some printed values are rounded to hundreds (3,100; 29,700; 19,400; 17,100; 18,000; 18,200); see `OPEN_ITEMS.md`.
