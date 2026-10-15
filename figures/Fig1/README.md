# Figure 1

Figure 1 shows how FastaLake selects a compact search database from a large protein
sequence reservoir by exact lookup of de novo peptide predictions, and how specific that
lookup is. Panel a is a schematic that compares the three database strategies and follows
one healthy-stool acquisition, and its scrambled-query control, from the 233.5-million
sequence reservoir to accepted peptides. Panels b to e give the sizes of all 726 selected
databases, the chance-match rate of scrambled queries by length, how many selected
proteins rest on a single query, and the donor-selected entrapment estimate.

The figure in the v33 manuscript is `expected/Fig1_v33.png` (`artwork/Figure1_v33.pdf`); v30 used `artwork/Figure1_v30.pdf`, kept as an input of the funnel build.
`NUMBERS.tsv` keeps the v28d line numbers and quotes and carries the v30 paragraph, wording,
printed value and status in the appended `v30_*` columns (values unchanged).

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | database strategies; worked acquisition and scrambled control as one funnel (v30), stage counts in the table | Results 51, 53, 59; legend 71 | one acquisition; 58 acquisitions pooled for the study-level lookup | artwork | artwork (`Panel a`, numbers recomputed; FastaLake column rebuilt by `artwork/build_fig1a_funnel.py` and pixel-checked) | artwork, see `ARTWORK_NUMBERS.md` | scramble_end_to_end_workflow, healthy_stool_stage1_lookup, healthy_stool_selection_sweep |
| b | entries per acquisition-specific FASTA, by cohort; human reference line | 55; legend 71 | 726 acquisitions (15 CAMPI, 58 healthy stool, 306 PREDICT, 347 CapScan) | none (medians) | `Panel b` | `acquisition_sizes_726.tsv`, `human_reference.json` | cohort_database_sizes |
| c | fraction of distinct authentic or scrambled queries with an exact reservoir match, by length | 59; legend 71 | distinct queries of one acquisition per seed and length; five scramble seeds | none (range over seeds) | `Panel c` | `length_by_arm.tsv`, `five_seed_summary.tsv` | database_selection_specificity |
| d | per acquisition, fraction of selected proteins with exactly one production query | 63; legend 71 | 726 acquisitions | none (medians) | `Panel d` | `acquisition_support_summary_726.tsv`, `worked_acquisition_prediction_support.json` | cohort_prediction_support |
| e | donor-selected entrapment paired estimate, 5-50 and 9-50 residues | 65; legend 71 | 58 acquisitions, each the mean of three seeds | none (medians) | `Panel e` | `donor_seed_means.tsv`, `donor_recount.tsv` | healthy_stool_donor_entrapment |

All source_data files are in `source_data/Fig1/`; their columns are described in
`source_data/Fig1/README.md`.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `Q_CUTOFF` | 0.01 (peptide q) | acceptance threshold used throughout (Results line 53; Methods line 727) | q 0.005 and 0.001 and PEP gates in `scramble_end_to_end.tsv`; 0.001 to 0.05 in `donor_recount.tsv` |
| `ENTRAPMENT_TIER` | `T1_within_group` | donor database from another acquisition of the same design (Methods line 727) | self-selected `T0_self` (five acquisitions) in `donor_seed_means.tsv`, Supplementary Fig. 5d |
| `ENTRAPMENT_LENGTHS` | 5 and 9 (to 50 residues) | 9-50 is the lookup query length; 5-50 the search range (Methods line 729) | none |
| `ENTRAPMENT_SEEDS` | 20260910, 20260911, 20260912 | three seeds averaged per acquisition (Methods line 729) | per-seed values in `donor_recount.tsv` |
| `SCRAMBLE_SEEDS` | 20260901, 20260913, 20260914, 20260915, 20260916 | Methods line 665 | terminal-preserving design with the same seeds, `scramble_end_to_end.tsv` |
| `WORKED_SEED` | 20260901 | seed propagated through the full workflow in panel a (Methods line 665) | the other nine workflows in `scramble_end_to_end.tsv` |
| `WORKED_ACQUISITION` | 20181231_QX7_PhGe_LC95_SA_WorkflowReplicates_2 | representative healthy-stool acquisition (Results line 53) | all 58 in `selection_grid_acquisitions.tsv` |
| `TOP_FRACTION_ARM` | `tp0.30` | top 30% of eligible predictions per acquisition (Methods line 557) | 5-100% in `selection_grid_acquisitions.tsv` (Supplementary Fig. 2) |
| `SUPPORT_POPULATION` | `all_selected` | all proteins in the selected FASTA (Methods line 677) | `identified_accessions`, `native_protein_q01` in the same table |
| `COHORTS`, `ROSTER` | CAMPI 15, healthy stool 58, PREDICT 306, CapScan 347 | v28f roster, 726 acquisitions (legend line 71) | 710-acquisition roster (CapScan 331) in the superseded Source Data copy |
| `SCRAMBLED_MIN_LENGTH_NO_MATCH` | 13 | ">= 13 residues" (Results line 59) | full curve in `length_by_arm.tsv` |
| `JITTER_SEED_BASE` | 610 | horizontal point jitter of panels b, d, e, from the producer | none; jitter does not change any value |
| `DPI` | 300 | producer PNG resolution and pixel-check resolution | none |

## Known limitations

- Panel a is not regenerated. It is Illustrator artwork changed by about ten successive
  PyMuPDF edit scripts (`ASSEMBLY.md`). Its typed numbers are recomputed from tables
  (`ARTWORK_NUMBERS.md`).
- Panel e in the v28f composite is the v10 render (18 September 2026) of the same plot
  code and data. It is checked as a crop of that render (`expected/panel_e.png`), because
  the producer saved only the four-axis row.
- The plotting functions apply the producers' style module values (0.6 pt lines, regular
  panel letters, tight bounding box) on top of `style/fl_style.py`, which is a later
  revision. See `panels.py` and `OPEN_ITEMS.md`.
- The worked-acquisition counts (panel a, Results 53, 59, 61) come from one acquisition;
  they are descriptive, not estimates over acquisitions.
- The count of 4,258 accepted peptides absent from the predictions is read from the
  upstream record and was re-derived independently from the raw predictions and Sage
  output on Slurm (`PROVENANCE.md`); those raw inputs are not in this repository.
- Several printed values are rounded to two or three significant figures in the text
  (12,600 for 12,606); `tools/check_numbers.py` does not yet accept that rounding
  (`OPEN_ITEMS.md`).
