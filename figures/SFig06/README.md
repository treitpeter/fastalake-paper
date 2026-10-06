# Supplementary Figure 6

The figure asks how much FastaLake identification depends on selecting the protein
database from the same spectra that are then scored. In panels a-c, databases were
selected from all de novo predictions (all-prediction database), from odd-scan
predictions (same-half database, the same half as the scored spectra) or from even-scan
predictions (disjoint database), and all three were searched against the same odd-scan
MS2 spectra of 726 acquisitions. Panels d,e show the same-spectrum benchmark of
public-sequence FastaLake against the matched metagenome search (296 PREDICT and 282
CapScan acquisitions). All counts are peptides of 9-50 residues (v30).

## Panels

Claims cite paragraphs of the v30 text (`visible_v30f.txt`, first field).

| panel | what it shows | claim (v30 paragraph) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | all-prediction / same-half / disjoint selection, identical odd-scan spectra searched | design of the spectrum-split experiment (49) | none | none (schematic drawn in code) | Panel a | none (ARTWORK_NUMBERS.md) | cohort_spectrum_split_selection |
| b | disjoint/same-half peptide-count ratio, 9-50 residues | medians 92.4%, 82.9%, 86.5%, 79.9%; fewer in every acquisition; retention 70.0-83.7% (49); 13.5% and 20.1% loss (64) | 726 acquisitions: 15 CAMPI, 58 healthy stool, 306 PREDICT, 347 CapScan | none; median and IQR | Panel b | split_half_acquisitions_726.tsv | cohort_spectrum_split_selection |
| c | disjoint/all-prediction peptide-count ratio, 9-50 residues | legend (231) | 726 acquisitions, as b | none; median and IQR | Panel c | split_half_acquisitions_726.tsv | cohort_spectrum_split_selection |
| d | peptides per acquisition, MG vs FL, same spectra, 9-50 residues | median gain 9.4% with 272 of 296 higher (64); 42.4% with 281 of 282 higher (65) | 296 PREDICT, 282 CapScan acquisitions | none; medians | Panel d | same_spectrum_benchmark_0_predict.tsv, same_spectrum_benchmark_1_capscan.tsv | clinical_public_only_comparison |
| e | fraction of MG peptide identities also found by FL, 9-50 residues | 87.7% (64) and 86.4% (65) | 296 PREDICT, 282 CapScan acquisitions | none; medians | Panel e | as d | clinical_public_only_comparison |

All five panels are one canvas (`panels.figure`, 180 x 186 mm), as in the v30 build
(ASSEMBLY.md).

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `MIN_LENGTH`, `MAX_LENGTH` | 9, 50 residues (the upper bound is applied upstream) | headline peptide population of the v30 text (legend line 231) | 5-50 residues: second length view in the source tables; Source Data |
| `Q_CUTOFF` | 0.01, Sage peptide q per search | acceptance rule; applied upstream | none in this figure |
| `SPLIT_ROSTER_N`, `N_CAPSCAN` | 726 acquisitions, 347 CapScan | v30 roster (legend line 231) | 710-acquisition run of v28f (no longer used) |
| `SPLIT_TABLE` | split_half_acquisitions_726.tsv | b,c input, parameter of `figure` | none |
| `BENCH_TABLES`, `BENCH_N` | same_spectrum_benchmark_{0_predict,1_capscan}.tsv; 296 PREDICT, 282 CapScan | d,e inputs, parameters of `figure` and `load_benchmark` | CapScan metagenome benchmark extended to 295 acquisitions (pending; OPEN_ITEMS.md) |
| `COHORT_ORDER` | CAMPI, healthy stool, PREDICT, CapScan | producer x order | none |
| `BENCH_TEST_ARM`, `BENCH_REF_ARM` | FL_public vs MG | public-sequence FastaLake against matched metagenome | FL_cohortMG rows in the same tables (Supplementary Fig. 8) |
| gain (d; `stat_bench`) | 100 x (median over acquisitions of FL/MG peptide-count ratio - 1) | producer definition of the text gains (9.4%, 42.4%) | ratio of medians; not reported |
| recovery (e; `stat_bench`) | median over acquisitions of shared / MG peptides | producer definition (87.7%, 86.4%) | none |
| `JITTER_SEED_BC` | 20260927 (+0 for b, +1 for c) | producer jitter of points | none (display only) |
| `JITTER_SEED_DE` | 4 | producer jitter of points | none (display only) |
| `FIG_W_MM`, `FIG_H_MM`, `PNG_DPI` | 180 x 186 mm, 300 dpi | v30 placed size; PNG rasterised from the PDF | none |
| `COL` | all-prediction grey, same-half FastaLake blue, disjoint green; d,e MG orange, FL blue | v30 colours by meaning | none |

## Known limitations

- Panels d,e (and the CapScan same-spectrum numbers of line 65) will change when the
  CapScan metagenome benchmark is extended to 295 acquisitions; the inputs are parameters
  (`BENCH_TABLES`, `BENCH_N`) so the figure can be redrawn from new tables.
- All values are descriptive acquisition-level medians and interquartile ranges.
  Acquisitions are nested in participants (PREDICT, CapScan), and no inference across
  acquisitions is made in this figure.
- The participant and specimen columns of the split table were removed (raw PREDICT
  identifiers). Only their counts per cohort are kept (`split_half_structure.tsv`).
- Counts per search use each search's own native q values; the count ratios do not
  measure whole-workflow FDR or identical-peptide retention (the latter is
  `fraction_self_identities_retained`).
