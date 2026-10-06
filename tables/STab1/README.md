# Supplementary Table 1

Supplementary Table 1 reports how clustering candidate proteins with MMseqs2/Linclust
before FastaLake parsimony changes the search result in 15 CAMPI acquisitions. Each row is
one acquisition and one arm: the unclustered search and six identity thresholds (99, 97,
95, 90, 75 and 50%), 105 rows and 21 columns. Counts are accepted peptides, searched
protein entries and compatible KEGG Orthology (KO) identifiers; ratios, losses, gains and
retained fractions are relative to the unclustered search of the same acquisition. The
notebook rebuilds the table from `source_data/STab1/`, checks it against the table as
produced (byte-identical TSV, cell-identical XLSX), repeats the producer's consistency
checks and prints the medians quoted in the Results.

The legend text for v28f is `expected/LEGEND.txt` (the v28d legend at line 475 is the
title only; see OPEN_ITEMS.md).

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| table | 105 acquisition-arm rows, 21 columns | "all 105 acquisition–arm configurations and paired counts are supplied in Supplementary Table 1" (695); legend title (475) | 15 acquisitions x 7 arms | none | Table | per_acquisition_identifications.tsv; per_acquisition_KO_coverage.tsv | campi_candidate_clustering (not in pipelines/) |
| medians | median paired peptide reduction and KO recovery per threshold | "reduced accepted peptide counts in all 15 acquisitions ... by a median of 0.9%, 2.2%, 2.5%, 4.1%, 7.2% and 8.8% ... retaining 96.4–99.0% of annotated KEGG functions" (89) | 15 acquisitions per threshold, paired with own unclustered search | none (descriptive median) | Checks against the table as produced and recomputed numbers | per_acquisition_identifications.tsv; per_acquisition_KO_coverage.tsv | campi_candidate_clustering |
| checks | 45 earlier rows, 60 new rows and an independent recount agree with the table | none printed | 105 rows | none | Producer checks | support_97_99_identification_comparison.tsv; check_new_rows_95_90_75_50.tsv; independent_recount_*.tsv | independent recount |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| ARM_ID | unclustered, 99, 97, 95, 90, 75, 50% (Linclust --min-seq-id) | thresholds of the sensitivity analysis (Methods, line 695) | 97/99% only in the earlier 45-row table (support_97_99_*) |
| LINCLUST_COVERAGE | 0.8 | Methods, line 695 | none run |
| LINCLUST_COV_MODE | 0 | Methods, line 695 | none run |
| LINCLUST_CLUSTER_MODE | 2 | Methods, line 695 | none run |
| PEPTIDE_Q_CUTOFF | 0.01 (peptide q) | Methods, line 695; LEGEND.txt | none run |
| PEPTIDE_LENGTH | 5 to 50 residues, I and L merged | LEGEND.txt; upstream peptide parser | none run |
| N_ACQUISITIONS | 15 | all CAMPI acquisitions (legend, line 475) | none |
| N_ROWS | 105 | 15 acquisitions x 7 arms (Methods, line 695) | none |
| FLOAT_TOL | 1e-12 | producer tolerance for fraction cross-checks | none |

The Linclust and peptide settings were applied upstream; in `panels.py` they are recorded
constants, not inputs to any computation in this repository.

## Known limitations

- The searches, clustering and KO mapping that produced the two `per_acquisition_*`
  tables are not re-run here and are not yet in `pipelines/` (OPEN_ITEMS.md, PIPELINE).
- By their names, `S11_Fraction1` to `S11_Fraction4` are four fractions of one CAMPI
  sample, so the 15 acquisitions are not 15 independent specimens. The medians are
  descriptive; no test is made.
- KO recovery uses one fixed peptide-to-KO mapping; it counts annotated function
  evidence, not function activity.
- The 50% threshold is a stress test of candidate clustering after reservoir lookup; it
  does not cluster the reservoir and does not calibrate whole-workflow FDR.
