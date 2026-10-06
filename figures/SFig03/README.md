# Supplementary Figure 3

CAMPI member support and identity recovery. The CAMPI SIHUMIx mock community has eight
known members and 15 public acquisitions. The figure shows how many peptides specific to
each member the FastaLake (FL) search accepted (a, b), how much of the reference-database
search each Sage arm recovers in the same acquisition (c), and the same recovery inside a
second search engine, MSFragger with Percolator (d). It supports the CAMPI statements on
v28d lines 83 and 89 and the legend on line 417. This unit builds the v30 figure: panels c,d
draw the nested 9-50-residue rows (the v28f figure drew 5-50), and the figure is redrawn in
the house style (Söhne, ink, arm colours as in main Figure 2, titles moved into the x-axis
labels). The v30 text for line 89 reads "In CAMPI, FastaLake yielded as many peptides as the
reference database (median ratio 1.00) and 0.94 times as many as MG, while recovering a median
85.3% of comparator peptide identities (Supplementary Fig. 3)".

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | pooled FL member-specific peptides per member, 5-50 residues, with acquisitions/15 | all eight members recovered; L. plantarum twelve pooled peptides (83) | 15 acquisitions; distinct peptides pooled over acquisitions | none (counts) | Panel a; v30 house style | CAMPI_member_support_recount.tsv (recounted from CAMPI_membership_*.tsv) | A14_campi_member_gap |
| b | as a, nested 9-50 residues | L. plantarum seven at 9-50; MG six of eight members at 9-50 (83) | 15 acquisitions; pooled distinct peptides | none (counts) | Panel b; v30 house style | CAMPI_member_support_recount.tsv | A14_campi_member_gap |
| c | Sage recovery of REF identities per acquisition, five arms, 9-50 | legend (417; v30 adds "c,d use peptides of 9–50 residues") | 15 acquisitions per arm | none (median bars) | Panel c; v30 house style | CAMPI_reference_recovery.tsv (recounted from CAMPI_membership_*.tsv) | A14_campi_member_gap |
| d | MSFragger recovery of REF and MG identities by FL, 9-50 | median ratios 1.00 and 0.94; median 85.3% recovered (89, v30 wording; v28d printed 0.99, 0.94, 84.6% at 5-50) | 15 acquisitions per comparator | none (median bars) | Panel d; v30 house style | CAMPI_within_MSFragger.tsv (checked against CAMPI_MSFragger_search_receipts.tsv) | campi_second_engine |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `MIN_LENGTHS` | (5, 9) | accepted 5-50 population and its nested 9-50 subset (Methods, line 685) | both lengths in all source tables; c and d draw 9 only |
| `MAX_LENGTH` | 50 | upper peptide length (Methods, line 685) | none |
| `PLOTTED_MIN_LENGTH_CD` | 9 | c and d draw the nested 9-50 population (len950 rebuild, as Comparisons in Methods use 9-50) | 5-50 rows (the v28f figure) in CAMPI_reference_recovery.tsv and CAMPI_within_MSFragger.tsv |
| `SPECIFIC_N_MEMBERS` | 1 | member-specific means the peptide occurs in exactly one of the eight reference proteomes (legend, line 417) | per-acquisition shared counts are not drawn; n_members column in CAMPI_membership_*.tsv |
| `N_ACQUISITIONS` | 15 | all CAMPI acquisitions (Methods, line 253) | none |
| `ARMS` | FL, FL_R, FL_R_M, REF, MG | the five Sage arms of the legend (line 417) | all five in every table |
| `MSF_COMPARATORS` | REF, MG | comparators of the MSFragger searches (line 89) | Sage versus MSFragger identities in second_engine_between_engine_identities.tsv (not drawn) |
| `ARMLAB` | FL; FL + Reference; FL + Reference + MG; Reference; MG | label edit FAR15 of the producer | none |
| `COL` | #2878a3, #df8f44, #619b68, #b35872, #8277ac | producer palette (producer step only) | v30 colours below |
| v30 arm colours (`build_figure_v30`) | FL `fl.FL` #0072B2, FL + Reference #009E73, FL + Reference + MG #CC79A7, Reference #7f7f7f, MG #E69F00; a, b and both d columns FL blue | colours by meaning as in main Figure 2 | producer palette in `COL` |
| `V30_PLACED_W_MM`, `V30_PLACED_H_MM`, `V30_PNG_DPI` | 160, 140, 300 | v30 figure drawn at placed size and rasterised at 300 dpi | none |
| `PLACED_W_MM`, `PLACED_H_MM` | 160, 140 | placed size in the manuscript | none |
| `BASE_FONT_PT`, `MIN_TEXT_PT` | 7.0, 5.0 | producer base size; smallest text at placed size is asserted >= 5 pt | none |
| `SLANT_SHEAR`, `N_SLANTED_LABELS` | 0.20, 16 | synthetic slant of the 16 taxon labels (label edit FAR16; Söhne has no italic) | none |
| `RENDER_DPI` | 400 | resolution of the producer-step display render (the v28f image) | producer PDF compared at 300 dpi; v30 PNG at 300 dpi |
| `XLIM_AB`, `YLIM_CD` | (0.8, 1e7), (0, 105) | producer axis limits | none |
| `CROP_BOXES` | see panels.py | per-panel cut-outs of the producer render, for display only | composites checked as a whole |

Medians are `numpy.median` over acquisitions; no random draw is made anywhere in this unit.
Upstream, Percolator ran with `--seed 913` per acquisition (campi_second_engine).

## Known limitations

- Member specificity is relative to the eight supplied reference proteomes, which contain
  surrogate strains; it is not universal species specificity.
- The tables behind a to c are recounted here from per-acquisition peptide membership
  tables. Panel d's table is checked against the accepted counts each search wrote, but
  the accepted MSFragger peptide lists themselves are not in this repository
  (OPEN_ITEMS.md, ZENODO).
- The v30 sentence for line 89 says "a median 85.3% of comparator peptide identities"
  without naming the comparator. At 9-50 residues the REF and MG medians (85.331% and
  85.324%) and the median over both (85.327%) all print as 85.3%.
- The "3- to 17-fold" range on line 83 holds at 5-50 residues (3.37 to 16.75). At 9-50 the
  range would be 3.38 to 22.6. The text does not name the length.
