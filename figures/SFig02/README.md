# Supplementary Figure 2

Supplementary Figure 2 tests the two choices that make peptide-guided database selection work.
Panels a,b measure how often a peptide query matches a reservoir protein by chance, using
composition-preserving scrambles of the predictions of one healthy-stool acquisition. Panels c,d
show how sensitive the selected databases and the accepted peptides are to the minimum lookup
length and to the top-scoring fraction. The continuation (panels e,f) reports an exploratory
comparison of bag-aware and mass-substitution lookup against exact lookup, and a separate fragment
audit. Production analyses used exact lookup.

Panels a-d are one canvas in the original producer, and so are e,f. Each canvas is regenerated
and pixel-checked as a whole. The unit also builds the v30 figure: both canvases redrawn in the
house style (Söhne, ink, Okabe-Ito colours by meaning, placed widths 160 and 155 mm, subpanel
titles removed); `build_v30` proves every plotted value identical to the producers, and both v30
images are pixel-identical to the v30 PNGs (`checks/SFig02_v30.json`,
`checks/SFig02_continued_v30.json`). All four checks are `IDENTICAL` (0 differing pixels), and the two
regenerated images placed in the manuscript are byte-identical to v28d media `image8.png` and
`image9.png`.

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | fraction of distinct queries matching at least one reservoir protein, by length 9-50, authentic and five scrambled seeds; inset lengths 9-13 | scrambled matches 8.1-8.9% at 9, 0.7-1.1% at 10, none at >=13 (59; legend 409) | 1 acquisition; sequence (distinct query) per arm, seed and length; 24,560 authentic queries | none (descriptive fractions) | Panels a-d; v30 house style | `selection_length_lookup.tsv`, `lookup_query_outcomes.tsv` | database_selection_specificity |
| b | distinct authentic queries per length (denominator of a) | legend 409 | 1 acquisition; sequence | none | Panels a-d; v30 house style | `selection_length_lookup.tsv` | database_selection_specificity |
| c | median accepted peptides (exclusive quartiles) and median selected proteins for minimum lookup lengths 6-12 at top 30% | 679; legend 409 | 58 acquisitions per setting, same acquisitions in all settings; acquisition | none (median, exclusive quartiles) | Panels a-d; v30 house style | `selection_grid_summary.tsv`, `selection_grid_acquisitions.tsv` | healthy_stool_selection_sweep |
| d | as c for top-scoring fractions 5-100% at minimum length nine | 679; legend 409 | as c | none | Panels a-d; v30 house style | as c | healthy_stool_selection_sweep |
| e | gained, lost and net peptide-acquisition identifications of six relaxed-lookup arms against exact lookup | 67, 561; legend 413 | 48 acquisitions (15 CAMPI, 32 healthy stool, 1 CapScan); peptide (peptide-acquisition identification) | none (descriptive sums) | Panels e,f; v30 house style | `bags_bounded48_per_acquisition.tsv` | bag_mass_lookup_search |
| f | gained composition-preserving rearrangement pairs with and without a diagnostic fragment ion | 67, 581; legend 413 | 215 same-spectrum pairs in 6 acquisitions; peptide | none (descriptive counts) | Panels e,f; v30 house style | `bags_fragment_audit_summary.tsv` | bag_mass_lookup_search |

The 93.0% coverage of the eight-ordering bound (v28d line 561) is computed from
`bags_bounded48_orderings.tsv` by `panels.py:stat_e_bags`. Design settings (seeds, lengths, score
fractions, q threshold) are in the Choices table, not in `NUMBERS.tsv`.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `SCRAMBLE_SEEDS` | 20260901, 20260913, 20260914, 20260915, 20260916 | fixed scramble realizations; 20260901 is the historical scramble reproduced byte for byte (v28d line 665) | terminal-preserving scrambles with the same seeds, v28d line 667 (not in this figure) |
| `MIN_LOOKUP_LENGTH`, `MAX_LOOKUP_LENGTH` | 9, 50 | production query length window (v28d line 557) | minimum lengths 6-12 in panel c |
| `TOP_SCORE_FRACTION` | 0.30 | production lookup threshold (v28d line 557) | 5-100% in panel d |
| `LENGTH_GRID` | 6-12 | sensitivity grid (v28d line 679) | none |
| `SCORE_GRID` | 0.05-1.00, 12 values | sensitivity grid (v28d line 679) | none |
| `QUARTILE_METHOD` | exclusive (`statistics.quantiles`, n=4) | as used by the upstream sweep; legend line 409 says exclusive-quartile ranges | none |
| `PEPTIDE_Q` | 0.01 | nominal peptide q threshold for gains and losses (v28d lines 413, 579) | none |
| `ORDERING_BOUNDS` | 2, 4, 8 | bag-ordering bounds of the e arms (v28d lines 413, 579) | unrestricted bags in the same panel |
| `FRAGMENT_TOLERANCE_PPM`, `FRAGMENT_CHARGES` | 20, (1, 2) | fragment-audit definition of a diagnostic ion (v28d line 581); applied upstream, documented here only | none |
| `SHORT_LENGTH_CUTOFF` | 13 | length from which no scrambled query matched (v28d line 59) | full curve in panel a |
| `EF_PLACED_WIDTH_IN`, `EF_TARGET_DPI` | 6.1 in, 400 dpi | re-export of the e,f vector figure at its placed width | none |
| `PDF_RASTER_DPI` | 300 | rasterisation of the a-d vector figure for the manuscript | none |
| `V30_W_ABCD`, `V30_H_ABCD`, `V30_W_EF`, `V30_H_EF` | 160.0, 135.2, 155.0, 117.1 | v30 canvases drawn at their placed sizes (mm) | none |
| `V30_MARKER_PT` | 3.2 | v30 marker size at placed size (producer used 6 on a larger canvas) | none |
| `V30_PNG_DPI` | 300 | v30 PNGs rasterised from the PDFs | none |

## Known limitations

- Panels a,b come from one representative healthy-stool acquisition and five scramble
  realizations. The ranges are ranges of five realizations, not confidence intervals.
- In panels c,d the same 58 acquisitions enter every setting; quartiles describe spread between
  acquisitions, and no test of differences between settings is made. The 19 arms include the
  older 97-pool reference arm at nine residues and 30%, so the grid has 18 distinct settings.
  The reference arm is in the source data but not drawn.
- Peptide counts in c,d strip modifications and keep isoleucine and leucine distinct, so they are
  not directly comparable with the primary peptide counts (v28d line 679).
- Panel e counts are summed peptide-acquisition identifications at nominal q, not distinct
  sequences, and the relaxed-lookup error rate is not established. Panel f compares two candidate
  sequences per spectrum; neither is known truth. The two panels use different acquisitions.
- The bar lengths of panel f (193, 22) are literals in the original plotting code. The original
  guarded them with an assertion against the audit table; that assertion is kept in
  `panels.py:panel_ef`.
- The original a-d and e,f producers used matplotlib's default typeface (DejaVu Sans), not the
  house Söhne style. `panels.py` therefore sets the producers' own rcParams and does not apply
  `fl_style` to these panels.
- `lookup_query_outcomes.tsv` (18.8 MB) is below the 20 MB limit and is kept in the repository so
  that the per-length counts of a,b are recounted from individual query outcomes.
