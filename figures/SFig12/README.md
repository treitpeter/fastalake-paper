# Supplementary Figure 12

Shows how much of the measured functional signal can be given a KEGG ortholog (KO) and, within
that, how much can be given a species label and how stable the leading species is. Panel a is the
KO-allocated fraction of accepted feature intensity per acquisition in three cohorts. Panels b and e
are healthy stool: the expected number of distinct leading species per function after rarefaction to
eight acquisitions, and leader agreement under two resolution filters. Panels c, d and f are CapScan:
per participant, the leading species' share of resolved signal (c), the resolved share of non-host
signal (d), and agreement of the leading species between intestinal regions and own stool under four
eligibility rules (f). In v30 all panels use direct peptide assignment; the original
protein-grouping series was removed from e and f.

In this unit the column `v28d_line` of NUMBERS.tsv and STATS.tsv holds the paragraph index of the
v30 text (`visible_v30f.txt`, first column), not a v28d line. Every quote was checked against that
text with `tools/check_numbers.py --manuscript-text` (one paragraph per line).

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | KO-allocated fraction of accepted feature intensity per acquisition, by cohort; black bar = median | legend (v30 paragraph 243) | 710 acquisitions: 58 HSTOOL, 306 PREDICT, 346 CapScan | none (descriptive) | Panel a | `source_data/SFig12/{HSTOOL,PREDICT,CAPSCAN}_KO_all.tsv`, `source_data/SFig12/KO_allocation_coverage.tsv` | peptide_centric_functional_allocation (not in `pipelines/`) |
| b | expected distinct leading species in eight resolved acquisitions, per function, by healthy-stool design | legend (243) | 3,862 function x design rows; acquisitions per design 10-15 | none (analytical rarefaction expectation) | Panel b | `source_data/SFig12/healthy_stool_leader_diversity.tsv` | peptide_centric_functional_allocation |
| c | per-participant ECDF of leader share within species-resolved signal | legend (243) | 21,168 function x participant cells, 1,467 functions, 15 participants | none (descriptive) | Panel c | `source_data/Fig5/contributor_cells.tsv` | peptide_centric_functional_allocation |
| d | per-participant ECDF of species-resolved share of non-host functional signal | median 29.1% (82); counts (243, 402) | as c | none (median over cells) | Panel d | `source_data/Fig5/contributor_cells.tsv` | peptide_centric_functional_allocation |
| e | healthy-stool equal-function leader agreement, direct peptide assignment, by design (technical, workflow) and filter (primary, well resolved) | legend (243) | 4 design x filter points; 694-768 functions per point | none (descriptive) | Panel e | `source_data/SFig12/healthy_stool_grouping_leader_sensitivity.tsv` (rows `method == simplified`, `subset == all_eligible`) | peptide_centric_functional_allocation |
| f | CapScan donor-balanced leader agreement, region vs own stool, under four eligibility rules, direct peptide assignment | legend (243); 46.4% and 57.6% cited with Fig. 5f (81) | 15 participants, 60 participant-region combinations | none (descriptive) | Panel f | `source_data/SFig12/CapScan_grouping_leader_sensitivity.tsv` (plotted, rows `A33_CAPSCAN__simplified`), `source_data/SFig12/capscan_contributor_region_cells_direct.tsv` and `source_data/Fig5/contributor_cells.tsv` (recomputed) | peptide_centric_functional_allocation |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `RAREFACTION_DEPTH` | 8 | legend, v30 paragraph 243; read from the source-data column name; the expectation is exact (hypergeometric), no draws | none |
| `ELIGIBILITY_FRACTION` | 0.8 | healthy stool: resolved leader in at least 80% of a design's acquisitions (Methods, v28d lines 249, 783) | panel b and e use per-design eligibility; functions eligible in every design are flagged (`eligible_in_all_five_designs`, 582 per design); see OPEN_ITEMS.md |
| `MIN_PARTICIPANTS_PER_KO` | 12 | CapScan: resolved leader in at least twelve of fifteen participants (Methods, v28d lines 249, 785) | none |
| `STRICT_RESOLVED_MIN` | 0.2 | strict rules: at least 20% species-resolved signal (v28d lines 249, 785) | the three strict rules in panel f |
| `STRICT_GAP_MIN` | 0.05 | strict rules: leader margin of at least 5 percentage points (v28d lines 249, 785) | the three strict rules in panel f |
| `COHORTS` | HSTOOL, PREDICT, CAPSCAN | producer x-axis order in panel a | none |
| `FIG_MM`, `FONT_PT`, `PNG_DPI` | 180 x 172 mm, 6 pt, 300 dpi | v30 build: authored at the placed size in fl_style (Söhne, Okabe-Ito) | none |

## Known limitations

- Panel b plots values computed upstream (A33 `derive_hstool_pc.py`). The per-acquisition leader
  lists needed to recompute the rarefaction are not in this repository.
- Panel e plots values computed upstream (A33 S10 summary); the pairwise inputs are not in this
  repository.
- Panel f is plotted from the producer's table and its four direct-assignment rows are recomputed
  here from the cell tables with the producer's `four_rules`; the two agree to 5.6e-17
  (floating-point summation order), and all counts are equal.
- The two sensitivity tables still contain the original protein-grouping rows, which are Source
  Data but are no longer drawn (OPEN_ITEMS.md).
- The `A33_` prefix in the `source` column of the panel f table is a producer label; those rows
  come from the 347-roster rerun.
