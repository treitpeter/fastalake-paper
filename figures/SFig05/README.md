# Supplementary Figure 5

Supplementary Figure 5 shows what each error control in the paper measures, and why
they are not pooled into one error rate. Panels a-c give the fusion-linked diagnostic of
the linked target/control fusion assay in four cohorts, panel d the explicit-peptide
entrapment controls with a database fixed before control construction, and panels e,f
directly counted false discoveries in known-truth simulations. The continuation page
splits the CAMPI fusion diagnostic into its paired and unresolved parts.

Both pages are regenerated from `source_data/SFig05/` with the v30 plotting code (the
v28f/v29 producers restyled in house style without changing any plotted value) and are
pixel-identical to the v30 renders (receipts in `checks/`). Every printed number covered
here is recomputed (`NUMBERS.tsv`: 61 MATCH, 1 OPEN).

## Panels

| panel | what it shows | claim (v30 paragraph) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | fusion-linked diagnostic, 5-50 residues, q≤0.01, paired and unresolved parts, one bar per cohort and seed | legend 227; 368 | acquisitions summed within cohort and seed: CAMPI 15, healthy stool 58, PREDICT 306, CapScan 331; 2 seeds | none (descriptive ratio) | Panels a-f | `fusion_pooled_current.tsv` (drawn); `fusion_acquisition_components.tsv` (recomputation) | healthy_stool_donor_entrapment; campi_fusion_375m |
| b | as a, 9-50 residues | legend 227 | as a | none | Panels a-f | as a | as a |
| c | diagnostic of the nested 5-8-residue population, one dot per seed, bar at the median | legend 227 | as a | none | Panels a-f | as a | as a |
| d | explicit-peptide paired estimate, self-selected and donor-selected, 5-50 and 9-50 residues | 30 (1.0%, 0.49%); legend 227; 373-374 | acquisitions: 5 self-selected, 58 donor-selected; each the mean of 3 seeds | none (median) | Panels a-f | `donor_entrapment_seed_means.tsv` (drawn); `donor_recount.tsv` (recomputation) | healthy_stool_donor_entrapment |
| e | realised FDP per mixture, baseline regime, three arms; means with 95% bootstrap intervals | 46; legend 227; 376 | synthetic mixtures, 12 | percentile bootstrap interval of the mean | Panels a-f | `simulation_mixture_outcomes.tsv`, `simulation_uncertainty.tsv`, `simulation_all_outcomes.tsv` | restricted_reference_simulation |
| f | as e, stress regime | 46; legend 227; 376 | synthetic mixtures, 12 | percentile bootstrap interval of the mean and of the paired self minus independent difference | Panels a-f | as e | restricted_reference_simulation |
| continued | CAMPI fusion diagnostic per seed and length window, paired (E + E₀ + 2B)/N and unresolved U/N | 45; legend 229; 367 | 15 acquisitions summed within seed; 2 seeds | none (descriptive ratio) | Panel continued | `fusion_decomposition_current375.tsv` (drawn); `fusion_acquisition_components.tsv` (recomputation) | campi_fusion_375m |

Panels a-f are one figure (160 x 189.3 mm; the continued page 180 x 90 mm); they are regenerated together by
`panels.figure_a_to_f`, which calls `panel_ab`, `panel_c`, `panel_d` and `panel_ef`.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `Q_CUTOFF` | 0.01 | nominal peptide q of all panels (legend 227, 229; Methods 367, 373, 376) | q 0.001, 0.005, 0.02, 0.05 in `fusion_pooled_current.tsv`, `donor_recount.tsv` |
| `FUSION_WINDOWS_AB` | 5 and 9 (to 50) | reporting populations of a,b (legend 227) | both drawn |
| `FUSION_MAX_RESIDUES` | 50 | upper search length (Methods 367) | none |
| `SHORT_WINDOW` | 5-8 | nested short population, same search q values (legend 227, Methods 367) | none |
| `FUSION_SEEDS` | 20260907, 20260908 | two positional-permutation control seeds (line 45, "two randomizations") | each seed shown separately |
| `DONOR_TIERS` | T0_self, T1_within_group | self-selected and donor-selected strata (legend 227, Methods 374) | the strata are kept separate |
| `DONOR_SEEDS_PER_ACQUISITION` | 3 | seeds averaged within target acquisition (Methods 374) | per-seed rows in `donor_recount.tsv` |
| `SIM_ARMS` | fixed, independent, self | the three ordinary searches per mixture (Methods 376) | fusion and appended-entrapment arms in `simulation_all_outcomes.tsv` |
| `SIM_REGIMES` | baseline, stress | Methods 376 | none |
| `BOOTSTRAP_DRAWS` | 20,000 | Methods 376 | t intervals also computed, in `simulation_uncertainty.tsv` |
| `BOOTSTRAP_SEED` | 2026091051 | seed of the producing collect step; not stated in the text | none |
| `BOOTSTRAP_LEVEL` | 0.95 | Methods 376 | none |
| `PAIRED_COLOR`, `UNRESOLVED_COLOR`, `DOT_COLOR` | FastaLake blue, light grey, FastaLake blue | v30 colours by meaning (orange is reserved for MG in the main figures) | v28f colours (blue/orange/green/purple) |
| `AF_SIZE_MM`, `CONT_SIZE_MM`, `LETTER_X_MM`, `PNG_DPI` | 160 x 189.3 mm, 180 x 90 mm, 1 and 82 mm, 300 dpi | v30 placed sizes, one letter x per column, PNG from the PDF | none |
| `NOMINAL_LINE_PCT` | 1 | dotted reference line (legend 227, 229); not a calibration claim | none |

## Known limitations

- The fusion controls of panels a-c for CapScan use 331 acquisitions, not the v28f roster
  of 347 (see `OPEN_ITEMS.md`).
- Panels a-c compare fusion controls built on different reservoirs (CAMPI broad UniProt,
  375,333,654 sequences; the other cohorts their own smaller reservoirs). The diagnostic
  is not a calibrated false-discovery estimate (Methods 368).
- The simulation counts errors only in a restricted synthetic setting (Methods 376).
- The 375,333,654 reservoir size printed in the continuation legend is not recomputed by
  this unit.
- Upstream analyses (fusion searches, donor searches, simulation) are not rerun here; this
  unit starts from their saved component tables.
