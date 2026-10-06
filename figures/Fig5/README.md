# Figure 5

Figure 5 shows that each CapScan participant keeps a characteristic microbial functional
profile across four intestinal sampling regions and stool, while the species that lead
individual functions often change between regions. It uses 347 LC-MS/MS acquisitions from
308 specimens of 15 participants; the cross-omics panel uses the 13 participants with DNA,
protein and metabolite profiles in all five sampling groups. The figure is unchanged from
v28d to v28f.

`Fig5.ipynb` redraws panels b-g with the original drawing code, composes them with the
panel-a artwork, checks every render pixel by pixel against the original producer's output
(all IDENTICAL with the Söhne fonts, see `checks/`), and recomputes every Figure 5 statistic
printed in the text (`NUMBERS.tsv`, `STATS.tsv`).

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | study design: four capsule regions and stool, acquisition counts, data layers | 347 acquisitions, 308 specimens, fifteen participants (147) | 347 acquisitions, 15 participants | artwork | artwork | artwork (numbers in `ARTWORK_NUMBERS.md`, checked against `design_347.tsv`) | capscan_participant_roster |
| b | human-assigned, microbial and unresolved shares of LFQ intensity; participant points; q versus stool | regions 1-4 exceed stool by 12.1, 11.1, 8.4 and 4.7 percentage points (149) | 15 participants (75 participant-region cells, 346 acquisitions) | exact paired sign-flip, bootstrap CI, BH over 4 | Panel b | `origin_fraction_cells.tsv`, `four_region_tests.tsv`, `taxonomy_signal.tsv`, `design_347.tsv` | capscan_functional_profiles |
| c | mean distance of capsule profiles to own versus other stools | all fifteen closer to own stool, mean difference 0.061 (155) | 15 participants, 75 profiles, 73,316 OGs | stool-label permutation (99,999), bootstrap CI | Panel c | `own_stool_participants.tsv`, `similarity_summary.tsv`; Zenodo `OG_other_unresolved.npz` | capscan_functional_profiles |
| d | variance partition with and without stool | participant 45.0%, sampling group 14.8%; capsules only 56.9% and 5.3% (155) | 75 and 60 participant-region profiles | descriptive | Panel d | `orthogonal_partition.tsv`; Zenodo `OG_other_unresolved.npz` | capscan_functional_profiles |
| e | leading annotated species for 24 functions, participant x region | 726 of 1,440 comparisons; 111 of 360 profiles (159) | 1,440 function-participant-region comparisons | descriptive | Panel e | `mosaic_donor_regions.tsv`, `mosaic_abundance.tsv`, `heatmap_function_key.tsv` | peptide_centric_downstream |
| f | leading-species agreement under direct peptide assignment, all eligible functions and strict subset (v30; production only) | 46.4% and 57.6% (157; v30 Results paragraph 81) | 15 participants | descriptive | Panel f | `participant_leaders_resolved.tsv`, `participant_leaders_strict_both_methods.tsv` | peptide_centric_downstream |
| g | Spearman correlations of region-adjusted distances, native, OG and KO | rho 0.51, 0.55, 0.38, all q<0.001 (161) | 13 participants, 65 profiles | Spearman, block permutation (1,999), BH over 6 native and 16 functional | Panel g | `native_correlations.tsv`, `functional_correlations.tsv`, `native_distances_long.tsv`, `functional_distances_long.tsv` | capscan_crossomics |
| text | leading-species share and Pielou evenness, fluid versus stool | 22.5% vs 12.7%; 0.623 vs 0.669; P=6.1 × 10⁻⁴ and 3.1 × 10⁻⁴ (159) | 15 participants | exact Wilcoxon signed-rank | Text statistic: leading-species share and evenness | `evenness_per_acquisition.tsv`; Zenodo `function_species.npz` | capscan_dominance_evenness |
| text | species-label coverage of non-host functional intensity | median 29.1% (159) | 21,168 function-participant cells | descriptive | Panel e | `contributor_cells.tsv` | peptide_centric_downstream |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `NP` | 15 | participants with acquisitions in all five sampling groups | 14-participant roster of earlier versions (superseded, see PROVENANCE.md) |
| `N_DNA_PARTICIPANTS` | 13 | participants with complete DNA, protein and metabolite profiles | 15-participant non-DNA family: Supplementary Fig. 11 and Source Data (Methods line 837) |
| `BOOT_SEED` | 20261003 | seed of the participant bootstrap (code; not stated in the text) | none |
| `BOOT_DRAWS` | 2000 | participant bootstrap resamples (Methods line 843) | none |
| `OWN_STOOL_SEED` | 20260910 | permutation seed (line 849) | none |
| `OWN_STOOL_PERMUTATIONS` | 99999 | stool-label permutations (line 849) | none |
| `PREVALENCES` | 0.0, 0.1, 0.2 | all OGs (primary) and two prevalence filters (line 847) | 10% and 20% results in `similarity_summary.tsv` |
| `CROSSOMICS_SEED` | 20260907 | block-permutation seed (code; not stated in the text) | none |
| `CROSSOMICS_PERMUTATIONS` | 1999 | whole-participant block permutations (line 835) | none |
| `FUNCTION_SCOPE` | matched12 | the matched 13-participant functional family (the label predates participant 15) | `displayed_mixed12_14` (15/13 participants) in `functional_correlations.tsv` |
| `NATIVE_ANALYSIS` | sensitivity_all_features | native profiles with all features, as in Fig. 5g | `primary_prevalence10` in `native_correlations.tsv` |
| `EVEN_BAD_LABELS` | Host, Mixed host/microbial, Conflicting species, No species annotation | excluded from named-species evenness (line 849) | none |
| `PANEL_F_MEAN_LABELS` | 46.4, 57.6 | asserted panel f mean labels, all eligible and strict subset under direct assignment (v30) | original protein grouping (42.5%, 58.3%) removed from the manuscript in v30 |
| `JITTER_SEED` | 5 | horizontal jitter of panel b points, display only | none |
| `PAGE_W_PT`, `PAGE_H_PT`, `Y0` | 518.74, 688.28, 140.0 pt | page geometry of the producer (183 x 242.8 mm) | none |

## Known limitations

- Panels b-g are drawn on one shared page, as by the original producer. The per-panel PNGs
  (`build/Fig5/panel_<x>.png`) are crops of that page written after all panels are drawn,
  with the boxes in `panels.PANEL_BOXES_PT`; neighbouring panels can show at crop edges.
- Panels c and d and the per-acquisition evenness metrics are recomputed only when the
  three Zenodo files (`source_data/Fig5/README.md`) are in `build/zenodo/Fig5/`. Without them
  the notebook still draws every panel and recomputes all other numbers.
- The tables behind panels e and f (`mosaic_donor_regions.tsv`, `participant_leaders_*.tsv`,
  `contributor_cells.tsv`) and the distance matrices behind panel g are outputs of the
  upstream analyses; their computation from peptide-level data is not repeated here.
- Pixel identity needs the Söhne TrueType files `soehne-buch.ttf` and `soehne-halbfett.ttf`
  in `FASTALAKE_FONT_DIR` (the files the producer used).
- Open items, including one Results sentence without a computation, are in `OPEN_ITEMS.md`.
