# Figure 4

Figure 4 compares public-sequence FastaLake (FL) databases with specimen-matched metagenome (MG) databases in PREDICT stool and in CapScan stool and intestinal fluid, explains why FL accepted spectra the MG search did not, and shows PREDICT molecular profiles: the effect of antibiotic exposure on the human-assigned signal and protein-group depth, host, microbial and metabolite features across the Enterotype Dysbiosis Score (EDS), and the agreement of patient-profile distances across DNA, RNA, protein and metabolite layers. This folder regenerates the v28f version of the figure, in which panels b,c show the 10-fold cross-fitted comparison at 9-50 residues (the former same-spectrum panels move to Supplementary Fig. 6d,e).

Run `python paper.py build Fig4` from the repository root. Every quantitative panel is redrawn from `source_data/Fig4/` by the original producer's code (`panels.py`) and is pixel-identical to the producer's own output at 300 dpi with the Söhne fonts (`checks/`; e, g and h checked outside the repository, OPEN_ITEMS). Panels e, g and h need controlled clinical inputs (`FASTALAKE_CONTROLLED_DIR`); without them those panels are skipped and the public aggregate statistics are printed instead.

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | PREDICT study design and matched measurements | 306 acquisitions with companion DNA, RNA, metabolites (121) | 306 acquisitions; 296 matched | artwork | artwork | artwork (`ARTWORK_NUMBERS.md`) | none |
| b | peptides of 9-50 residues per acquisition, MG vs FL, 10-fold cross-fitted | cross-fitted gains −3.5%, −5.5%, 28.3% with CIs; regions 21-36% (123, 133) | 296 PREDICT (97 participant clusters), 56 CapScan stool (15), 226 intestinal fluid (14) acquisitions | median gain with participant-cluster bootstrap CI (no test); mixed model for line 129 (unadjusted Wald z and 3-df Wald; Methods line 229 describes LRT and Holm, OPEN_ITEMS) | Panel b | `crossfit_per_acquisition.tsv` | capscan_crossfit_10fold |
| c | MG-search peptide identities recovered by cross-fitted FL, 9-50 residues | 78.4%, 72.0%, 75.4% (v30 paragraph 65; not in v28d) | as b | none | Panel c | `crossfit_per_acquisition.tsv` | capscan_crossfit_10fold |
| d | MG-search outcome of PREDICT spectra accepted only by FL | 29% absent, about 69% present but top MG target failed (131) | 6,689,188 spectra from 296 acquisitions | none | Panel d | `spectrum_audit_per_acquisition.tsv` | predict_fl_only_spectrum_audit |
| e | human-assigned LFQ share off and on antibiotics | +11.1 pp (95% CI 5.6-17.2, q = 0.00068) (137) | 42 patients, off and on | two-sided paired sign-flip test; bootstrap CI; BH over 2 endpoints | Panel e | controlled `paired_antibiotic_profiles.tsv`; `clinical_paired_statistics.tsv` | predict_clinical_refit |
| f | mean origin composition of LFQ signal by antibiotic state | legend (127) | 42 patients (aggregate) | none | Panel f | `clinical_composition_means.tsv` | predict_clinical_refit |
| g | quantified protein groups off and on antibiotics | −2,214 (95% CI 998-3,469, q = 0.0013) (137) | 42 patients, off and on | as e | Panel g | controlled `paired_antibiotic_profiles.tsv`; `clinical_paired_statistics.tsv` | predict_clinical_refit |
| h | 19 host-protein, microbial-pathway and metabolite features across 85 patients ordered by EDS (top track labelled "Enterotype dysbiosis score (0–1)" in v30) | display (127, 137) | 85 patients | none | Panel h | `heatmap_features.tsv`, `heatmap_z.tsv`, `heatmap_tracks.tsv`; controlled `heatmap_85_patient_order.tsv` | predict_eds_screen, predict_eds_feature_selection |
| i | patient-profile distance correlations, native, eggNOG OG, KEGG KO | ρ = 0.44, 0.45, 0.44, all q < 0.001 (139) | 77 patients (2,926 pairs) | Mantel-type Spearman correlation, 1,999 permutations, BH | Panel i | `crossomics_native_correlations.tsv`, `crossomics_functional_correlations.tsv`, `crossomics_distances_predict77.tsv` | predict_crossomics_refit |

Text numbers of the Figure 4 paragraphs that no panel prints (lines 121-133: same-spectrum benchmark, gains by CapScan sampling site, mixed model, DNA amounts, database sizes) are recomputed in the section "Text numbers" from `same_spectrum_benchmark.tsv`, `crossfit_per_acquisition.tsv` and `capscan_dna_yield_by_acquisition.tsv`. Numbers in these paragraphs that cite another figure or table are covered by that unit: spectrum-split losses 11.4% and 17.4% (Supplementary Fig. 6), curated specimen MG 9.1%, 8.2%, 286 (Supplementary Fig. 8), tertile gains 80%, 69%, 50% (Supplementary Fig. 9), rifaximin sensitivity (Supplementary Table 3), metabolite-enzyme pairs (Supplementary Fig. 10, Supplementary Table 4), per-function correlations (Supplementary Table 5).

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `LENGTH_WINDOW_FIG4BC` (`L_BC`) | 9-50 residues | v28f Fig. 4b,c show the cross-fitted comparison at the length the text reports (line 123) | 5-50 residues: same-spectrum in Supplementary Fig. 6d,e; cross-fitted 5-50 in Source Data |
| `LENGTH_WINDOW_SAME_SPECTRUM` (`L_SAME`) | 5-50 residues | production search length window (v28d line 123); used for the v28d comparison values | v30 text quotes 9-50 (`TEXT_LENGTH`) |
| `BOOT_DRAWS_CROSSFIT`, `BOOT_SEED_CROSSFIT` | 10,000; 20261004 | participant cluster bootstrap of median gains (Methods lines 267, 777) | separate per-region stream gives stool upper bound −2.53 instead of −2.26 (line 777) |
| `CI_PERCENTILES` | 2.5, 97.5 | percentile 95% bootstrap intervals (lines 267, 275, 777) | none |
| `BOOT_SENSITIVITY_SEEDS` | 20261005-20261024 | repository Monte Carlo check of the stool CI upper bound (not in the manuscript) | result in OPEN_ITEMS and PROVENANCE.md |
| `CROSSFIT_STRATA` | PREDICT, CapScan all, stool, intestinal fluid | stream order of the original bootstrap; changing it changes the CIs | none |
| `MIXED_MODEL_LENGTH` | 9-50 residues | v30 Results paragraph 69 (1.40-fold, 95% CI 1.23-1.60, P = 4 × 10⁻⁷; regions P = 0.63) | 5-50: fold 1.42 (1.24-1.63), P = 2.9 × 10^-7, regions P = 0.61, as printed in v28d line 129 (printed by the notebook, Panel b cell) |
| `SAME_VS_CROSSFIT_LENGTH` | 9-50 residues | v30 paragraph 70: same-spectrum gains exceed cross-fitted gains in all 296 PREDICT and all 282 CapScan acquisitions | 5-50: 295 of 296 and 281 of 282 (v28d line 133; printed by the notebook) |
| `TEXT_LENGTH` | 9 (9-50 residues) | v30 paragraphs 64, 65, 70 print same-spectrum gains, recoveries and fluid/stool yields at 9-50 | 5-50 values (v28d lines 123, 129) printed by the notebook, Text numbers cell |
| `H_EDS_TRACK_LABEL` | Enterotype dysbiosis score (0–1) | v30 Fig. 4h track label (Peter, 6 Oct 2026) | v28f printed "EDS (0–1)" |
| `CLIN_BOOT_DRAWS`, `CLIN_SIGNFLIP_DRAWS`, `CLIN_SEED` | 50,000; 200,000; 20260907 | paired antibiotic inference of clinical_refit.py | rifaximin definitions: Supplementary Table 3 |
| `MANTEL_PERMUTATIONS`, `MANTEL_SEED` | 1,999; 20260907 | cross-omics distance correlations | prevalence-filtered profiles (≥10%) in Supplementary Fig. 11 and Source Data |
| `Q_STAR` | 0.05 | single-star threshold in Fig. 4i, as drawn by the producer | legend line 127 defines three tiers (OPEN_ITEMS AUTHOR) |
| `BOX_WHIS` | 1.5 IQR | Fig. 4e,g whiskers (legend line 127) | none |
| `HEATMAP_CLIP` | ±2 | display range of within-feature z scores in Fig. 4h | unclipped values in `heatmap_z.tsv` |
| `CAPSCAN_REGIONS` / `REGION_LABEL` | duodenum, jejunum, ileum, colon = regions 1-4 | device-defined sampling regions; source labels kept, manuscript names regions 1-4 | none |

## Known limitations

- Panels b and c are drawn by one producer call on one transparent page the size of the composite; they are checked together (`checks/panel_bc.json`).
- Panel a is artwork assembled over several steps (ASSEMBLY.md); the original SVG's author step was not located.
- Panels e, g and h need controlled inputs. Their producer renders and check receipts are not in `expected/` and `checks/` until Peter decides (OPEN_ITEMS AUTHOR); the stored notebook outputs come from a run without controlled inputs. The published composites (`Fig4_v28f.*`, `Fig4_v30.*`) and the vector panel-a artwork (`panel_a.pdf`, which carries the whole page) show panels e, g and h, so they are held with the controlled inputs (`$FASTALAKE_CONTROLLED_DIR/Fig4_expected/`, `.../Fig4_artwork/`, 6 Oct 2026) and are not in this public repository; `artwork/panel_a.png` is the panel-a band cropped from the v30 composite at 300 dpi. Without the controlled directory the v30 composite checks print SKIPPED.
- The fold count (ten) is a design parameter of the upstream cross-fit (capscan_crossfit_10fold); no code here uses it.
- The colour bar of panel i reuses panel-h layout geometry (computed without any data value) and two page constants measured by the fig_v27 composer (`I_DXI`, `PAGE_H_V27`).
- The Fig. 4h feature selection (19 of 20 features with any positive value) comes from an earlier exploratory screen whose selection code was not located (OPEN_ITEMS).
- The v28f and v30 composites are not rebuilt by the notebook; ASSEMBLY.md documents each step and `expected/Fig4_v28f.*` and `expected/Fig4_v30.*` hold the results. v30 differs from v28f only in the panel-h track label (`checks/v30_vs_v28f.json`).
- v30 sync (2026-10-06): `NUMBERS.tsv` carries the four v30 columns (`v30_para`, `v30_quote`, `v30_printed_value`, `v30_status`; paragraph index of the v30 visible text). The regenerated b,c and h are compared with the v30 composite (`check_against_v30`); both differ from it only by rendering, not content (OPEN_ITEMS PIXEL).
