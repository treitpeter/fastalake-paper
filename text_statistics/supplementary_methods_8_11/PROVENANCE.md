# Provenance: Supplementary Methods 8-11

`LAB_ROOT` is the lab project directory. `C2` abbreviates
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15`;
`PKG` abbreviates `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/figure_sd/source_data/package/FastaLake_Source_Data`.

## Panels

This unit has no figure panels. The statistics come from these producers.

| statistic | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| CapScan regional tests, own-stool test (347) | `C2/scripts/capscan_fig5bcd.py` | `C2/fig5/fig5bcd_347/four_region_tests.tsv`, `similarity_summary.tsv`, `own_stool_participants.tsv` | 885df8ef9cef7037...9736, ac61627e0eb03c78...9a39, 67f22f573d418bde...ff0e | 6063598 | Python 3.12.4, numpy 2.2.6, scipy 1.16.3, pandas 2.3.3 |
| CapScan cross-omics (347) | `C2/fig5g/crossomics_capscan.py` | `C2/fig5g/crossomics_347/{functional,native}_correlations.tsv` | 4c3468fc11199d32...e3f8, 0004f35d52242258...5305 | 6063849 | same |
| CapScan controlled distance (347) | `.../final/reruns_347/A4_controlled_distance/code/distance_capscan347.py` | `.../A4_controlled_distance/controlled_176_comparisons_347.tsv` | e0757c104a4e1452 (prefix; equals source_data/SFig11 copy to 1.1e-16) | 6064864, 6064865 | same |
| PREDICT cross-omics, individual functions, EDS screen, held-out EDS prediction | `.../iteration1/fix_eggnog_native/downstream/verified50/B/{crossomics,individual_functions,eds_biochemistry,nested_eds}.py` | `PKG/Biology_refit/{cross_refit,individual_refit,eds_refit,nested_refit}/` (identical to verified50 outputs) | nested_refit/model_metrics.tsv df5b688227c6083d...e7f | 6048399 (nested refit) | Python 3.12.4, numpy 2.2.6, scipy 1.16.3, statsmodels 0.14.5, scikit-learn 1.7.2 |
| Antibiotic endpoints | `.../iteration1/fix_rifaximin/rifaximin_sensitivity.py` | `source_data/STab3/endpoints_by_definition.tsv` (STab3 unit) | d4273f1aefad9ca3...f3f | see tables/STab3/PROVENANCE.md | see tables/STab3 |
| Histidine and terminal families (347) | `LAB_ROOT/MM_review_analysis/A24_4_7_terminal_inference/frozen/search.py`, `A28_4_7_histidine_inference/frozen/collect_revision4.py` | `C2/fig6/A24_fixed347/results/CAPSCAN/complete_statistical_families.tsv`, `C2/fig6/A28_fixed347/results_revision4/full_CAPSCAN/complete_statistical_families.tsv` | see figures/Fig6/PROVENANCE.md | 6058149, 6058108 | Python 3.10.9 (Anaconda 2023.03), numpy 1.26.4, scipy 1.12.0 |
| Non-specific searches | `LAB_ROOT/MM_review_analysis/A05_4_7_unspecific`, `C2/fig6/A05`, `LAB_ROOT/MM_review_analysis/A13_4_7_predict_unspecific` | per-acquisition `inputs/*__unspecific.json` and `results/*/COMPLETE.json` | not hashed individually (653 files); summary 96f06a076849ce5d...0bf | 6058023 (participant 15); production see A05/A13 status | Sage 0.14.6 |
| Healthy-stool grouping counts | grouping receipt `.../revision_2026-09-24/delivery/source_data/compute_records/HSTOOL_grouping_COMPLETE.json`; KO refit receipt `.../jess_integration_2026-09-29/source_data/current/HSTOOL_KO_refit/upstream_receipts/COMPLETE.json` | the receipts themselves | cd1f1bbb0756190220d0151a756edb8df073d28f929056a1177d51be598a1fb8; 1c09c40e6dd48cb87e6a4229287aa2b517c410456567827a8408f6283ba3d253 | not recorded | not recorded |

This unit's own computation ran on Slurm: `prepare_source_data.py` (srun, 4 CPUs, 1.4 GB peak, 9 s),
`write_tables.py` and the notebook (srun, 4 CPUs, 0.28 GB peak, 10 s), with the repository environment
(Python 3.12.4, numpy 2.2.6, scipy 1.16.3, pandas 2.3.3, statsmodels 0.14.5, scikit-learn 1.7.2).

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| none | not applicable | not applicable | not applicable | not applicable | text unit; the census covers embedded figure images only |

Independent Codex statistics reimplementation (`CODEX_TASKS_2026-10-05/out_stats/REPORT.md`), compared,
not reused:

| statistic | Codex value | this unit | agree? |
|---|---|---|---|
| CapScan regional mean differences, 15 participants | 12.10865315, 11.12783135, 8.392889373, 4.669726646 | identical to the producer table (max difference 0) | yes |
| CapScan regional exact BH q (32,768 patterns) | 0.002766927083 (three), 0.03796386719 | identical | yes |
| Own-stool advantage, P | 0.06073191429, 15 of 15, P = 1e-05 | 0.06073191428677263, 15 of 15, P = 1e-05 | yes |
| Antibiotic q (human share, protein groups) | 0.00061, 0.00124 (own random stream); exact enumeration 0.000655, 0.00117 | producer 0.00068, 0.00126 (family structure recomputed here; tests in tables/STab3) | same conclusion; differences are Monte Carlo stream differences |

Codex statistics register (`out_statsreport/STATISTICS_REGISTER.tsv`, `MISSING.md`) rows S055-S106 and
S110-S125 list missing seeds, implementations, statistics and family sizes for these sections; STATS.tsv
fills them from the code where the code answers, and OPEN_ITEMS.md lists each as TEXT-GAP.

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| own_stool_distance_matrices.tsv | `C2/fig5/mc/profiles/CAPSCAN_noclust/OG_other_unresolved.npz` + `_axes.json`, `C2/roster/design_347.tsv` | e277cb600c7ee2b6f6e94b0bbd0f86a32a90878cd1b082045c77ed55b7b946a1; 3c3544a93f32ad1eacd55e3f1a9758265cf724fa3004aed0f13c3999e81af955; 30f70233143a370755d1bc4d1dd8f236d94a2fe4b9b27cc10c64e94602461000 | `prepare_source_data.py:own_stool` (producer code; diagonal and row means asserted equal to `own_stool_participants.tsv`). Package sheet Biology_refit/capscan_refit holds the superseded 14-participant values (53,670 / 44,129 OGs); not used |
| capscan_dna_runs.tsv | `LAB_ROOT/analyses/capscan_metag_analysis_2026-09-10/manifests/runs.tsv`, `aggregate/metaG_{OG,KO,Pathway}.tsv`, `.../cohort_reconciliation/crossomics/native/current331_plus16_all/COMPLETE.json` | 0edba605...6306; 82bc95cf...279c; 17b81e78...eab; 8d13d06f...763e; f51e05e9...8726 | `prepare_source_data.py:dna_runs` (columns selected; flags computed). Not in the package |
| predict_native_feature_summary.tsv | `.../jess_integration_2026-09-29/source_data/current/Crossomics/full/native_SOURCE_DATA_feature_summary.tsv` = package sheet `current/Crossomics/full/native_SOURCE_DATA_feature_summary.tsv` | 9b3d206a69e2e5c4153bee1316510190f02558af768214988f84307c4357233f | `prepare_source_data.py:native_features` (PREDICT DNA and RNA rows; values agree with the package) |
| eds_prediction_model_metrics.tsv, _permutations, _increment, _selected_alphas | package `Biology_refit/nested_refit/{model_metrics,permutation_results,increment_over_baseline,selected_alphas}.tsv` | df5b6882...e7f; 16e32d66...0a0; c5ba12b2...b75; 127cbda3...8bc | byte-identical copies |
| eds_prediction_heldout.tsv | package `Biology_refit/nested_refit/heldout_scores.tsv` | 50806b84ddd080042bad284ce70c56b20572bec0320c91ab946598af60601060 | `prepare_source_data.py:eds_prediction` (specimen and enterotype dropped, row position added; values agree) |
| eds_prediction_outer_folds.tsv, _fit_audit, _distance_reference | package `Biology_refit/nested_refit/{outer_folds,fit_boundary_audit,distance_reference_folds}.tsv` | f8047ef3...e36; cf2115b8...752; f3127cfe...eb | patient and specimen ids dropped; row lists replaced by counts |
| nonspecific_search_summary.tsv | per-acquisition search configurations and completion receipts (A05, C2 A05, A13) | not hashed individually | `prepare_source_data.py:nonspecific` (counts; settings asserted identical). Not in the package |
| hstool_grouping_receipt.tsv | the two healthy-stool receipts above | cd1f1bbb...fb8; 1c09c40e...253 | `prepare_source_data.py:hstool_receipt` (values read). Not in the package |

`SOURCE_DATA_MAP_reanchored.tsv` rows for these sections: SM-F5c0/1/2 still point to the 14-participant
`Biology_refit/capscan_refit/similarity_summary.tsv` (72,761 / 53,670 / 44,129 OGs); the v28d text prints
the 347-roster values (73,316 / 53,023 / 44,330), which the Fig5 unit recomputes. M-XOMICS (120 functional,
32 native results) agrees; here the 120 and 32 are counted from the PREDICT rows (SFig11) plus the 347
CapScan rows (Fig5). ST5-n (30,581) and M-EDS rows agree with the STab5 and SFig10 units.
