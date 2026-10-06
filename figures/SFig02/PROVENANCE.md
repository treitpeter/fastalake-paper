# Provenance: Supplementary Figure 2

`REV` below stands for
`LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24`
and `CR3` for
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03`.

## Panels

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a-d (with title) | `REV/supplements/revision_2026-09-24/build_revision.py` (sha256 `e636bd85164896f53ec1d13ef671754ff800597d2d955a3eb01f5151dced1bbe`), selection block lines 61-87 | `REV/supplements/revision_2026-09-24/delivery/figures/SELECTION_SPECIFICITY.pdf` | `f3b4d4a1c620eebbbc979ce2ec46fce2ff2728aa6db79accaa0a45c06496c7bc` | not recorded (run 24 Sep 2026, 14:32) | Matplotlib 3.10.8 (PDF metadata); interpreter not recorded; DejaVu Sans |
| a-d (manuscript image) | `REV/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/prepare_supplement_figures.py` (sha256 `5c52b4b89037916a4503e3b77e8454522186a83c0fdff48cced5f1776141ce4d`), with `SUPPLEMENT_FIGURE_SOURCES.json` key "1" (sha256 `69943ed802e30b3182777be1f2a8347afeb7b306fdc82ffff18320a5431c5f15`) | `REV/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/release_candidate/supplementary_figures/Supplementary_Figure_01.png` (and `.pdf`) | png `9dacf0331290f5719b51bc50aff870730fa8318c797d236040d3bb1195be1154`; pdf `b4148cb2513c2972c650bf2eb3eea59948144aa261727b4107920536212a3d6e` | not recorded (run 25 Sep 2026, 15:51) | MuPDF 1.27.1 (PDF metadata) |
| e,f (vector) | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/bags_prominence_review_2026-09-30/integrate_v18.py` (sha256 `610b761e1a7be82e6257ddd53d1beb8dfba61111dfd95e2932190540aa4944d8`) executing lines 46-51 of `audit_and_plot.py` (sha256 `a2de46ea6bb2940443ce2a4555286275553c87a8b65d6d4dfa044691dcd86849`) with letters e,f | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/bags_prominence_review_2026-09-30/Supplementary_Figure2ef_bags.pdf` | `1f84736efe2a2efa85d270246dd72763ba8a33a5cb8a760e011b2cc5d4a88a00` | not recorded (run 30 Sep 2026, 12:40) | Matplotlib 3.10.8 (PDF metadata); DejaVu Sans |
| e,f (manuscript image) | `CR3/iteration3/source_data/figures/reexport.py` (sha256 `cab43851449afc661d5846044573f1a9dec78850be646e1daa43e40714d6f0d1`) | `CR3/iteration3/source_data/figures/out/SupplementaryFigure2_continued_v21_reexport.png` | `f4231f5f090a014bf19e764d156067577c597df499d6efdbb92f31b3705ce806` | 6048949 (control-only run 6048946) | CPython 3.12.4, PyMuPDF 1.27.1, Pillow 12.0.0 (miniconda interpreter named in `reexport.sbatch`) |

Regeneration in this repository (CPython 3.12.4, matplotlib 3.10.7, PyMuPDF 1.27.1, pandas
2.3.3, numpy 2.2.6, Slurm, 4 cores, 0.44 GB peak RSS, 6 s): all four pixel checks `IDENTICAL`
(`checks/*.json`), and both manuscript images byte-identical to the originals.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| a-d | Supplementary Figure 2, `word/media/image8.png` | BYTE_MATCH to `Supplementary_Figure_01.png`; producer_status COPY_ONLY; cited script `v15_MM_coauthor_2026-09-25/build_figures_v15.py:36` only copies the image; "original plotting/rendering producer unresolved" | written by `prepare_supplement_figures.py` (line 15) from the PDF written by `build_revision.py` (line 51, selection block); regenerated here pixel-identical and byte-identical | partly | same image and same origin file; this unit resolves the producer the census left open |
| e,f | Supplementary Figure 2, `word/media/image9.png` | BYTE_MATCH; WRITE_SITE_HASH_RECEIPT, `reexport.py:126`, receipt `REEXPORT_VALIDATION.json`; "producer was not rerun" | same write site; vector input written by `integrate_v18.py` (exec of `audit_and_plot.py` plotting block); both steps regenerated here, byte-identical | yes | this unit adds the rerun the census did not do |

`REBUILT_CANDIDATES.tsv` lists no rebuilt candidate for Supplementary Figure 2.

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| `selection_length_lookup.tsv` | `REV/supplements/revision_2026-09-24/source_data/original_S01/63594d955d35_MM_review_analysis__A02_4_2_lookup_specificity__results__length_by_arm.tsv` (= `LAB_ROOT/MM_review_analysis/A02_4_2_lookup_specificity/results/length_by_arm.tsv`, job 5920360) | `63594d955d35bd564a32fe68317274537519376858d3bc0daee1b09cc3b1cc29` | byte copy. Package `Figure1/length_by_arm.tsv` is byte-identical; package `Supplementary_Figure_02/selection_length_lookup.tsv` (rewritten by `build_revision.py`) has the same values, differing by at most 1.1e-16 in five float columns |
| `lookup_query_outcomes.tsv` | `LAB_ROOT/MM_review_analysis/A02_4_2_lookup_specificity/results/query_source_data.tsv` (job 5920360) | `64182cf62ebb7d9992a0ffa2541d51590ead3a0f780dca2717f1e157e9100eba` | byte copy; not in the Source Data package. Recounted to every row of `selection_length_lookup.tsv` by `panels.py:stat_ab_lookup` |
| `selection_grid_acquisitions.tsv` | sheet `Per_acquisition` of `REV/supplements/revision_2026-09-24/source_data/original_S23/68712926f1b0_SupplementaryTable_selection_sensitivity.xlsx` | `68712926f1b020c19d313a35767031b589f07a0a19dca5a9104dd68764800008` (workbook) | sheet written as TSV (pandas); reads back exactly. Byte-identical to package `Supplementary_Figure_02/selection_grid_acquisitions.tsv`. All 1,102 rows equal the corresponding rows of the upstream sweep output `LAB_ROOT/canonical/publication/analyses/length_score_sweep_hstool_2026-09-06/tables/per_sample_by_arm.tsv` (sha256 `ba1aebd162a2d0c3c122e3482af53ea3c73a2eafbcdfa4fc8fd42f0970ca55e5`, written by `sweep_tables.py`) in peptides, selected proteins, protein groups, seeds and decoy fraction; that table has two further k-mer arms not used here |
| `selection_grid_summary.tsv` | sheet `Arm_summary` of the same workbook | as above | as above; byte-identical to package `Supplementary_Figure_02/selection_grid_summary.tsv`; medians and quartiles recomputed from the 1,102 rows by `panels.py:stat_cd_grid` (largest relative difference 8e-17) |
| `bags_bounded48_per_acquisition.tsv` | `rows` of the 48 `LAB_ROOT/MM_review_analysis/A04_4_2_expansion/alphanovo_bags_2026-09-23/bounded_sweep_2026-09-25/outputs/*/*/SEARCHES_COMPLETE.json` (48 Slurm jobs, ids 6009200 to 6009754) | sha256 of the concatenated per-file sha256 values in sorted path order: `cc766aff3411a4c55a0407b76ca54e3e11b7af72c85d57618d7a32e444755217` | rows concatenated as `audit_and_plot.py` line 16-21 did; byte-identical to package `Bag_aware_lookup/bounded48_per_acquisition.tsv` |
| `bags_bounded48_orderings.tsv` | `bag_orderings_hist` and `lookup_eligible_predictions` of the 48 `.../bounded_sweep_2026-09-25/outputs/*/*/BOUNDS.json` (written by `run_bounded.py`, sha256 `1173ca7e46eda476b01921cf8aa435484febb364fb6887cff02288506c751029`) | concatenated per-file sha256, sorted: `8ad04dfe86553969a8249ecd86b9e5be08042921c9f466f908e39eed26e0c7eb` | histograms written one row per acquisition and ordering count; not in the Source Data package (it was the untraced source of the 908,037 / 976,861 figure on v28d line 561) |
| `bags_fragment_audit_summary.tsv` | `LAB_ROOT/MM_review_analysis/A04_4_2_expansion/alphanovo_bags_2026-09-23/gain_loss_2026-09-24/fragment_validation/audit/LOCAL_DIAGNOSTIC_SUMMARY.tsv` | `d9ed02c344d5b22427667505e80a420ddc59f7aa3e1d21e5d1f36391efea3596` | the one PREDICT acquisition file name replaced by `PREDICT_acq_01`; otherwise unchanged. The six-row subset used by panel f equals package `Bag_aware_lookup/fragment_order_subset.tsv` after the same replacement |

The extraction ran on Slurm with a script kept outside the repository (work directory,
`figures/SFig02/extract_source_data.py`, receipt `EXTRACT_RECEIPT.json`). Package files not used:
`Supplementary_Figure_02/selection_length_matched_pool58.tsv` (a row subset of the summary that
`panel_abcd` selects itself), `Supplementary_Figure_02/reference58/` (roster audit of the 58-pool
rerun, not plotted), `Bag_aware_lookup/bounded48_summary.tsv` (recomputed by `bags_summary`, equal),
`Bag_aware_lookup/experiment_status.tsv` and `Bag_aware_lookup_source_data.xlsx` (execution
snapshot and workbook copy of the same rows).

## v30 update (2026-10-06)

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a-f v30 | `FIN/v30/figures_style/SFig02/build_sfig02_style.py` (sha256 a8c181e1794b80dc386b028ea09d6909eea3f8f8e4520cfdc1432ef47bc1166b) with `FIN/v30/figures_style/_shared_sf01_03/restyle_common.py` (0c7ced18eb696d77db2639e5121fb66399debe60684097d15cd59cff816addb1); it imports this unit's `panels.py` | `FIN/v30/figures_style/SFig02/SupplementaryFigure2{,_continued}.pdf/.png`; receipts `values_check.json`, `pdf_audit.json` | see ASSEMBLY.md, v30 | not recorded | paper environment, matplotlib 3.10.7, PyMuPDF 1.27.1 |

`FIN` = `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`.
The redraw and helper code were moved into `panels.py` unchanged apart from repository-relative
paths and `fl_style` from `style/`. Notebook run on Slurm (job 6075247, 25 s): both v30 PNGs and
all four v28f checks IDENTICAL. Source data and numbers are unchanged; the v30 legend text
(v30f) is unchanged from v29 (LEGEND_CHANGES.md: none needed).
