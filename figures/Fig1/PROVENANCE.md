# Figure 1: provenance

`FINAL` = `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`.
`FAR` = `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/four_author_integration_2026-10-02/actual_FAR_261002`.
`V10` = `LAB_ROOT/MM_review_analysis/figure_restoration_2026-09-16/Figure1`.

## Panels

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a (v30) | `FINAL/v29/figures/Fig1_v2/build_fig1a_funnel.py` (sha256 `de4b12164674b5eddb6e28f4f55ed4a7373b04d086fefa7d13b754c36e8abee0`), ported as `artwork/build_fig1a_funnel.py` (sha256 `60d054a8ef50554a7f98bcc506ff69188f1654cf66b1c25766f30c0a66cfd20f`, paths only); inputs `FINAL/fig1a_perfecter/Figure1.pdf` (8bc5a81b...), `FINAL/fig_v27/Figure1.pdf` (432699de...), `FINAL/v29/figures/Fig1/Figure1.pdf` (cb28f0d2...) | `FINAL/v29/figures/Fig1_v2/Figure1.pdf` and `Figure1.png` (receipt `_build_receipt.json`, sha256 `5eaf0fe166cb117fc9653e8c155f717eb7a9a10e1b4d9c191c67839e8f99b3f6`) | PDF `0f686c53a29d7c360ace01aceae7c2918d8aed450f40311c4016c35e4e6cc23f`; PNG `b1481ff689fb4439824ec5ec4f567b03e2d14860fa2445e6dfbabb94cd4fcb95` | producer: not recorded (CHANGES.md: run on Slurm); repository rebuild: 6075277 (notebook) | producer: miniconda python, PyMuPDF; rebuild: CPython 3.12.4, PyMuPDF 1.27.1 |
| a (v28f base) | chain of edit scripts ending in `FINAL/fig1a_perfecter/build_fig1a.py` (sha256 `26714201117a52ccf6ddc0c0aeec77a45b7df110dbaab61a02c18f305a13ccfa`); see ASSEMBLY.md | `FINAL/fig1a_perfecter/Figure1.pdf` (panel a = CropBox 0, 0, 518.74, 372 pt) | `8bc5a81b489f1ab5fb7ec23f695c0cc60e92d26346faff4255973acb952c5128` | not recorded (CHANGES.md: run with miniconda python) | PyMuPDF 1.27 |
| b | `FINAL/fig1bd_n15/code/build_panels_bd_n15.py` mode `n15` (sha256 `45c06b21f0252e53feb057d74b3d2d92b5665628b97464c1e7ddb4ef108b73c5`), style `LAB_ROOT/figure_style/fl_style.py` (sha256 `58b64c343f095cbace09c3e9766e3e78060c1c8cf522b075c4047300322164c3`) | `FINAL/fig1bd_n15/panel_b.pdf` | `1fcc8730a6f64fb6320ad1027f10e08256bfdd6de11a2f22ea0d2e302a3e916c` | 6064327 | anaconda 2023.03, CPython 3.10.9, matplotlib 3.10.8, numpy 1.26.4, pandas 2.3.3 |
| c | `FAR/build_FAR12.py` (sha256 `757b492fe3224100eb4319dea63e6032384ea767079e52834cfe66513dd1e8af`), same style module | `FAR/FAR12_figure1c_applied/Figure1c.pdf` | `ecff68781c968ced8e2ef9439f2226383fb7d114a44e774bfc900a1349289671` | not recorded | not recorded |
| d | as b | `FINAL/fig1bd_n15/panel_d.pdf` | `88cfb3382ef98b9d8fcee01d7617c7b143e0a8fc7109819d3591a307753562fe` | 6064327 | as b |
| e | `V10/panels_bcd_v10/build_panels_bcd_v10.py` (sha256 `770b5add04801b4773a80431db3b9caf4ef976ca31b91aee10476b3f16d351c3` per its receipt); the panel-e code is unchanged in `build_panels_bd_n15.py` | `V10/panels_bcd_v10/Figure1_panels_bcd_v10.pdf` (four-axis row; panel e cropped) | `1dc4f23dde58a06724099b56609909b77e0753daa498763853119d322636c064` | 5966808 | anaconda 2023.03, CPython 3.10.9 with user-site matplotlib (log) |

Regeneration here: `panels.py` with CPython 3.12.4, matplotlib 3.10.7, numpy 2.2.6, pandas
2.3.3, PyMuPDF 1.27.1 and the Söhne fonts. Result at 300 dpi (`checks/*.json`):

| check | new | expected | verdict |
|---|---|---|---|
| panel_b | build/Fig1/panel_b.pdf | expected/panel_b.pdf | IDENTICAL, 0 differing pixels |
| panel_c | build/Fig1/panel_c.pdf | expected/panel_c.pdf | IDENTICAL, 0 differing pixels |
| panel_d | build/Fig1/panel_d.pdf | expected/panel_d.pdf | IDENTICAL, 0 differing pixels |
| panel_e | build/Fig1/panel_e.png (crop of the regenerated row) | expected/panel_e.png (same crop of the v10 row) | IDENTICAL, 0 differing pixels |
| panel_e_expected_derivation | crop re-derived from expected/panels_bcde_row_v10.pdf | expected/panel_e.png | IDENTICAL |
| Fig1_v28f | artwork/Figure1_v28f.pdf at 300 dpi | expected/Fig1_v28f.png | IDENTICAL |
| Fig1_v30 | artwork/Figure1_v30.pdf at 300 dpi | expected/Fig1_v30.png | IDENTICAL, 0 differing pixels |
| Fig1_v30_rebuild | build/Fig1/funnel/Figure1.pdf (artwork/build_fig1a_funnel.py rerun, Slurm job 6075277) | expected/Fig1_v30.png | IDENTICAL, 0 differing pixels; the rebuilt 300 dpi PNG is also byte-identical to the producer PNG (b1481ff6...) |

The panel checks are pixel-exact only with the producers' style values, which `panels.py`
applies (see its module comment and `OPEN_ITEMS.md`). With `style/fl_style.py` alone,
panels b, c and d differ in about 2,600 to 3,000 pixels each, all on axis spines, ticks and
the panel letter (1 pt against 0.6 pt lines; bold against regular letter).

Independent check of the 4,258 accepted peptides absent from the predictions (Results line
63): the predictions CSV (sha256 `8826920da56b6ed2fa443f519270c830f02294e40f36a7c497741f64b13cb15d`)
and the Sage result `results.sage.tsv` (sha256
`92339fc9c159bd00ba24c0edd07e3955120962a98b1d0951c3a67dc012a391d4`) of the worked
acquisition were read directly with the selection rules of A09 `analysis.py` (target,
rank 1, peptide q <= 0.01, modifications stripped, I to L) on Slurm job 6068960
(2026-10-06): 189,806 prediction rows, 141,766 distinct predictions, 29,663 accepted
peptides, 4,258 not predicted. Same as the A09 record.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| Figure 1 as embedded in v28d | Figure 1, `word/media/image1.png`, sha256 26ec1aaf... | BYTE_MATCH to `FINAL/fig_v27/Figure1.png`; write site `FINAL/fig_v27/build_fig1_3_v27.py:72`; producer not rerun | same file and script (ASSEMBLY step 8) | yes | this is the superseded image: old panel a, CapScan n = 331 in b and d |
| Figure 1 v28f candidate PDF | REBUILT_CANDIDATES.tsv, `FINAL/fig1a_perfecter/Figure1.pdf`, sha256 8bc5a81b... | write route `FINAL/fig1a_perfecter/build_fig1a.py:523`, inspected, not rerun | same file (ASSEMBLY step 9); shipped as `artwork/Figure1_v28f.pdf` | yes | Codex did not rerun it; neither did this unit. Panels b to e inside it are reproduced pixel-exactly from source_data |
| Figure 1 v28f candidate PNG | REBUILT_CANDIDATES.tsv, `FINAL/fig1a_perfecter/Figure1.png`, sha256 6e92d894... | same write route | same; shipped as `expected/Fig1_v28f.png`; equals the 300 dpi raster of the PDF | yes | none |

The Codex statistics register rows S001, S002 (entrapment) and S131 to S133 (scrambled
fractions) agree with `STATS.tsv`. The gaps it lists (entrapment seed values and summary
implementation; software for the scrambled fractions) are filled in `STATS.tsv` from the
code. The Codex number recheck (`out_recheck/REPORT.md`) has no row for a Figure 1 panel
or for Results lines 51 to 65; nothing to reconcile there.

## v30 artwork files

| file | origin (LAB_ROOT/...) | sha256 | note |
|---|---|---|---|
| artwork/Figure1_v30.pdf | `FINAL/v29/figures/Fig1_v2/Figure1.pdf` | `0f686c53a29d7c360ace01aceae7c2918d8aed450f40311c4016c35e4e6cc23f` | byte-identical copy; the v30 manuscript figure |
| expected/Fig1_v30.png | `FINAL/v29/figures/Fig1_v2/Figure1.png` | `b1481ff689fb4439824ec5ec4f567b03e2d14860fa2445e6dfbabb94cd4fcb95` | byte-identical copy |
| artwork/Figure1_v27_art_source.pdf | `FINAL/fig_v27/Figure1.pdf` | `432699dea4b0f12142de3ef76c96e9adffa8bb26b3ac7b7debaf336a66e6140d` | build input (clipped icon and peptide bars) |
| artwork/Figure1_v29_intermediate.pdf | `FINAL/v29/figures/Fig1/Figure1.pdf` | `cb28f0d27cfcae47f432e7357aa0ad5f9b5b2fee6daed1774f60efff8d57a968` | build input (outside-column assertion only) |
| artwork/panel_a.pdf | crop of artwork/Figure1_v30.pdf, CropBox (0, 0, 518.74, 372) pt, PyMuPDF `set_cropbox` | `26dd41a212b8e277cdf7666ab770afa0c27f107a22b89c742522bfc61ef6ae7c` | replaces the v28f crop (92a13d81...) |

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| acquisition_sizes_726.tsv | `FINAL/fig1bd_n15/acquisition_sizes_726.tsv`; replaces package sheet `Figure1/acquisition_sizes.tsv` (710 rows, CapScan 331). The 710 shared rows are byte-identical in both. | `fd9ade3bb4ebd2ee6c129f69819b9697cf20abd7c10c77a12644bef25396e386` | copy (`figures/Fig1/prepare_source_data.py`) |
| acquisition_support_summary_726.tsv | `FINAL/fig1bd_n15/acquisition_support_summary_726.tsv`; replaces package sheet `Figure1/acquisition_support_summary.tsv` (2,130 rows, CapScan 331) | `6f81da0a03893bcdf2ed5406ddb0e52eb16791cfd46f55a3253b31c65ac7dec2` | copy |
| length_by_arm.tsv | `LAB_ROOT/MM_review_analysis/A02_4_2_lookup_specificity/results/length_by_arm.tsv`; byte-identical to package sheet `Figure1/length_by_arm.tsv` | `63594d955d35bd564a32fe68317274537519376858d3bc0daee1b09cc3b1cc29` | copy |
| five_seed_summary.tsv | `LAB_ROOT/MM_review_analysis/A02_4_2_lookup_specificity/results/five_seed_summary.tsv`; byte-identical to package sheet `Figure1/five_seed_summary.tsv` | `ebb3de14cd36ab997a1e0e2245e882de30f822b3b06d78cbf42e8af79bf73cd6` | copy |
| donor_seed_means.tsv | `LAB_ROOT/MM_review_analysis/A11_4_2_control_reconciliation/results/donor_seed_means.tsv`; byte-identical to package sheet `Figure1/donor_seed_means.tsv` | `5c6b0af7c7224b7f1968b82a8a2c22ecae4d25679f691ef04c0957644c27ea66` | copy |
| donor_recount.tsv | `LAB_ROOT/MM_review_analysis/A11_4_2_control_reconciliation/results/donor_recount.tsv` (per-seed values; not in the package) | `4e51dd75c37deed3f13c84ad2407032b01ae16d877efac9d1b60a174aafb5ea8` | copy |
| scramble_end_to_end.tsv | package sheet `Review_v16/Terminal_scramble/scramble_end_to_end.tsv` (= `LAB_ROOT/review/v15_integration/evidence/scramble_end_to_end.tsv`) | `844f92521718387c35fea524b5cc8502cac8847e0064dfea04c27637293ce865` | copy |
| selection_grid_acquisitions.tsv | package sheet `Supplementary_Figure_02/selection_grid_acquisitions.tsv` | `8b472997023e29d4e24b5a9f7a06d9aceda6f45fc4c45391825c075d85aa3425` | copy |
| stage1_stats.json | package sheet `Supplementary_Figure_02/reference58/stage1_stats.json` | `2be7eafbd02937b6df5ac05e1845c1a4d1d64203bfc62bd77ebb23739a21fcd8` | copy |
| human_reference.json | `LAB_ROOT/MM_review_analysis/A01_4_1_database_sizes/results/human_reference.json` (not in the package) | `b705dfc4aba8718c83939a812e3de6614ce0a72eb9260c3e0c7b6ff632cf1079` | absolute path reduced to file name, `mtime_ns` dropped (`prepare_source_data.py:reduce_record`) |
| worked_acquisition_prediction_support.json | `LAB_ROOT/MM_review_analysis/A09_4_1_prediction_support/results/HSTOOL/20181231_QX7_PhGe_LC95_SA_WorkflowReplicates_2/COMPLETE.json` (package map row RES-17 says "NOT IN THIS ARCHIVE") | `f4a34959914caa7dacbb3a3a9f85445d214cfbf6edf10f98134ba92e9c0fb066` | as above |

Agreement with `SOURCE_DATA_MAP_reanchored.tsv`: rows RES-01 to RES-16, RES-18, RES-19,
FIG1-c, FIG1-d and FIG1-e give the same values as this unit. FIG1-b points to the 710-row
sheet, which this unit replaces with the 726-row table. RES-17 has no source in the
package; this unit adds it. RES-10k (ratio 11,428.6) agrees with
`stat_a_reservoir` (11,428.598...). The 36,588 and 5,841 typed in panel a have no map
row; both are traced in `ARTWORK_NUMBERS.md`.
