# Provenance: SFig07

`J` below is `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/joint_search`.
`A21` is `LAB_ROOT/MM_review_analysis/A21_7_hstool_quant_comparison`.
`CQR` is `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/clinical_quant_review`.
`PKG` is the Source Data package folder `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Figure_07`.

`FIN` is `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`.

## Panels

v30 figure (current):

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a-d (one canvas), v30 house style | `FIN/v30/figures_style/SFig07/build_sfig07_style.py` (sha256 648e38efeb5cb4d6093e73d7352ce15da5e7e03867216fcddeed6c0acafd2a5c) | `FIN/v30/figures_style/SFig07/SupplementaryFigure7.png`; `.pdf` | ab999c1d525baddfd84e2357c5271faa12180e9fea793c00f4994ca7739c19e3; c364e69ef304c0695f21faa2b11bd414d76fe3c4e997830701d2c5eb1494e4f4 | srun via `FIN/v30/figures_style/_work_sf07_09/run.sh` (job id not recorded) | `fastalake-paper_v28_work/env` python; fl_style with Söhne from the lab styleguide directory |
| a-d, 9-50 producer drawing (v29) | `FIN/len950/SFig07/panels_len950.py` (37898232f850eae47d91734eb1ef5c93bdf638ef4cfc22533393be3bb8891b0a), driven by `FIN/len950/SFig07/build_sfig07_len950.py` (da9c0a216ca479325e6f91f9268b9aa5efae3a4915103e88bd5c1be24c595e22) | `FIN/len950/SFig07/Supplementary_Figure7_len950.png` | 7bb2e76578f216a49c49cb2d44bd0430952a65e3175dd9701371c57a9e7af8c2 | inputs from 6072769 (`hstool_recount.py`, 2a497356...) and array 6072755 (`clin_recount.py`, 7a9c331b...) | as above |

The repository render of the v30 figure is IDENTICAL to the v30 PNG and to the v30 PDF at
300 dpi. The rows below describe the v28d producer, whose drawing code is kept as
`panels.py:render` and run unchanged inside `render_v30` for the value comparison.

v28d producer (superseded):


| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a, b, c, d (one canvas) | `J/code/build_supp_fig7.py` (sha256 56d24cf74641dc8ad4c42c26e49982cb5bfa61457d392e86060c51dfa1ed0de0) | `J/figure/Supplementary_Figure7_revised.png` | b88113c0c9aa418ef52e86ba1309c30d53e6e447fa14002087f34950fc98c47b | none recorded for the render (run interactively 2026-10-03); its inputs from jobs 6048637 (array 0-2, extraction) and 6048638 (comparison) | matplotlib 3.10.7 (PNG and PDF metadata); font `soehne-buch.ttf` (FIGURE_VALIDATION.json) |
| a, b, c, d (PDF) | same | `J/figure/Supplementary_Figure7_revised.pdf` | a01bc44b02c623f87aa1c90db897a65f21b035a381b58bb7f192065b998f20c9 | as above | as above |
| a, b, c, d (SVG, not used here) | same | `J/figure/Supplementary_Figure7_revised.svg` | 9ac62c16ca38e0f7f8c71629e0fd83b63a42e14b43f413a4d1a6b2b9d2bdfae8 | as above | as above |
| a (study-level columns) | `A21/frozen/compare_revision2.py` (sha256 4e959c46a78392ae086e9d1669f27efa1c57ec0d014da89123cd98d83ad7dd47) | `A21/results/identity_comparison.tsv` | c5a863acb448d698d0a89b440e5c0a42064e55a5475b0de3d793f6f836287090 | 5928104 (compare_revision2) after 5928103 (extract_revision2) | CPython 3.10.9, numpy 1.26.4, pandas 2.3.3 (anaconda 2023.03) |
| a (acquisition-q column), c, d | `J/code/extract_fileq.py` (b816c095ef001c1d7ae9b161a7403738861bab46e131b08eb9511e54d4249990), `J/code/compare_fileq.py` (fbbaaeb90fbdc5fb221efc6dd3efee71006d946ba7f822b9307929817088a4d6) | `J/out/compare/identity_per_acquisition.tsv`, `J/out/compare/microbial_group_pairs.tsv` | 56c7271c5092a44bf665f2b6ba441da79e463694149a51d1bfb877ec6e28c23a, e8d58907fb5564149344d76b09e6dbe8b3dd3062e83084ed14da083143471329 | 6048637 (array 0-2; peak RSS 1.9-2.7 GB, 3-3.5 min each), 6048638 (peak RSS 1.5 GB, 1 min 18 s) | CPython 3.12.4, numpy 2.2.6, pandas 2.3.3 (miniconda3) |
| b | `CQR/code/summarize_v2.py` (ca0b5eb5eedcf1a4cf3f6a9716f1cc3b67b86ee29748854f36c13fb6b01b3537) | `CQR/PREDICT/per_participant.tsv`, `CQR/CAPSCAN/per_participant.tsv` | b392e6c677b8f4157b3703287bb0535e121979a47f248475efc8c45263aaa4e2, 68e2b5b21d71c7c43c317c727639af18104ee0fd74781cccca6e0248bc2fcf8e | 5948001 (PREDICT v2), 5944018 (CapScan; 5948002 v2 control gave the identical file) | not recorded in the producer receipt |

The panel b table the producer read
(`LAB_ROOT/MM_review_analysis/figure_restoration_2026-09-16/illustrator_2026-09-24/vector_redraw/source_data/common_search_participants.tsv`,
sha256 80244f3615800861a915d8698a8a958aa17ce265ea9b84674132924e131a9a15) is the row
concatenation of the two `per_participant.tsv` files above with a `cohort` column (checked:
all 108 rows, all numeric columns equal). Its participant column holds raw PREDICT patient
identifiers and is not used here.

Pixel results of this unit for v30 (notebook, fonts present): composite PNG and PDF are
IDENTICAL to the v30 producer output. The v30 producer read the pseudonymised panel b table
in its repository row order, so the earlier v28d draw-order difference does not arise.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| Supplementary Figure 7, all panels | `word/media/image15.png`, sha256 b88113c0... | BYTE_MATCH; producer WRITE_SITE_HASH_RECEIPT `J/code/build_supp_fig7.py:73`, receipt `J/figure/FIGURE_VALIDATION.json`; other candidate `J/code/build_proposals.py:79` | `J/code/build_supp_fig7.py` writes `J/figure/Supplementary_Figure7_revised.png`; re-executed here from repository source data | yes | Codex did not rerun the producer; this unit did (pixel results above). `build_proposals.py` only proposes text changes and does not draw the figure. Codex lists no rebuilt candidate for this figure in REBUILT_CANDIDATES.tsv; v28f keeps the v28d image. |

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| `identity_comparison.tsv` | `A21/results/identity_comparison.tsv` (producer input); package sheet `PKG/common_search_identity.tsv` | c5a863acb448d698d0a89b440e5c0a42064e55a5475b0de3d793f6f836287090 | copied unchanged. The package sheet (ea77408b...) has the same 464 rows; some float fields differ in the last printed digit (for example 1.0685236588987155 vs ...157), no plotted value changes |
| `identity_per_acquisition.tsv` (9-50) | `FIN/len950/SFig07/build/source_data_len950/identity_per_acquisition.tsv` (min_length 9 rows of `FIN/len950/SFig07/work/hstool/identity_per_acquisition_len.tsv`, `hstool_recount.py`, job 6072769) | 95e1ea4bb104f6cd41f978db5ef790cbb0b71bdfde42d0420ead16be03413b67 | copied unchanged. Replaces the 5-50 table (`J/out/compare/identity_per_acquisition.tsv`, 56c7271c...); the recount's 5-50 rows equal that table (len950 control check). Package sheet: not yet re-exported at 9-50 |
| `identity_summary_by_length.tsv` | `FIN/len950/SFig07/work/hstool/identity_summary_len.tsv` | 47f2bc3f60d1e0a587eddbdf2e043aebede2dd774adfa8c91eef0c9d13720f87 | copied unchanged; its 5-50 rows equal the earlier `identity_summary.tsv` (f4b9f08e...) |
| `completeness_by_length.tsv` | `FIN/len950/SFig07/work/hstool/completeness_len.tsv` | cb7f9080323d48c26254c6398d68db495e575825f81d136df62a6ad59a797a48 | copied unchanged; its 5-50 rows equal the `features` rows of the earlier `completeness_all.tsv` (0d73c378...) |
| `microbial_group_pairs.tsv` | `J/out/compare/microbial_group_pairs.tsv` = `PKG/per_acquisition_confidence/microbial_group_pairs.tsv` | e8d58907fb5564149344d76b09e6dbe8b3dd3062e83084ed14da083143471329 | copied unchanged; byte-identical. `PKG/Supplementary_Figure7cd_pair_rows.tsv` (e96ecb67...) is the producer's re-export of the same rows |
| `microbial_group_pair_summary.tsv` | `J/out/compare/microbial_group_pair_summary.tsv` = `PKG/per_acquisition_confidence/microbial_group_pair_summary.tsv` | 776895eadb4f32753d6521237ba2322ced6fa7831ed83d33218c274cac3caf28 | copied unchanged; byte-identical; every row recomputed by `panels.py:consistency` (all agree) |
| `common_search_participants.tsv` (9-50 peptides) | `FIN/len950/SFig07/build/source_data_len950/common_search_participants.tsv` | 4d53cd6c590d2f354972ba834c087bee23def3e719c3eefe823734fefe2b672c | copied unchanged. Built by `build_sfig07_len950.py` from the pseudonymised table (`PKG/common_search_participants.tsv`, c1cb93b5...) by replacing the peptide columns with the 9-50 recount (array 6072755, aggregated as `summarize_v2.py`); gained/lost peptides emptied; feature and protein-group columns unchanged |
| `common_search_pair_CV.tsv` | `LAB_ROOT/MM_review_analysis/figure_restoration_2026-09-16/illustrator_2026-09-24/vector_redraw/source_data/common_search_pair_CV.tsv` = `PKG/common_search_pair_CV.tsv` | d0fd92b1481b938e3de5967496b80c47b234a976e16524ba4d240020b2c7f1e2 | copied unchanged; byte-identical; used only for the producer's consistency check |
| `technical_agreement_summary.tsv` | `A21/results/leaders/technical_agreement_summary.tsv` (written by `A21/frozen/leaders_revision2.py`, sha256 d3d905d0..., job 5928573) = `PKG/technical_agreement_summary.tsv` | c9d59a24319a1106ab0e1e7895b50133040d78a28faca01918c1895bf7632751 | copied unchanged; byte-identical |
| `arm_definitions.tsv` | `PKG/per_acquisition_confidence/arm_definitions.tsv` | 18e520c0595b21296a832d96bf09330cff2fba3c5c8ba56b63dffe0248e1c93c | copied unchanged |

Not copied: `PKG/Supplementary_Figure7_plotted_values.tsv` (8fab291f...) and
`PKG/Supplementary_Figure7a_per_acquisition_column.tsv` (d2b9cab9...) are producer outputs;
the notebook writes the equivalent `build/SFig07/plotted_values.tsv`.

`SOURCE_DATA_MAP_reanchored.tsv` rows RES-57, RES-57b-e, RES-58, RES-58b-c, RES-59, RES-59b-c,
RES-59s1-s4, RES-60, RES-60b-d, RES-61, RES-61b-j, RES-62, SF7-a, SF7-b and SF7-cd point to the
tables above; every value they record is reproduced in `NUMBERS.tsv`. No map row covers the
5.9% PSM PEP result on line 115.
