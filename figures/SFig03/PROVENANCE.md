# Provenance: Supplementary Figure 3

`FIN` = `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`.
`I5` = `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration5`.
`REV` = `LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/supplements/revision_2026-09-24`.
`SDPKG` = `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Figure_03`.

## Panels

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a-d (one figure) | `I5/figures/build/sf3.py` (sha256 170af12a5088afa9dda2698e3aeaf73ea32f9ba77289be88cba1cda1a1fe3357) with `I5/figures/build/common.py` (d4fb7de9b65e7db67bc1f5caf0cbb28c2e95f55aaaad6d6b1c51ff56a64ba7fe), launcher `I5/figures/build/run.sbatch` | `I5/figures/out/SupplementaryFigure3_R5_400dpi.png`; vector `I5/figures/out/SupplementaryFigure3_R5.pdf`; receipt `I5/figures/check/SF3_BUILD.json` | PNG cb1411f832b003a1fa441800964ddf11792b6728b68b4b469f6a6c7765cdd2bd; PDF 8b7a2b38b9c854b7debe1d309af3772f50a161e752a9119473bd44a76f57c7fd | 6052419, 6052430, 6052434 (each wrote the same PNG hash; 6052434 is the last write, `I5/figures/logs/`) | CPython 3.12.4, matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3, PyMuPDF 1.27.1; Söhne Buch TTF sha256 db966efe63405323b2d1c2cfb3f104491d61ad56bf3a95fe3bfa712965e40ed8 |

| a-d v30 (one figure) | `FIN/len950/SFig03/build_sfig03_len950.py` (sha256 eedb5c068da6be09c81bab0cb1e066f9ca0ccbfd01fa7a8737a34a22e32a0cd4) with `FIN/len950/SFig03/panels_len950.py` (124613e9688c444eac5819e2342984407532772198b444e7db92b24787827d32), then `FIN/v30/figures_style/SFig03/build_sfig03_style.py` (efe16c6ab5f5288f06f9d25f7031e2f7c96df83e91c023b3eb7dd8b8de182696) with `FIN/v30/figures_style/_shared_sf01_03/restyle_common.py` (0c7ced18eb696d77db2639e5121fb66399debe60684097d15cd59cff816addb1) | `FIN/len950/SFig03/SupplementaryFigure3_len950.pdf`; `FIN/v30/figures_style/SFig03/SupplementaryFigure3.pdf` and `.png`; receipts `values_check.json`, `pdf_audit.json` | len950 PDF 8976aefa66807bc0073ef973f959e14619f7a07de82105b33b0f2231c7336fd9; v30 PDF 215d8d30e9bd0315d91c82112a6eea8cafd4bfb84d60f43566e1c2719ee4a1a3; v30 PNG ebc18430635566789888c1ec031579ca0148aa48faf6953bdc8337d67ccaefd5 | not recorded (run with `srun` via `FIN/v30/figures_style/_shared_sf01_03/run.sh`) | paper environment (`fastalake-paper_v28_work/env`), matplotlib 3.10.7, PyMuPDF 1.27.1 |
sf3.py is a redraw at placed size of the CAMPI_RECOVERY block of `REV/build_revision.py`
(sha256 e636bd85164896f53ec1d13ef671754ff800597d2d955a3eb01f5151dced1bbe): same inputs,
marks, colours and jitter, with the approved label edits FAR15 and FAR16 written natively.

Reproduction here: `panels.py` holds the sf3.py drawing code and the common.py helpers
(init, slant, render_png, min_text_pt) unchanged apart from repository-relative paths and
the font location. Run on Slurm (job 6069001 failed on a crop-box bound in the new
per-panel PDF cut, fixed; the final build ran as job 6069043, 12 s wall time, peak RSS 0.82 GB; the same job also ran the notebook without fonts, which passes and skips the pixel checks).
Result: composite PNG and composite PDF at 300 dpi IDENTICAL to the producer output, and
each panel crop IDENTICAL (`checks/`). The producer's plotted-values table
(`I5/figures/plotted_values/SupplementaryFigure3_plotted_values.tsv`, sha256
5019ea9096173bbeec37b25fce9fd898aa1a1733e5397552a0d7d87467b80614) is in
`expected/plotted_values.tsv` with its `source` column reduced to file names (lab paths
removed; no other change); the notebook's plotted values equal it exactly.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| Supplementary Figure 3 (a-d) | Supplementary Figure 3, `word/media/image10.png`, sha256 cb1411f8...cdd2bd | origin BYTE_MATCH to `I5/figures/out/SupplementaryFigure3_R5_400dpi.png`; producer WRITE_SITE_HASH_RECEIPT at `sf3.py:68` with receipt `SF3_BUILD.json`; limit "producer was not rerun" | same file and producer; producer now re-run from repository source_data with pixel-identical output | yes | Codex REBUILT_CANDIDATES.tsv has no Supplementary Figure 3 row: no competing rebuild exists. The census also lists byte copies in docx working folders; they are copies, not producers. |

Codex number recheck (`out_recheck/REPORT.md`, second-engine section) recounted the frozen
CAMPI Percolator files independently: 15 acquisitions per arm, median FL/MG 0.936656 and
FL/REF 0.986541. This unit gives 0.9366555214583429 and 0.9865405953012688: agree.
The Codex statistics register has no row for v28d lines 83, 89, 253, 417 or 685 statistics
of this figure (descriptive medians and counts only); no test is involved.
`SOURCE_DATA_MAP_reanchored.tsv` rows RES-40 (0.986541, VERIFIED), RES-40b (84.6%,
LOCATED only) and SF3 (VERIFIED): RES-40b is now computed by `panels.py:stat_line89_msfragger`.

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| CAMPI_member_support_recount.tsv | `SDPKG/CAMPI_member_support_recount.tsv` = `REV/delivery/source_data/CAMPI_member_support_recount.tsv` (written by `REV/build_revision.py` member recount from the membership tables) | a230cf25e366bf279891b80875be88ddc8a7838673e4faac74fa8c61a9d6d107 | byte copy; package sheet agrees (identical). Also equals `LAB_ROOT/MM_review_analysis/A14_4_3_campi_gap/results/member_summary.tsv` (pooled counts, acquisitions_ge1), and is recounted exactly by `panels.py:stat_recount_membership` |
| CAMPI_member_acquisition_counts.tsv | `SDPKG/CAMPI_member_acquisition_counts.tsv` = REV delivery copy | 4e995436231bf70ab90cfd44e68b56b3b3d0757f015a21967470b96eca138486 | byte copy; recounted exactly |
| CAMPI_reference_recovery.tsv | `SDPKG/CAMPI_reference_recovery.tsv` = REV delivery copy, a pandas re-export of `LAB_ROOT/MM_review_analysis/A14_4_3_campi_gap/results/paired_reference_recovery.tsv` (c446dabab7ad3c39c0d055424f60e19c50411eaebe82ee06fae52c54a4f74919, values equal) | 7f2d7a7586c57b47332ba62d7348aceac33aea79dbed3218c46d2dc9d4cb8551 | byte copy; recounted exactly (ratios to 1e-12) |
| CAMPI_within_MSFragger.tsv | `SDPKG/CAMPI_within_MSFragger.tsv` = REV delivery copy, a pandas re-export of `LAB_ROOT/MM_review_analysis/A07_4_2_second_engine/results/within_MSFragger_comparisons.tsv` (a89e4461af2d44d24efedb006c288bcb41b2c6c57ab9bb685d2fad8eec95fcdd; bytes differ in float formatting only, values equal) | 13490970bb7873ff2dd349f8de1256e656329bc2c8433fc448702c8c98ae5c47 | byte copy; every count checked against CAMPI_MSFragger_search_receipts.tsv |
| second_engine_between_engine_identities.tsv | `SDPKG/second_engine_between_engine_identities.tsv` = `LAB_ROOT/MM_review_analysis/A07_4_2_second_engine/results/between_engine_identities.tsv` | d455994adac00eb735a6955bd13ab265de2ca2bdbc2c04d5f42a4dd6c6274501 | byte copy; Source Data only |
| CAMPI_MSFragger_search_receipts.tsv | `LAB_ROOT/MM_review_analysis/A07_4_2_second_engine/results/<sample>__<arm>/COMPLETE.json` (45 files; not a package sheet) | per file, not listed; e.g. S01__FL_UNIPROT e7b8ffe1b9f3d57595645781711bce0ed6243e79ae6152fedca7dbe584e21d37 | fields sample, arm, status, confidence, accepted counts, target entries, raw and FASTA sha256, job id copied without paths (one-off export in the work area) |
| CAMPI_membership_<sample>.tsv (15) | `REV/delivery/source_data/CAMPI_membership/<sample>.tsv` = `LAB_ROOT/MM_review_analysis/A14_4_3_campi_gap/results/<sample>/peptide_membership.tsv` (not a package sheet) | listed in `source_data/SFig03/README.md` (identical) | byte copies, renamed with a `CAMPI_membership_` prefix |

The Source Data package for Supplementary Figure 3 holds the first five tables; all five
agree byte for byte. `SOURCE_DATA_MAP_reanchored.tsv` row SF3 names the same five files.


## v30 update (2026-10-06)

The v30 figure supersedes the v28f row above. `panels.py` now sets `PLOTTED_MIN_LENGTH_CD = 9`
(the only change of `panels_len950.py` against the earlier `panels.py`), and holds the v30 redraw
(`build_figure_v30`, `build_v30`) and the restyle helpers (`fingerprint`, `compare`, `texts`,
`align_letters`) moved unchanged from `build_sfig03_style.py` and `restyle_common.py` apart from
repository-relative paths and `fl_style` from `style/`. Run on Slurm (job 6075286, notebook via
`tools/run_notebooks.py`, 25 s): the 9-50 producer PDF and the v30 PNG are IDENTICAL to
`expected/` (`checks/SFig03_len950.json`, `checks/SFig03_v30.json`). The v28f expected images,
panel crops and their checks were removed. Source data are unchanged.
