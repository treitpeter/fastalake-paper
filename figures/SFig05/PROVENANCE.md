# Provenance: Supplementary Figure 5

`RC` below stands for
`LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24`.

`FINAL` stands for `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`.

## Panels

The v30 pages (first two rows) restyle the v28f/v29 pages (remaining rows) without
changing any plotted value: the v30 build ran the previous producer unchanged, recorded
every axes' plotted arrays, labels, limits and ticks, and found the restyled figure
identical (`FINAL/v30/figures_style/SFig05/CHANGES.md`, `build_report.json`).

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a-f (v30) | `FINAL/v30/figures_style/SFig05/build_sfig05_style.py`, sha256 `4d03dcc2b0b7508eb83009e495d41f07f187636dcf84be4cbbb307cee4a663c5`, helper `FINAL/v30/figures_style/_sf456_common.py` (sha256 `c2a03a2c...`), run by `run_sf456.sh` | `FINAL/v30/figures_style/SFig05/SupplementaryFigure5.pdf` and `.png` | pdf `ffad3ec1135bcceb414a06d5ab23f92557b5385c56c65fb39290d0e887a352d5`; png `cc3d154ddd9ad10cdd041cbe88052740737ffe0774a7189328fa788e50f4bcf4` | not recorded | CPython 3.12.4, matplotlib 3.10.7, PyMuPDF 1.27.1 (repository environment) |
| continued (v30) | as a-f (v30) | `FINAL/v30/figures_style/SFig05/SupplementaryFigure5_continued.pdf` and `.png` | pdf `1189ae87bbf80921c45884ed4894ed7a566dd369c33f2e8a5473edf2ce13f078`; png `5cfa2801bed4b38229081f391fbe447c8f255ea6eeb84a4fd107d7d24fcc89ff` | not recorded | as above |
| a-f, step 1 (draw, v28f) | `RC/supplements/revision_2026-09-24/build_revision.py`, error-control block (lines 89-119; `save()` lines 32-50), sha256 `e636bd85164896f53ec1d13ef671754ff800597d2d955a3eb01f5151dced1bbe` | `RC/supplements/revision_2026-09-24/delivery/figures/ERROR_CONTROL_ARCHITECTURE.pdf` (v28f; no longer in `expected/`) | `a7da2f7f0585184c4028e2e5a9f373eef6214082b3f883803f5505d807ec2bf8` | not recorded (run outside Slurm) | matplotlib 3.10.8 (PDF metadata); DejaVu Sans |
| a-f, step 2 (title removal, 300 dpi) | `RC/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/prepare_supplement_figures.py` (key `4` of `SUPPLEMENT_FIGURE_SOURCES.json`), sha256 `5c52b4b89037916a4503e3b77e8454522186a83c0fdff48cced5f1776141ce4d` | `RC/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/release_candidate/supplementary_figures/Supplementary_Figure_04.png` and `.pdf` (v28f) | png `3d5c34c63457a1ea29dae8bb43fcd9002f62c474e34d69a3ee407b07e70b5eb8`; pdf `30d281cb5da7dd19ce5db805eca976c6410995039dd7cfdbaae04c0cc405f95d` | not recorded | PyMuPDF (version not recorded) |
| continued | `RC/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/package_release.py` lines 14-19, sha256 `17ff4cbbe7384650d2300f26b69b5d9798fa1010d5b5e20722815785592d0500` | `RC/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/release_candidate/supplementary_additions/Supplementary_Figure_4_decomposition_proposed.png` and `.pdf` (v28f) | png `e139bbe288933b63cdbc3bc04c6b7cd45d6b53d4d42528a821c10b41202c4c4e`; pdf `ca751f51588fab1a77106219b3afbe8768026c4d5cd4250acf7dba4cbe36ac2e` | not recorded | matplotlib 3.10.8 (PNG metadata); Söhne Buch TrueType (`soehne-buch.ttf`) |

The figure was numbered Supplementary Figure 4 in these producer runs; it is
Supplementary Figure 5 in v28d.

Regeneration in this repository (`panels.figure_a_to_f`, `panels.figure_continued`;
CPython 3.12.4, matplotlib 3.10.7, PyMuPDF 1.27.1, pandas 2.3.3, numpy 2.2.6, scipy
1.16.3), headless notebook build on Slurm (2 cores, 16 GB requested, 4 s). Receipts with
input and output hashes are in `checks/`. Results:

| page | pixel check | note |
|---|---|---|
| a-f | IDENTICAL, 0 differing pixels against `expected/SFig05_v30.png` | the repository PDF also rasterises identically to the v30 PDF at 300 dpi |
| continued | IDENTICAL, 0 differing pixels against `expected/SFig05_continued_v30.png` | as above |

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| a-f | Supplementary Figure 5, `word/media/image12.png`, sha256 `3d5c34c6...` | BYTE_MATCH; producer COPY_ONLY (`v15_MM_coauthor_2026-09-25/build_figures_v15.py:36` only copies); "original plotting/rendering producer unresolved" | drawn by `build_revision.py` (error-control block), title removed and rasterised by `prepare_supplement_figures.py`; confirmed by regenerating a byte-identical PNG from source data | agree on bytes; this unit resolves the producer Codex left open | the copy chain Codex found (`coauthor_package`, `embedded_media/supp_04.png`) ends at the `release_candidate` PNG written by `prepare_supplement_figures.py` |
| continued | Supplementary Figure 5, `word/media/image13.png`, sha256 `e139bbe2...` | BYTE_MATCH; WRITE_SITE `package_release.py:19`; "producer was not rerun" | `package_release.py` lines 14-19; rerun here, pixel-identical | agree | rerun closes Codex's "not rerun" limit |

`REBUILT_CANDIDATES.tsv` lists no candidate for Supplementary Figure 5, and the figure is
not among the v28f rebuilt figures: v28f uses the v28d images unchanged.

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| `fusion_pooled_current.tsv` | package sheet `Supplementary_Figure_05/fusion_pooled_current.tsv`; written by `build_revision.py` from `RC/supplements/revision_2026-09-24/source_data/original_S24/08af3135cab4_fusion_summary.tsv` (human cohorts, CAMPI_MOCK rows dropped) and `.../current_CAMPI_fusion/pooled_by_seed.tsv` | `0f90dc0a3003ae3f5097832ec9c5cd55404647d065593adda72d04bc77ba121e` (equal to `RC/supplements/revision_2026-09-24/delivery/source_data/fusion_pooled_current.tsv`) | copy; values agree with `panels.stat_ac_fusion_pooled` to 1e-16 |
| `fusion_acquisition_components.tsv` | package sheet `Supplementary_Figure_05/fusion_acquisition_components.tsv`; written by `build_revision.py` from `.../current_CAMPI_fusion/independent_per_acquisition.tsv` and `.../original_S24/2d33c301d1a1_fusion_recount.tsv` | `84243602ce1bdac6cf5270e110a6be9d2e9d751ce95525ef89eec8ba0bd5961c` | copy |
| `fusion_decomposition_current375.tsv` | `RC/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/FUSION_RECOMPUTED.tsv`, written by `verify_evidence.py` (sha256 `85dfaba502066b6e404e1d8812ae6a04b9418edabfd65045fe675ee2eac69df2`) from `LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/campi_fusion_rerun/results/independent_per_acquisition.tsv` | `1b13ec931e8b0bdcd18e07a1c8341327662ec323688853682ae7964159e2374c` | copy of the producer's input; package sheet (`2a134617...`) agrees to ≤ 2e-15; `panels.stat_cont_fusion_decomposition` reproduces it from the components (max difference 1.8e-15) |
| `donor_entrapment_seed_means.tsv` | package sheet `Supplementary_Figure_05/donor_entrapment_seed_means.tsv`; = `RC/supplements/revision_2026-09-24/source_data/original_S24/5c6b0af7c722_donor_seed_means.tsv` re-exported | `f1b4658e64b3a0e4b854ad38fb5b38145b30f6ec0438d6891bdaa0703a43581a` | copy; `panels.stat_d_donor_seed_means` reproduces it (max difference 7e-12 in counts, 1e-16 in fractions) |
| `donor_recount.tsv` | `RC/supplements/revision_2026-09-24/source_data/original_S24/4e51dd75c37d_donor_recount.tsv` (A11 `results/donor_recount.tsv`) | `4e51dd75c37deed3f13c84ad2407032b01ae16d877efac9d1b60a174aafb5ea8` | leading `index` column removed with `cut -f2-` |
| `simulation_mixture_outcomes.tsv` | package sheet `Supplementary_Figure_05/simulation_mixture_outcomes.tsv`; = `.../original_S06/3027bfcc1c6e_ordinary_outcomes.tsv` re-exported | `783e4aa37060ded3d937be84495cfdb2c38d08fcfb07411b1ca963c9722c8387` | copy |
| `simulation_all_outcomes.tsv` | `LAB_ROOT/analyses/fdr_reference_benchmark_2026-09-10/runs/validation/collected/primary.tsv` | `e6427d549e37a18d60d26447c962c274f08c2b41a58594a66edbadc7029067a1` | copy |
| `simulation_uncertainty.tsv` | package sheet `Supplementary_Figure_05/simulation_uncertainty.tsv`; = `LAB_ROOT/analyses/fdr_reference_benchmark_2026-09-10/runs/validation/collected/uncertainty.tsv` (sha256 `7492d210...`) re-exported | `f267eee37eaa0cf29e122b82a27914abdb309a8edc31a8d4fd9c397d20b0e0ec` | copy; `panels.stat_ef_simulation_uncertainty` reproduces it to 1e-16 with the original seed |
| `simulation_design.tsv` | `LAB_ROOT/analyses/fdr_reference_benchmark_2026-09-10/common.py` line 34 (sha256 `8899096b...`), `generate.py` lines 47 and 60 (sha256 `f317c410...`), `LAB_ROOT/analyses/fdr_validation_2026-09-08/results/campi_reference_v2/REF_targets.fasta` (sha256 `1dccc896...`, 29,673 records) | per row | written by hand from those lines; FASTA records counted with `grep -c '^>'` |

`SOURCE_DATA_MAP_reanchored.tsv` rows SF5-ac, SF5-d, SF5-ef, SF5-cont and RES-42 point to
these package sheets; all agree with the tables here. No n = 15 replacement table exists
for Supplementary Figure 5.

## Upstream code (not run here)

| analysis | code (LAB_ROOT/...) | sha256 |
|---|---|---|
| healthy_stool_donor_entrapment | `MM_review_analysis/A11_4_2_control_reconciliation/frozen/analysis.py` (copied as `RC/supplements/revision_2026-09-24/code/evidence/CONTROL.py`, lines 72-79 compute seed means and fusion summaries) | CONTROL.py `51ed28a4b419c4c1963cc36ba6e1c37218f460f73bceeee6bb10868aa2c7b3d5` |
| campi_fusion_375m | `00_PACKAGE_2026-09-17/cleanup_2026-09-22/campi_fusion_rerun/audit_fusion.py`; Slurm jobs 5984954-5984958 (`status/SUBMITTED.json`) | `a5289adc1cdd415e7529c51b27cdb0d8fd48710c0139b68bc7a844f36991433d` |
| restricted_reference_simulation | `analyses/fdr_reference_benchmark_2026-09-10/collect.py` (bootstrap), `evaluate.py`, `generate.py`, `common.py`; original interpreter per `job.sbatch` CPython 3.12.4, numpy 2.2.6, scipy 1.16.3 | collect.py `bd884e5ee02caa94c8913e3f37908d1b56eb8a2863899de2084bef6ae2c86198` |
