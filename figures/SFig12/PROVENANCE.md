# Provenance: Supplementary Figure 12

`LAB_ROOT/.../` abbreviates
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/`.
`REV24/` abbreviates
`LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/supplements/revision_2026-09-24/`.
`A33/` abbreviates `LAB_ROOT/MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/`.
`C2/` abbreviates `LAB_ROOT/.../coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig5e_f/`.

## Panels

All six panels are drawn by one producer run into one figure.

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| v30 figure (a-f, one canvas) | `LAB_ROOT/.../final/v30/figures/SFig12/build_sfig12.py` (sha256 5fc81f1c74dc6f772ae7d50fc0621de3c67c077e09a72bad2d04d19c2b622902) | `LAB_ROOT/.../final/v30/figures/SFig12/SupplementaryFigure12.png`; `.pdf`; `plotted_summary.tsv` | png e644ae023b54c93e6926e4041b7dfdfcb2fbd566f6d3b02a0d02e48f3a6c211c; pdf fc58c6f14d523a36deefc4a15b11207748fc0c17eba6cfada064ebf8ede13c35; summary 548d5790055d8754f15f15d050ef4520802abe2f32000968e0e84ac9f910681b | not recorded | fastalake-paper_v28_work env (Python 3.12), matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3; Söhne via fl_style |
| a-f (v28f figure; superseded) | `REV24/build_revision.py` (sha256 e636bd85...), sliced verbatim into `LAB_ROOT/.../final/suppfigs_n15/SF12/builder/sf12_build.py` (sha256 79e12e74...), run by `.../SF12/alloc347/build_n15a.sh` (sha256 d0231fdd...) | `.../SF12/alloc347/build_n15a/figures/FUNCTIONAL_ATTRIBUTION.pdf` | fdf84699794903254cda6c3c360e4e1b05fd14fa99a0ac8013dc91b85b4b84be | 6066033 (post347: check and rebuild) | CPython 3.12.4, matplotlib 3.10.7, pandas 2.3.3, numpy 2.2.6 |
| v28f composite (suptitle removed, 300 dpi; superseded) | `LAB_ROOT/.../final/suppfigs_n15/SF12/export.py` (sha256 e1538b19...) | `.../SF12/alloc347/out/SupplementaryFigure12_n15a.pdf`; `..._n15a_300dpi.png` | f9657339c91bfaf61b71d25d9a7c2bd77fea32453f12c76cdba60a5c07e4a4a4; 0ff82388b06d441badd4ddfb7242a7cd65bb218c91db511dcf7f9ef0a117a825 | 6066033 | as above, PyMuPDF 1.27.1, Pillow 12.0.0 |
| a, input matrix CAPSCAN | `A33/revised/build_variants.py` (byte copy `.../SF12/alloc347/producer/build_variants.py`, sha256 57ea24e6...), runner `.../SF12/alloc347/stage.sh` on the C2 grouping run of the 347 roster | `.../SF12/alloc347/run347/CAPSCAN_KO_all.tsv.gz` | 69a953d61e60764bef720835bce495ab0387c114a7ea39c3d78e517dd7ff9a86 | 6066031 (run347); control 6065988 (repro331 reproduced the 331-roster A33 matrices) | CPython 3.12.4 |
| a, input matrix PREDICT | `A33/revised/build_variants.py` via `A33/jobs/repaired_predict.sbatch` | `A33/inputs/matrices_repaired/PREDICT_KO_all.tsv.gz` (= `REV24/source_data/allocation/PREDICT_KO_all.tsv.gz`) | f452039de3d0e8b07ec6f3ef7116a29104becb771ae2405d465a55d5ab5268b8 | 5973734 | CPython 3.12.4 |
| a, input matrix HSTOOL | `A33/revised/build_variants.py` | `A33/inputs/matrices/HSTOOL_KO_all.tsv.gz` (= `REV24/source_data/allocation/HSTOOL_KO_all.tsv.gz`) | c3533c1b8846f487a4a109fd8344b4b36f0ee53f2164b4d02a3c5f80a104791a | not recorded (OPEN_ITEMS.md) | CPython 3.12.4 |
| c, d, f inputs | `C2/scripts/fig5ef_run.py` | `C2/run347/CAPSCAN__{simplified,original}/source_data/Figure5/` | see Source data below | 6058098 (run347); control 6058096 (control331 reproduced the 331 contributor table byte for byte) | CPython 3.12.4 |
| f, summary table | `LAB_ROOT/.../final/suppfigs_n15/SF12/builder/s10c_347.py` (sha256 11dd1164...) calling `four_rules` of `.../SF12/builder/check_s10c_eligibility.orig.py` (sha256 db84d6e6...) | `.../SF12/inputs_n15/s10c_347/s10c_rule_comparison_347.tsv` | 8dae7133ceec728595b340693b31c6e6766ce0cb1ef82dea387dc0c3fd1aa088 | not recorded (log `.../SF12/logs/s10c_347.log`) | CPython 3.12.4 |
| b | `A33/revised/derive_hstool_pc.py` (function_tier_summary.tsv) | `REV24/source_data/original_S09/8b8c865efed0_function_tier_summary.tsv` | 93ed4ead5ea85f8d4fd8bf09750bb07dc8a66c0fe8565129f20c85988f496147 | not recorded | CPython 3.12.4 |
| e | A33 S10 summary (render lineage `A33/render/Supp_9_10_19/build_supp_9_10_19.py`) | `REV24/source_data/original_S10/945c8cc37ffd_SUMMARY.tsv` | fb4e1e3c34a233f627aad13fd78e33af77c74d9b785959ce7f934cbd8d589716 | not recorded | CPython 3.12.4 |

The repository environment has the same versions as the producer interpreter (CPython 3.12.4,
matplotlib 3.10.7, pandas 2.3.3, numpy 2.2.6, PyMuPDF 1.27.1, Pillow 12.0.0). The notebook
render is byte-identical to the producer's 300 dpi composite and pixel-identical at 300 dpi to both
producer PDFs and to each panel cut from the composite (`checks/*.json`, 0 differing pixels).

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| v28d Supplementary Figure 12 (docx `word/media/image20.png`, sha256 1f8fffe3...) | Supplementary Figure 12, `word/media/image20.png` | BYTE_MATCH to copies in the 25 Sept coauthor package; producer_status COPY_ONLY; pixel-identical to `.../SF12/out/SupplementaryFigure12_old_300dpi.png`; "original plotting/rendering producer unresolved" | drawn by `REV24/build_revision.py` (FUNCTIONAL_ATTRIBUTION block) on the 24 Sept inputs; suptitle removed and rasterised at 300 dpi (reconstructed by `.../SF12/export.py`) | partly | same file identity; Codex's pixel match to `_old_300dpi.png` is the evidence for this unit's attribution. The 24 Sept producer is identified here; the tool that removed the suptitle for v28d is not, so the post-step remains a reconstruction (OPEN_ITEMS.md) |
| v28f candidate PNG (sha256 0ff82388...) | REBUILT_CANDIDATES.tsv, `.../SF12/alloc347/out/SupplementaryFigure12_n15a_300dpi.png` | WRITE_ROUTE_INSPECTED_NOT_RERUN, `export.py:14`; differs from all v28d media | the same file; rerun here from repository code, byte-identical | yes | this repository reruns the write route |
| v28f candidate PDF (sha256 f9657339...) | REBUILT_CANDIDATES.tsv, `.../SF12/alloc347/out/SupplementaryFigure12_n15a.pdf` | WRITE_ROUTE_INSPECTED_NOT_RERUN | the same file; rerun here, pixel-identical at 300 dpi | yes | |
| superseded candidate (panel a still on 331) | REBUILT_CANDIDATES.tsv, `.../SF12/out/SupplementaryFigure12_n15_300dpi.png` (sha256 107b0d6d...) | rebuilt candidate | not used; panel a is not on the 347 roster | yes | the decision names `alloc347/out` as the target |

## Source data

The two original-grouping cell tables (`capscan_contributor_cells_original.tsv`, `capscan_contributor_region_cells_original.tsv`) were removed in v30 together with the original series of panel f.

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| `source_data/SFig12/HSTOOL_KO_all.tsv` | `REV24/source_data/allocation/HSTOOL_KO_all.tsv.gz` | c3533c1b8846f487a4a109fd8344b4b36f0ee53f2164b4d02a3c5f80a104791a | gunzip only; no package sheet (matrix underlying package sheet `Supplementary_Figure_12/KO_allocation_coverage.tsv`) |
| `source_data/SFig12/PREDICT_KO_all.tsv` | `REV24/source_data/allocation/PREDICT_KO_all.tsv.gz` | f452039de3d0e8b07ec6f3ef7116a29104becb771ae2405d465a55d5ab5268b8 | gunzip only |
| `source_data/SFig12/CAPSCAN_KO_all.tsv` | `.../SF12/alloc347/run347/CAPSCAN_KO_all.tsv.gz` | 69a953d61e60764bef720835bce495ab0387c114a7ea39c3d78e517dd7ff9a86 | gunzip only |
| `source_data/SFig12/KO_allocation_coverage.tsv` | `.../SF12/alloc347/build_n15a/source_data/KO_allocation_coverage.tsv` | 05a9a9d09028d3835706be3fad3d4f71903c59288cba70056fd92744b1648b01 | copy; re-derived byte for byte by `panels.py:coverage_from_matrices` + `write_coverage`. Package sheet `Supplementary_Figure_12/KO_allocation_coverage.tsv` (sha256 95223a27..., 330 CapScan rows) is the 331-roster version: HSTOOL and PREDICT rows agree, CapScan rows differ (346 rows, shared acquisitions within 0.51 percentage points); the package sheet must be replaced |
| `source_data/SFig12/healthy_stool_leader_diversity.tsv` | `.../SF12/alloc347/build_n15a/source_data/healthy_stool_leader_diversity.tsv` (export of `REV24/source_data/original_S09/8b8c865efed0_function_tier_summary.tsv`) | 93ed4ead5ea85f8d4fd8bf09750bb07dc8a66c0fe8565129f20c85988f496147 | copy; identical to package sheet `Supplementary_Figure_12/healthy_stool_leader_diversity.tsv` (same sha256) |
| `source_data/SFig12/healthy_stool_grouping_leader_sensitivity.tsv` | `.../SF12/alloc347/build_n15a/source_data/healthy_stool_grouping_leader_sensitivity.tsv` (export of `REV24/source_data/original_S10/945c8cc37ffd_SUMMARY.tsv`) | fb4e1e3c34a233f627aad13fd78e33af77c74d9b785959ce7f934cbd8d589716 | copy; identical to package sheet `Supplementary_Figure_12/healthy_stool_grouping_leader_sensitivity.tsv` (same sha256) |
| `source_data/SFig12/CapScan_grouping_leader_sensitivity.tsv` | `.../SF12/alloc347/build_n15a/source_data/CapScan_grouping_leader_sensitivity.tsv` (A33_ rows of `.../SF12/inputs_n15/s10c_347/s10c_rule_comparison_347.tsv`) | 593d0251c4db0947172fd0437b342908546ea1910cac6cdfa18dd801fe5af9b2 | copy; recomputed by `panels.py:stat_f_recompute` (agreement within 5.6e-17, counts equal). Package sheet `Supplementary_Figure_12/CapScan_grouping_leader_sensitivity.tsv` (sha256 61f61928...) is the 331-roster version (14 donors, 56 donor-regions); must be replaced |
| `source_data/SFig12/capscan_contributor_region_cells_direct.tsv` | `C2/run347/CAPSCAN__simplified/source_data/Figure5/contributor_region_cells.tsv` | e7cb87a78397858a6487ab4c81475f3dfdc7b465cf9fe8643963af4f6e6ee78c | copy; no package sheet |
| `source_data/Fig5/contributor_cells.tsv` (owned by Fig5) | `C2/run347/CAPSCAN__simplified/source_data/Figure5/contributor_cells.tsv` | 6a704174763ce561aa2189fba7b137f3dc480edeebf7f52bdd16cabbb492dd84 | read, not copied. Package sheet `Supplementary_Figure_12/CapScan_function_contributors.tsv` (sha256 d08a7052..., 20,352 cells) is the 331-roster version; the v28f producer's re-export of this table (`build_n15a/source_data/CapScan_function_contributors.tsv`, sha256 cc28f7db...) has the same rows written by pandas with shortest round-trip floats |

Source-data map (`SOURCE_DATA_MAP_reanchored.tsv`) rows for this figure: RES-87 (still 29.4%,
recomputed 29.36 on the 331 table) and SF12-cd (still 20,352) are anchored to the 331-roster
table; the v28d text and this repository give 29.1% and 21,168. RES-85/85b (Fig. 5f, cited with
Supplementary Fig. 12) still record 41.0/44.0 and 56.7/58.2; v28d prints 42.5/45.6 and 58.3/60.0,
which this repository recomputes. SF12-a, SF12-b and SF12-ef locate the files without values.
