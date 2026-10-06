# Provenance: Supplementary Figure 11

`SF11` = `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/suppfigs_n15/SF11`.
`BR` = `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit`.

The producer is `SF11/setup_and_build.sh n15` (sha256 c68090765bc122bc9b22683b1b122ad2b714f8dd922ecf735a2a043e33832a25).
It runs `SF11/builder/build_figures.py` (34971b8e...; = `BR/claude_review_2026-10-03/iteration4/figures/build/d12_build_figures.py`
716b36ee... plus `SF11/builder/SF11_BUILDER.diff` 1f5b82e3..., which changes only the input path and
the CapScan title), `SF11/builder/build_supplement.py` (47cd8d79...; = `d12_build_supplement.py`
b51bf4d1... plus the same diff: input paths and footnote) and `SF11/export.py` (d57efa19...). The
same script with the old inputs (`setup_and_build.sh old`) reproduces the v28d embedded image byte for
byte (sha256 fa6b628a7a1f6f1436e5122ee3d8368d7730115553bed42fc56ba9d4dd543e23; `SF11/logs/build_old.log`).

## Panels

Rows below the v30 row describe the superseded v28f build (panels a-f in four producer figures); they are kept as history.

`FINAL/` stands for `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/`.

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| v30 figure (a-f, one canvas) | `FINAL/v30/figures/SFig11/build_sfig11.py` (sha256 b283403ffaa26ea9a6d0f7b5f530e12304f7cfc01418097e0d80e334c08ebd82) | `FINAL/v30/figures/SFig11/SupplementaryFigure11.png`; `.pdf`; `plotted_values.tsv` | png cb1e56e0c91cf3f2512a519525581f2bce4dedd4e5bb76af3720d438495b1438; pdf d03ac4bc92f79b3e3e73662eb353d9bbc80aa83222369f9681564d6132015df2; values 2597c033404714843f4266b5f731bda03d2b75adb6739b0ccf318129e3ded999 | not recorded | fastalake-paper_v28_work env (Python 3.12), matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3; Söhne via fl_style |
| a | `SF11/builder/build_figures.py` (`pyramid`) | `SF11/shadow_n15/B/figures/Supplementary_crossomics_PREDICT.pdf` | 332f4d537a228a3caac8c46ba493def2407c0df947ffde2370cc0b707388653c | not recorded (see OPEN_ITEMS.md) | Python 3.12.4, matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3, PyMuPDF 1.27.1; font soehne-buch.ttf |
| b | `SF11/builder/build_figures.py` (`pyramid`) | `SF11/shadow_n15/B/figures/Supplementary_crossomics_CAPSCAN.pdf` | 94fe00b85009e6e727a9f3a71f21848e0bfcc12a5ca26299c66ac676ccf4f1af | not recorded | as a |
| c, d | `SF11/builder/build_supplement.py` | `SF11/shadow_n15/B/figures/supp_native_prevalence.pdf` | a0e6d1a4ea23ccbcf19e267d60d239b1c0fd813cd35a4bf68ce713e4c10ff71a | not recorded | as a |
| e, f | `SF11/builder/build_supplement.py` | `SF11/shadow_n15/B/figures/supp_native_methods.pdf` | 9e07e4bf20f44ed7149373cd74f769f71a00e8a21cdc35688a900b7ef78c2fbc | not recorded | as a |
| composite PDF | `SF11/builder/build_supplement.py` | `SF11/out/SupplementaryFigure11_n15.pdf` | 4294fcf6a82e04c6bdbd505cc70e7979f88b451a4b34c36556307b77bc3347a8 | not recorded | as a |
| composite PNG | `SF11/export.py` | `SF11/out/SupplementaryFigure11_n15_400dpi.png` | 9b0504c56345884d6fe5bb6f6ff664f394a9e002e6d9e3f7046e42a326222db2 | not recorded | as a; Pillow 12.0.0 |

Repository result: every panel PDF and the composite PDF rasterise to IDENTICAL pixels at 300 dpi,
and the regenerated composite PNG has the same sha256 as the producer's (`checks/*.json`). The
plotted-value tables written by the notebook equal the producer's
(`expected/all_plotted_crossomics_values.tsv`, SF11 rows; `expected/supplementary_plotted_values.tsv`,
byte-identical).

Upstream statistics (all recomputed in the notebook from the distance matrices; maximum absolute
difference to the producer tables 1.1e-16 for rho and q, 0 for every P value, 0 changes in q<0.05
classification over 252 tests):

| analysis | code (LAB_ROOT/...) | sha256 | Slurm job | outputs used |
|---|---|---|---|---|
| A84 PREDICT functional and native | `BR/crossomics.py` with `check_prevalence.py` (b05f65da...) and `native_pyramids.py` (64ba668b...) | 5d0322e7f1fbe3c4f73001a2840d6fae1572035fe69d97dbb8f9072342e86404 | 6048399 (verified50 profiles); production native outputs 6044312 | `BR/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/cross_refit/{functional,native}_correlations.tsv`, `*_distances.npz` |
| A82 CapScan functional and native, 347 | `BR/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig5g/crossomics_capscan.py` | 0ab9a2465b4308b145dcc75d81423601a5fc50ebe4ef16f5aa9e02b79f9591d2 | 6063849 | `.../fig5g/crossomics_347/{functional,native}_correlations.tsv`, `*_distances.npz` |
| A85 controlled grid, PREDICT | `BR/distance_sensitivity.py` with `experiment.py` (382d8f70...) | f23a5468d394757c4ba73eee691bde3d19f0896a1716dfe6b37f1d7205992d65 | 6044329 | `BR/distance_refit/controlled_176_comparisons.tsv` (PREDICT rows), `PREDICT_distances.npz` |
| A85 controlled grid, CapScan 347 | `BR/claude_review_2026-10-03/final/reruns_347/A4_controlled_distance/code/distance_capscan347.py`, `collect.py` (36f8cff0...) | d81c1293a9da1f045bf66b970716d88383a018e41b4dc2dd084fe363552a9fc3 | 6064864 (run), 6064865 (collect) | `.../A4_controlled_distance/run347/CAPSCAN_distances.npz`, `controlled_176_comparisons_347.tsv` |
| A86 native retention | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/native_pyramids.py` | 64ba668b266b9be42bb0a3f273878a2b93daba53022cc3c5a70ab04d4d56c8da | 5990407 | `.../Crossomics/native/SOURCE_DATA_retained_signal_by_profile.tsv` |

The A84, A82 and A85 launchers call the lab's miniconda Python 3.12.4 interpreter, which today has
numpy 2.2.6, scipy 1.16.3, pandas 2.3.3 and statsmodels 0.14.5 (the versions of
`environment/requirements.lock.txt`); package versions at the time of each run were not recorded.
The interpreter of the A86 run (job 5990407) is not recorded.

## Codex census reconciliation
| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| v28d embedded image | Supplementary Figure 11, `word/media/image19.png`, sha256 fa6b628a... | BYTE_MATCH to `iteration4/figures/out/SupplementaryFigure11_D12_400dpi.png`; write site `iteration4/figures/build/validate.py:111`; producer not rerun | v28d image produced by `d12_build_figures.py` + `d12_build_supplement.py` + validate.py step 6; reproduced byte for byte by `SF11/setup_and_build.sh old` (`SF11/out/SupplementaryFigure11_old_400dpi.png`) | yes | Codex lists the old reproduction among its byte-identical files; the census does not rerun producers, this unit relies on the recorded reproduction |
| v28f image | REBUILT_CANDIDATES.tsv rows for `SF11/out/SupplementaryFigure11_n15_400dpi.png` (9b0504c5...) and `SF11/out/SupplementaryFigure11_n15.pdf` (4294fcf6...) | REBUILT_CANDIDATE, write route `SF11/export.py:11` inspected, not rerun; differs from all v28d media | same files; producer route reproduced here: regenerated PNG has the same sha256 9b0504c5..., panels IDENTICAL | yes | a single rebuilt candidate for this figure, no selection needed |

## Source data
| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| `source_data/SFig11/functional_correlations_PREDICT.tsv` | `SF11/inputs_n15/functional_correlations.tsv` (PREDICT lines; = verified50 `cross_refit/functional_correlations.tsv` = package `Biology_refit/cross_refit/functional_correlations.tsv`, sha256 4c05b200...) | 882c36b28646498e5a7d1516b2af6c7f80397e95f3881f0fd1bcef169d11b6b7 | header and PREDICT lines copied as text (`prepare_source_data.py:split_rows`) |
| `source_data/SFig11/native_correlations_PREDICT.tsv` | `SF11/inputs_n15/native_correlations.tsv` (PREDICT lines; = package `Biology_refit/cross_refit/native_correlations.tsv` PREDICT rows) | 55d6a8c2aa0f79df9b7d7bbd0a0e845184a64ac2cf8ffc3e5f53992903abd5ea | as above |
| `source_data/Fig5/functional_correlations.tsv`, `native_correlations.tsv` (read, owned by Fig5) | CapScan lines of the same two producer inputs = `.../fig5g/crossomics_347/*_correlations.tsv` | 4c3468fc..., 0004f35d... | equality asserted by `prepare_source_data.py`; package sheets still hold 12-participant CapScan rows (do not agree) |
| `source_data/SFig11/controlled_176_comparisons.tsv` | `SF11/inputs_n15/controlled_176_comparisons.tsv` (PREDICT rows = `BR/distance_refit/controlled_176_comparisons.tsv` = package `Biology_refit/distance_refit/controlled_176_comparisons.tsv` PREDICT rows; CapScan rows = `A4_controlled_distance/controlled_176_comparisons_347.tsv`) | 1705248887b442a01f3c52b36c584b3107db38d5f8e5c53067c32081d506ffca | byte copy. `controlled_176_comparisons_347.tsv` (e0757c10...) has the same values but its PREDICT floats were re-written with one digit fewer; the producer read the `inputs_n15` copy. Package CapScan rows (12 participants) do not agree |
| `source_data/SFig11/distances_PREDICT_functional.tsv` | verified50 `cross_refit/functional_distances.npz` (= package `Biology_refit/cross_refit/functional_distances.npz`) | a8e0d5a9f7f2aa52f9dfbfad87a6ac720934a5458bb6f341bbabc49bd53648d0 | upper triangles, keys `PREDICT__displayed77__<level>__0.0__n77__<layer>`, repr floats, exact round trip asserted |
| `source_data/SFig11/distances_PREDICT_native.tsv` | verified50 `cross_refit/native_distances.npz` (= production `BR/cross_refit/native_distances.npz` = package sheet) | 8928b6d572273e1763cf3422fcd5fa28fd381fbbfce0e92ab5bd9d06efa65e51 | as above, keys `PREDICT__<analysis>__<layer>` |
| `source_data/SFig11/distances_PREDICT_controlled.tsv` | `BR/distance_refit/PREDICT_distances.npz` (= package sheet) | f596d528271810b8ad604b48c5a1f1b2f218870ba3b76d83d1c8661c2ed4c7ee | as above, 20 controlled keys; `wuyts_recipe` omitted |
| `source_data/SFig11/distances_CAPSCAN_native_prevalence10.tsv` | `.../fig5g/crossomics_347/native_distances.npz` | d7b0bb2573b8f984b1e026ca5f678f74c597384789d556ecb2d5bd69e5b4c0e3 | as above, keys `CAPSCAN__primary_prevalence10__<layer>`; the all-observed keys equal `source_data/Fig5/native_distances_long.tsv` (asserted) |
| `source_data/Fig5/functional_distances_long.tsv` (read, owned by Fig5) | `.../fig5g/crossomics_347/functional_distances.npz`, keys `CAPSCAN__matched12__<level>__0.0__n65__<layer>` | 7b9ba5caa5b78dd767ccb5c8087053d047379fbf5ea3760beabad5e694cb0508 | equality asserted by `prepare_source_data.py` |
| `source_data/SFig11/distances_CAPSCAN_controlled.tsv` | `.../A4_controlled_distance/run347/CAPSCAN_distances.npz` | 60c011d1de9e7851d2100ca1e7bf2001bd55091dba68f96af2ea749edf2322f6 | as above, 32 controlled keys; `wuyts_recipe` omitted |
| `source_data/SFig11/predict_native_retained_signal.tsv` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/SOURCE_DATA_retained_signal_by_profile.tsv` (not in the package) | 761a6b1242d432ae9f5b6fb1cf30221b5cd9b45a8cc492996aeaac57b8125524 | PREDICT `primary_prevalence10` DNA and RNA rows; patient label replaced by row position |

`SOURCE_DATA_MAP_reanchored.tsv` rows: RES-80 (line 139, q 0.0237391, VERIFIED) agrees with the
recomputed q 0.02373913043478261. SF11 (176 estimates) and SF11-perm (minimum P 0.0005) agree. The
SF11 row points to package sheets whose CapScan rows are still at 12 participants
(OPEN_ITEMS.md, ZENODO).
