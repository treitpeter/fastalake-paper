# Provenance: Supplementary Figure 10

In this file `IT2` abbreviates `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2` and `IT1` abbreviates `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B`. `LAB_ROOT/.../final/` abbreviates `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/`.

## Panels
| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| figure (v30, single panel) | LAB_ROOT/.../final/v30/figures/SFig10/build_sfig10.py (sha256 86bb5e3018da4553707483c34b97e77c35b7b0394c2b67095e27c56e1a0ab197), run on Slurm in the fastalake-paper_v28_work environment | LAB_ROOT/.../final/v30/figures/SFig10/SupplementaryFigure10.png; .pdf; plotted_values.tsv | png 9b85a672a7c39d87371ec10e5d5a8c5e7361590e150a33303d8c13bebb5a6fd8; pdf 7093d6897d946f253380df3a0f7431dc92bde1e85c3237783a1da427433d909b; plotted values de6c8bb6eef41fec138f0253497f0175c4bc785621c54033b73a7c50a9800d62 | not recorded in the v30 folder | CPython 3.12 (fastalake-paper_v28_work env), matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3; Söhne via fl_style |
| a and b (one canvas, v28f; superseded) | IT2/consolidate/suppfig10_merged/build_suppfig10_merged.py (sha256 eee87d7883394e6abfdb087eacab5ae1bc516fe8a478713c745bf318677a5ea0), launcher run.sbatch | IT2/consolidate/suppfig10_merged/out/SupplementaryFigure10_v20_verified50_notitles.png; .pdf | png 7cc12746de9e6907762fb8bd7dd5ea7f00ca25a8b5f8aae64c923475621b334a; pdf 9e7c863de2d6c0fe2b7fe56a2ede7d8290a0e94d3bd55048c62ee74361891144 | 6048828 (wrote the figure, then stopped at a comparison with a missing earlier render), 6048832 (complete run; wrote SUPPFIG10_VALIDATION.json, sha256 1ff488ae80f97c8f6536de930ce6a438b23aa49ef92f70a16a88e37455694632) | CPython 3.12.4 (miniconda), matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3, Pillow 12.0.0, PyMuPDF; typeface soehne-buch.ttf |
| a statistics, data for a and b | IT2/cr003_followup/scripts/build_approved_eds_cr003.py (sha256 0917829d274e9cc06b7f2715198ebabcc0f69a2e2cae0a5415b7219d49b871fc), run as IT2/cr003_followup/shadow/verified50/eds/build_approved_eds.py by IT2/cr003_followup/scripts/build_figs.sh | IT2/cr003_followup/shadow/verified50/eds/EDS_approved/source_data/ (BUILD_VALIDATION.json status NUMERIC_PASS_VISUAL_REVIEW_PENDING) | per file below | 6048652 | CPython 3.12.4, numpy 2.2.6, pandas 2.3.3, scipy 1.16.3 |
| screen behind a (candidate pairs, family q values) | IT1/eds_biochemistry.py (sha256 988ec3220b9ee42efac383bcbf4b5931d424e52f472878c5e0c973900c6fb686) running IT1/eds_refit/analyze.py (9925977db928a866b73ef424b34de45a99d04cc571c160b83f099e5d496b1a1a) and IT1/eds_refit/individual_followup.py (86ecade6b17e9291a38b57b72cf671af01a7dfa8fb627d5f2da0f52eb85fadfb) | IT1/eds_refit/augmented_pair_screen.tsv, all_EDS_associations.tsv | d7e59ab759d62aa5b2068ac5fbe0b774ece288ad3894bd633cfd219cac2a3967; dcd15919dde69cce09d3fdd3e96c5517e35f1786b824bb5575ca7dc0131008dd | 6048399 | CPython 3.12.4, scipy 1.16.3, statsmodels (BH cross-check) |

The repository notebook was run with the same interpreter and package versions
(`environment/requirements.lock.txt`). With FASTALAKE_FONT_DIR holding the Söhne fonts and
the controlled input present (the pseudonymised covariate table
8a4844354b3dfcdc8087e19def097e25dcfebfc59fabf1f8bfbad918f06baa68 and the v30 images, all in
FASTALAKE_CONTROLLED_DIR), `checks/SFig10_v30.json` (PNG) and `checks/SFig10_v30_pdf.json`
(PDF at 300 dpi) are IDENTICAL, 0 differing pixels. The plotted-values table rebuilt in memory
has the sha256 of the v30 build's plotted_values.tsv.

## Codex census reconciliation
| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| Supplementary Figure 10 v28f (a and b; superseded by v30) | word/media/image18.png, sha256 7cc12746... | BYTE_MATCH; WRITE_SITE_HASH_RECEIPT at IT2/consolidate/suppfig10_merged/build_suppfig10_merged.py:123; receipt SUPPFIG10_VALIDATION.json; producer not rerun | same write site and output; producer logic ported to panels.py and rerun here, output pixel-identical (PNG and PDF) | yes | REBUILT_CANDIDATES.tsv lists no alternative render for this figure |

## Source data
| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| patients.tsv | IT2/cr003_followup/shadow/verified50/eds/EDS_approved/source_data/patient_order_and_covariates.tsv; pseudonyms from Source Data package FastaLake_Source_Data/EDS_focused_20261003/source_data/patient_order_and_covariates.tsv; ko_total_signal from IT1/eds_refit/input_data/KO_profile.tsv.gz | 6d8b883529d90dbc682898c2025f01f62e781ba9a38c9ac072aac3a64ddbadc9; 8a4844354b3dfcdc8087e19def097e25dcfebfc59fabf1f8bfbad918f06baa68; 164fe6cde20cec1e6cc6803c7b2ad50a3f373f02e411e211769670b0f39cec02 | figures/SFig10/prepare_source_data.py: clinical group, antibiotic state and identifiers removed; package agrees (same rows, only patient_id pseudonymised) |
| focused_values.tsv | EDS_approved/source_data/focused_values.npz = package EDS_focused_20261003/source_data/focused_values.npz | e5900ff12c7f7e540f31a0bb1796e0be1bb6341216c436e9039145dd752ddc31 (identical in package) | npz to long table, sample ids replaced by position; exact round trip asserted |
| focused_features.tsv | EDS_approved/source_data/focused_features.tsv (identical in package) | 7945c40a934016cbfab1952c42b02cdf2854fb923b3815b166b308a81e352aae | copied |
| legacy_85_values.tsv | EDS_approved/source_data/legacy_85_values.npz (identical in package) | 11c171206d200b9051c0f785041288bbc621be61695baae1bfa8eaffb88b7b37 | npz to long table, sample ids replaced by position |
| legacy_features.tsv | EDS_approved/source_data/legacy_features.tsv (identical in package) | b883e1659743f5437c3d1d3981abd07a36d1ae41ee5d87b6321596fe90968e51 | copied |
| legacy_259_values.tsv | EDS_approved/source_data/legacy_259_heatmap_values.npz (identical in package) | b937bb50fc4c870173cd85638041805e8c5700ded180616be9cf91ef22c2f042 | npz to long table, sample ids replaced by observation index; the column table legacy_259_heatmap_columns.tsv (identifiers and clinical fields) is not copied |
| augmented_pair_screen.tsv | EDS_approved/source_data/augmented_pair_screen.tsv = IT1/eds_refit/augmented_pair_screen.tsv (identical in package) | d7e59ab759d62aa5b2068ac5fbe0b774ece288ad3894bd633cfd219cac2a3967 | copied |
| all_EDS_associations.tsv | EDS_approved/source_data/all_EDS_associations.tsv = IT1/eds_refit/all_EDS_associations.tsv (identical in package) | dcd15919dde69cce09d3fdd3e96c5517e35f1786b824bb5575ca7dc0131008dd | copied |
| selected_pair_summary.tsv | EDS_approved/source_data/selected_pair_summary.tsv (identical in package) | 42e89ff32dc0a4d3977feec112f3e1e827cb5940d82df4062d45df7f9cb8e8c1 | copied; recomputed by stat_a_unadjusted, stat_a_joint_adjusted, stat_family_q (agree) |
| separate_adjustments_all_models.tsv | EDS_approved/source_data/separate_adjustments_all_models.tsv (identical in package) | fe72b674538a021d48f14fc66c50d8aaa9e454960ca636969276b99a2133917b | copied; recomputed by stat_a_adjustment_grid, maximum absolute difference 0.0 over 192 estimates |
| producer_plotted_ab.tsv | Source Data package FastaLake_Source_Data/Supplementary_Figure_10/Supplementary_Figure10_plotted_values.tsv (= IT2/consolidate/suppfig10_merged/out/SupplementaryFigure10_plotted_values.tsv) | e08b1d79a049d1f0049c0d60ab2b5eea095723d0f518d14fb4232daac8c6755b | figures/SFig10/prepare_source_data.py:heatmap_rows: header and panel a and b rows copied byte for byte (2,125 rows); the 255 a_track rows (EDS, human share, per-patient antibiotic state) withheld. panels.py:render rebuilds the whole sheet in memory with the controlled input and its sha256 equals the origin; panels.py:stat_clip rebuilds the heatmap rows without it and they equal this file exactly |

Source Data map (SOURCE_DATA_MAP_reanchored.tsv) row SF10 points to the plotted-values
sheet, focused_values.npz and legacy_85_values.npz with status LOCATED; this unit
verifies all three (identical bytes or identical rebuilt sha256).
