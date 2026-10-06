# Provenance: Figure 3

## Panels
| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a | none (artwork extracted from the v6.1 manuscript figure) | MM_review_analysis/figure_restoration_2026-09-16/Figure3/panels_61/a_design.pdf | 1503196bb523aee098e7f899533431760b9d10aab2fb8ea8f9dd7bd44470cb07 | none | cairo 1.18.4 (PDF producer field) |
| b-e | MM_review_analysis/figure_restoration_2026-09-16/Figure3/build_figure3_cdef.py (sha256 2ba6574e5ba4fb588cbc5ada0c261a735771a1719955beea29dffa424a87b9c0), style FastaLake/figure_style/fl_style.py of 2026-08-02 (sha256 58b64c343f095cbace09c3e9766e3e78060c1c8cf522b075c4047300322164c3) | MM_review_analysis/figure_restoration_2026-09-16/Figure3/panels_cdef/Figure3_panels_cdef.pdf (and .png) | cf1fba8b3512a4ffc0af9c7725e15848dcd3814ea5c98e1845350287c56cb517 (png f71b5dec70a88d881e6614ff7f8b1a1398f62cd74e4f45cae8989eb234fa477a) | 5948528 (first render 5948289) | anaconda 2023.03 CPython 3.10 with user-site matplotlib 3.10.8 (PDF creator field), pandas 2.3.3, numpy 1.26.4 |
| f | review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/cr003_followup/scripts/build_fig3f_cr003.py, argument verified50 (sha256 65b68dac85e1b69860b491136b3e949c9891ae5178f4d76023099d7d948f8409) | .../iteration2/cr003_followup/shadow/verified50/fig3/Figure3f_proposed.pdf (and .png) | de6466793d5091295aee8c71a6872339ddfce9f4df2c983f1c31ac56a7e92d35 (png cdd3b56bd2d0f639b5f377f693b215d1a34ec6de4080397df013e4086616cce1) | 6048652 | miniconda3 CPython 3.12.4, matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3, PyMuPDF 1.27.1 |

Regenerated here with CPython 3.12.4, matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3, PyMuPDF
1.27.1 (`environment/requirements.lock.txt`). Pixel checks at 300 dpi (`checks/panel_bcde.json`,
`checks/panel_f.json`): both IDENTICAL, 0 differing pixels. The panel f matrices and the 20 label
summaries also equal the producer's receipt (`expected/panel_f_producer_receipt.json`, extracted
from `.../shadow/verified50/fig3/FIGURE_VALIDATION.json`, sha256 4281dfe6...).

Composite: `expected/Fig3_v28f.pdf` and `.png` are `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/fig_v27/Figure3.pdf` and `Figure3.png`, written by `final/fig_v27/build_fig1_3_v27.py`. The full chain is in `ASSEMBLY.md`.

## v30 (9-50 residues, relabelled legend)
| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| b-e (v30) | review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/len950/Fig3/build_fig3_len950.py (sha256 827f4359f90540866c566170eb97a518bc0d07d7083c5463b03aabfa52167313) with panels_len950.py (sha256 54e4ee83acffe4c5e662b329285a5294b8e6f541795e1da106c805158ff10230) | .../final/len950/Fig3/Figure3_len950_panels_bcde.pdf (and .png) | 27247c89159e291d249d94db4663abc77b39bd944dddd15bda3c1c31c01e14e1 | not recorded | not recorded (same interpreter as the repository env per the producer's import of repo panels.py) |
| composite (v30) | .../final/fig3_labels/relabel.py (sha256 a7993917f97dea43255f31fe549410dc2a8b0a631ffcbcd0af5b8b234c547519), then build_fig3_len950.py, then .../final/v29/figures/_tools/flatten.py (sha256 d73f26fc54d59e42e3fa439b11582f1c4158ff9e309f5ec69bab26dda2ab8cd1) | .../final/len950/Fig3/Figure3_len950.pdf; .../final/v29/figures/Fig3/Figure3.png and Figure3.pdf | pdf a4f9c36a41187c8fd86e4e593ca2143acda83103fa7703df730fbbc5ab36bf9c; png 13589443fb69fa47dfde7a592c3ea3221548f63f4cc2ff07d4287cfc3753ec9d; flattened pdf d6b5c2174750d64ff96067e76bc58c3e2bb0cc593d69d0ccbab0b44c891fa60e | not recorded | PyMuPDF (version not recorded) |

v30 values: `LAB_ROOT/.../final/len950/Fig3/NEW_VALUES.tsv` (written by `new_values.py`, sha256 2faa8dc3...) and
`LAB_ROOT/v28f_9-50_FINAL_REVIEW.md`; v30 wording from `LAB_ROOT/.../final/v30/visible_v30f.txt`. Every
NEW_VALUES.tsv b-e median equals the notebook's recomputed value at full precision.

Regenerated here (Slurm job 6075222, `tools/run_notebooks.py figures/Fig3`): `checks/panel_bcde.json` IDENTICAL
against the v30 producer render, `checks/panel_f.json` IDENTICAL, `checks/composite_v30.json` IDENTICAL against
the v30 figure PNG (0 differing pixels each).

## Codex census reconciliation
| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| whole Figure 3 | Figure 3, word/media/image3.png, sha256 6d72728a... | BYTE_MATCH to final/fig_v27/Figure3.png; producer WRITE_SITE final/fig_v27/build_fig1_3_v27.py:72; "producer was not rerun" | same final write site; `expected/Fig3_v28f.png` is byte-identical | yes | Codex lists no Figure 3 row in REBUILT_CANDIDATES.tsv, and this unit found no rebuilt Figure 3 newer than fig_v27 |
| b-e | (inside image3) | not attributed below the composite | b-e producer build_figure3_cdef.py, rerun here: IDENTICAL | adds | text in the composite was later re-typeset (ASSEMBLY steps 8, 9) |
| f | (inside image3) | not attributed below the composite | f producer build_fig3f_cr003.py verified50, rerun here: IDENTICAL | adds | texts removed and keys added later (ASSEMBLY steps 12, 13) |
| a | (inside image3) | not attributed below the composite | artwork panels_61/a_design.pdf, no drawing script | adds | ASSEMBLY step 3 is untraced (OPEN_ITEMS) |

## Source data
| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| source_data/Fig3/transfer_pool_comparisons.tsv | MM_review_analysis/A15_4_4_transfer_pool/results/all_comparisons.tsv. Not in the Source Data package (only medians are, `Figure3/comparison_summary.tsv`); recomputing that sheet from this table with A15 `report.py`'s grouping gives all 36 rows equal (largest absolute difference 1.4e-14) | a7618438c577a20b6037bf801b4daa89c181dcafb724066cfada229b06035947 | byte-identical copy |
| source_data/Fig3/protein_group_pairwise.tsv | MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/HSTOOL__simplified/source_data/Figure3/pyramid_pairwise_local_candidate.tsv. Not in the package (only medians are, `Figure3/reproducibility_display_summary.tsv`); the 10 "Protein groups" rows of that sheet equal the medians and pair counts recomputed here (largest absolute difference 3.6e-15) | b45b6398514eebef4d92ac080f5df139f75fb2efcde79aebccfc475cab787317 | byte-identical copy |
| source_data/Fig3/ko_pairwise.tsv | .../claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/hstool_precision_refit/pairwise_native_features.tsv; byte-identical to package `HSTOOL_KO_refit/revised_comparison/pairwise_native_features.tsv`; its 20 profile-design medians equal package `native_feature_summary.tsv` exactly | 2cf176c091f92f89162acfb091c3e3da6ca04886380f564c158040ac521a01f3 | byte-identical copy |
| source_data/Fig3/design.tsv | MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/HSTOOL__simplified/source_data/Figure3/design.tsv (= canonical/publication_data/arms/noclust_2026-08-07/_figroot/publication/figures/Fig4/tables/design.tsv); byte-identical to package `HSTOOL_KO_refit/revised_comparison/design.tsv` | 01ce64a542cfd664eba560aa99c03fc04546045951d74986b44527e55bc6c2c3 | byte-identical copy |

Package sheet that does not agree with the figure: the 10 "Functions (KO)" rows of package
`Figure3/reproducibility_display_summary.tsv` (labelled host and microbial; for example r 0.905
technical and 0.765 inter-individual) equal the medians of the superseded KO rows of
`protein_group_pairwise.tsv`, not the plotted and quoted KO values (OPEN_ITEMS).

Source Data map (`SOURCE_DATA_MAP_reanchored.tsv`): RES-48, RES-51, RES-52 and RES-54 point to
the sheets above and agree with the values recomputed here. FIG3-f is LOCATED only, and the KO
values on line 109 have no correct map row (OPEN_ITEMS).
