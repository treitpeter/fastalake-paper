# Provenance: Supplementary Figure 4

## Panels

All seven panels are drawn into one composite. The v30 figure (first row) is the v29
figure below restyled without changing any plotted value: the v30 build compared every
axes' plotted arrays, labels, limits and ticks with the v29 producer run unchanged
(all identical) and the plotted-values table with v29 (identical); see
LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig04/CHANGES.md and values_check.json. Per-panel images
are crops of the composite at boxes computed from the drawn axes; there is no per-panel
producer render.

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a-g (v30 composite) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig04/build_sfig04_style.py (sha256 b950b406...) with _sf456_common.py (sha256 c2a03a2c...), run by run_sf456.sh | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig04/SupplementaryFigure4.pdf and .png | pdf 93e4d0a7a183e58d5256ede1e5204273d45221c1093f9f23742520fd01f36506; png fb3641ee8ed2f9f5fc152696d528686b25e76344832764498d7b3fbdf03de684 | not recorded (OPEN_ITEMS.md) | CPython 3.12.4, matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3, PyMuPDF 1.27.1 |
| a-g (v29 composite) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/suppfigs_n15/SF04/build/sf4.py + common.py (runner rebuild_panel_c.sh) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/suppfigs_n15/SF04/n347_panel_c/out/SupplementaryFigure4_R5_400dpi.png | 5f36158aefa697484f0786bf8721ab624f094d7f51304fc82a931e823e46fe5d | 6065948 | CPython 3.12.4, matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3, PyMuPDF 1.27.1 |
| a-g (v29 vector) | as above | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/suppfigs_n15/SF04/n347_panel_c/out/SupplementaryFigure4_R5.pdf | 5eb65e23c4bb0df038b58ae1bf670731bcadd85a89c83dd867644fef8a3332e5 | 6065948 | as above |
| plotted values | as above | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/suppfigs_n15/SF04/n347_panel_c/plotted_values/SupplementaryFigure4_plotted_values.tsv | e1a216dda1d258a84ae3e457e9b0d8fe67bb542c975039fae30d3ac0726a6c2c | 6065948 | as above |

Lineage of the producer: `sf4.py` and `common.py` are copies of
LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration5/v23/figfix/build/sf4.py (sha256 dd8e05895c8af0a79d3d2ac22827d071040a4f34cc0add2b005b613749a2af9f)
and common.py (sha256 75b58c7fa87711fab2a8532cac2e0098ecc586f557380a3100a85e405b2593db),
which produced the v28d image. The copy adds, behind SF4_MODE=new, only input changes:
the 16 participant-15 CapScan FL_public rows in c, the 347-roster CapScan grouping in g,
and the tick label "CapScan (347)". In old mode the copy reproduced the v28d image
byte for byte (LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/suppfigs_n15/SF04/logs/repro_old.log).

The font used by the producer and by this unit: Söhne Buch TrueType, `soehne-buch.ttf`,
sha256 db966efe63405323b2d1c2cfb3f104491d61ad56bf3a95fe3bfa712965e40ed8.

Upstream jobs: M09 searches, Slurm arrays 6001594-6001598 and 6002227 (receipts in
LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7/receipts); CAMPI member counts 6007364; 16 CapScan
FL_cohortMG searches and 347 regrouping 6064772-6064775
(LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/reruns_347/A2_minlen7/JOBS.txt); 16 CapScan FL_public searches and collector
6065307-6065309 (LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/public347/minlen7_public347/JOBS.txt).

This unit's checks (on Slurm, 2 cores, 7 s): the 300 dpi composite PNG is
byte-identical to the v30 PNG (`checks/SFig04_v30.json`, IDENTICAL); the composite PDF
rasterised at 300 dpi is IDENTICAL to the v30 PDF; the plotted-values table written by
`panels.build_figure` is byte-identical to the v29/v30 table.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| v28d composite | Supplementary Figure 4, word/media/image11.png, sha256 75d7c08e086b... | BYTE_MATCH to iteration5/v23/figfix/out/SupplementaryFigure4_R5_400dpi.png; write site iteration5/v23/figfix/build/sf4.py:138; receipt iteration5/v23/figfix/check/SF4_BUILD.json (sha256 41255e74...) | same producer lineage; the v28d image is the 331-roster render | yes | the v28d image is not the v28f target |
| v28f composite | REBUILT_CANDIDATES.tsv, suppfigs_n15/SF04/n347_panel_c/out/SupplementaryFigure4_R5_400dpi.png (5f36158a...) and .pdf (5eb65e23...) | REBUILT_CANDIDATE, write route sf4.py:157 inspected, not rerun | selected target (decision of the authors); reproduced here byte for byte | yes | |
| v30 composite | no census row (v30 postdates the census) | none | build_sfig04_style.py, reproduced here byte for byte | not applicable | restyle of the v28f/v29 composite above |
| superseded candidate | REBUILT_CANDIDATES.tsv, suppfigs_n15/SF04/n347/out/SupplementaryFigure4_R5_400dpi.png (06c4f9d7...) and .pdf (e88b4802...) | REBUILT_CANDIDATE | superseded: panel c still on 331 | yes | not used |

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| benchmark_matched_per_acquisition.tsv | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7/wave2/peptide_summary/matched_comparisons.tsv; package sheet a_b_benchmark_per_acquisition (agrees) | 3e2ae7b0af35cbe896a59160ab3bda36452ff517b040ce6c5bf6880a0ac23452 | text lines with arm FL_public and comparator MG (prepare_source_data.py:stage_source) |
| per_acquisition_length_counts.tsv | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7/wave2/peptide_summary/per_acquisition_wide_claude.tsv; package sheet c_per_arm_change (agrees at 331) | 4ba1f5dcd266997502f4e96de7dfc2897e524310df3b43baf37067cb650c058e | byte copy |
| per_acquisition_length_counts_participant15.tsv | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/public347/minlen7_public347/per_acquisition_new16.tsv; no package sheet (package is 331) | 53c749ed5a550d77cae7f97a2d4b183708dea2f31a0ab228ac5aac826e74223d | byte copy |
| campi_member_counts.tsv | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7/members/MEMBER_COUNTS_BY_LENGTH.tsv; package sheet d_member_counts (agrees) | 7033b5a9e91b29446e095a0227bfb89dee83b460fb6bef3480c7ffee63ce0f49 | byte copy |
| campi_fusion_pooled_q01.tsv | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7/fusion375/POOLED_q01.tsv; package sheet e_fusion (agrees) | 669a785968dca936fe9ed9acdc45bf205558ca711c845bf77459374e4e3bc073 | byte copy |
| campi_fusion_per_acquisition_q01.tsv | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7/fusion375/receipts/summary_00.tsv .. summary_29.tsv; no package sheet | multiple files | text lines with q_threshold 0.01, files in name order |
| healthy_strategies_per_acquisition.tsv | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7/wave2/healthy_strategies/per_acquisition.tsv; package sheet f_healthy (agrees) | 8f8e7b225b9c9652c382e2a3bf875079a230d0f685f06c281be4f3f069aa598a | byte copy |
| grouping_study_groups.tsv | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7/grouping/COMPARE_{CAMPI_FL_REF,CAMPI_FL_REF_MG,CAMPI_FL_UNIPROT,CAMPI_MG,CAMPI_REF,HSTOOL_noclust,PREDICT_noclust}.json and LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/reruns_347/A2_minlen7/run/grouping/COMPARE_CAPSCAN_noclust.json (sha256 41eb205719dc7f4dcc696161b97142fa388742c7663d530c48d4a343563899e2); package sheet g_grouping (seven rows agree; CapScan differs by design: package 174,673 -> 175,282 at 331, here 177,717 -> 178,357 at 347) | per file in BUILD record SF4_BUILD.json | one row per JSON: task, arms_used, comparison.acquisitions, comparison.inferred_study_groups |
| search_inventory.tsv | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/M09_short_peptides/run_2026-09-24_minlen7/tasks.tsv, LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/reruns_347/A2_minlen7/run/tasks.tsv (new indices 331-346 only), LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/public347/minlen7_public347/tasks.tsv, each joined to receipts/receipt_<index>.json; no package sheet (Source Data gap G24) | multiple files | config and path columns dropped |

Package workbook: LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Figure_04/LENGTH_SENSITIVITY_source_data.xlsx,
sha256 b74bac2564ea799f2b536c5bce6495fae933e91606f25da0feacabb96c539345 (identical to the M09
builder's workbook). Sheet-by-sheet comparison: `checks/source_data_vs_package.tsv`.
