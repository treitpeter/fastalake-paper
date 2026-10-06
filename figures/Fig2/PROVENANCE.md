# Provenance: Figure 2

## Panels (v30)

All eight panels and the composite come from one producer run of the v30 (v29 relabel of the len950)
build, which reads this repository's `source_data/` (`Fig2/`, `SFig03/CAMPI_membership_*.tsv`,
`Fig1/acquisition_sizes_726.tsv`) and `artwork/bottle_61.pdf`.

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a-h and composite | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v29/figures/Fig2/build_fig2_len950.py (sha256 e45b08327ec98300908234cfd082149e745f451e662112e8d893648c8c03076e) with panels_len950.py (sha256 5dd8ccf9eadb7857d8aef504ec8a1e030bc8a795c7bd432376e21a7f02c64fa6) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v29/figures/Fig2/Figure2.pdf / Figure2.png (BUILD_LOG.json) | PDF d74e1be857204c462eabba8fd96442007890c4a18c5fad66631a8259b495a336; PNG ec359e7ff04f0ae2e8434e2da64da108736bb2dc8631539b2ce0d66bbaf2965f | not recorded (CHANGES.md: "Run on Slurm"; BUILD_LOG.json has no job id) | this repository's environment (fastalake-paper_v28_work/env) |

The len950 predecessor (LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/len950/Fig2/, same code with "Identified"/"identified" axis titles) is
documented in LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/len950/Fig2/TEXT_EDITS.md and NEW_VALUES.tsv; its values equal the v30 values here.

Regeneration here: Slurm job 6075249 (p.hpcl92/p.hpcl8/p.hpcl93, 2 cores, 16 GB; notebook 33 s),
repository environment with the Söhne fonts. Composite IDENTICAL to the v30 PNG (0 differing pixels;
same sha256). The v30 producer kept no per-panel renders, so there are no per-panel checks; the
v28f per-panel receipts (5-50 canvas, all IDENTICAL) were superseded.

---

# Provenance: Figure 2 (v28f, history)

## Panels

All eight panels come from one producer run that drew one canvas. The producer ran on the login
node without Slurm on 2026-09-22 (its README: "No new search or Slurm job was needed"); the canvas PNG
was refreshed at 300 dpi from the same PDF at 2026-09-22T10:19 UTC (FINAL_REVIEW.json, `pixel_identical_to_pdf_render: true`).
Producer letters differ from v28d letters; the mapping is in `panels.py` and `ASSEMBLY.md`.

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a (producer a) | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/figure2_redesign/build_campi_figure.py (sha256 dc20ed0b1b9f7edcf5d95d902907d9f3674228a24c0298864441de4939200400) | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/figure2_redesign/CAMPI_Figure2_review.pdf / .png | PDF b8c042347e03c06f1621bb9d035cb51e869daa04df602490d1b599084ccf3632; PNG 1a8fdc60742fd424cfc70a74b6988bcfaa2bb4d0988cc8227453ed8ac6131d60 | none (login node) | anaconda 2023.03 CPython 3.10 with user site-packages: matplotlib 3.10.8, numpy 1.26.4, pandas 2.3.3, PyMuPDF 1.27.1; house style LAB_ROOT/figure_style/fl_style.py (sha256 58b64c343f095cbace09c3e9766e3e78060c1c8cf522b075c4047300322164c3) |
| b (producer b) | same | same | same | same | same |
| c (producer c) | same | same | same | same | same |
| d (producer d-f, left) | same | same | same | same | same |
| e (producer d-f, right) | same | same | same | same | same |
| f (producer g) | same | same | same | same | same |
| g (producer h) | same | same | same | same | same |
| h (producer i) | same | same | same | same | same |

Regeneration here: repository environment (CPython 3.12.4, matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3,
PyMuPDF 1.27.1) with the Söhne fonts. Canvas and all eight panel crops are IDENTICAL (0 differing pixels,
`checks/*.json`). The only change needed was to restore the producer-era line weight (0.6 pt): the
repository's style module is a later version with a 1 pt floor, which made 0.35% of pixels (axis spines
and ticks) differ in a first trial; `panels.py:_producer_rc` and `producer_style` restore the producer's values.

`expected/panel_<x>.png` are crops of the producer's own 300 dpi PNG with the boxes in
`panels.PANEL_BOXES_MM` (`panels.crop_panel`), made once on 2026-10-06; the regenerated panels are cropped
with the same boxes.

The community artwork `artwork/bottle_61.pdf` (sha256 cb70ae529313330d33cb36a6a74730e58c8ba87af2a477a255f1e6595b4f40a1,
cairo 1.18.4, 2026-09-16) comes from LAB_ROOT/MM_review_analysis/figure_restoration_2026-09-16/Figure2/panel_a/bottle_61.pdf.
Its build script says the bottle was "reused from the 6.1 vector" (an earlier figure version); the only
record of its extraction is Slurm log `bottle_5948234.log` (page size and font only). No generating
script was found.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| whole figure | Figure 2 / word/media/image2.png (sha256 deb175c906964d33...) | BYTE_MATCH to final/fig_v27/Figure2.png; producer status COPY_ONLY (fig_v27/build_fig1_3_v27.py:75 only copies); "original plotting/rendering producer unresolved" | byte identity agrees (`expected/Fig2_v28f.png` has the same sha256). Plotting producer: build_campi_figure.py, reproduced pixel-identically here (step 0). Rendering: the PNG is a 300 dpi render of fig_clean/Figure2.pdf (verified IDENTICAL). Edits between: six PDF steps (`ASSEMBLY.md`), step 6 without a saved script | yes for the file; this unit resolves the producer | the Codex census did not trace back through fig_clean, FAR16, FAR15, fix_current_figures, FIG2F and r12 to the plotting script |

## Source data

All copies are byte-identical to their origin unless the transformation column says otherwise.

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| Fig2/paired_reference_recovery.tsv | LAB_ROOT/MM_review_analysis/A14_4_3_campi_gap/results/paired_reference_recovery.tsv; package FastaLake_Source_Data/Figure2/paired_reference_recovery.tsv (identical) | c446dabab7ad3c39c0d055424f60e19c50411eaebe82ee06fae52c54a4f74919 | copy |
| Fig2/recovery_summary.tsv | LAB_ROOT/MM_review_analysis/A14_4_3_campi_gap/results/recovery_summary.tsv; package Figure2/recovery_summary.tsv (identical) | 74658e932aea9f40a62245b6f3153f578289343182a02581ecca83736c060b13 | copy |
| Fig2/member_summary.tsv | LAB_ROOT/MM_review_analysis/A14_4_3_campi_gap/results/member_summary.tsv; package Figure2/member_summary.tsv (identical) | bdf6d5f0ded8bc8455f95bd1ac8c4ab52af38294ac7030959a269715873591b8 | copy |
| Fig2/stage_summary.tsv | LAB_ROOT/MM_review_analysis/A14_4_3_campi_gap/results/stage_summary.tsv; package Figure2/stage_summary.tsv (identical) | bfd51133b19b7c4bb6a5594ae4efe84aac9f8b6a1c632ff9c25f9874b587bd05 | copy |
| Fig2/length_by_member_class_a14.tsv | LAB_ROOT/MM_review_analysis/figure_restoration_2026-09-16/Figure2/source_data/length_by_member_class_a14.tsv (derive_upset_a14.py, Slurm 5948521); package Figure2/ (identical) | 76216e6553833ed793497c320d2298ce641bb8b81820a997c54d33dc2a926c9c | copy |
| Fig2/upset_intersections_FL_9_50.tsv | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/len950/Fig2/upset_intersections_FL_9_50.tsv (new_values_fig2.py, `panels_len950.upset_at` recount of the 15 membership tables at 9-50) | 1f7e4d0e4b01a59eb11a9059641a26274134e48adf83f85db032092c3004035b | copy; the notebook recounts it and asserts equality |
| Fig2/upset_intersections_a14.tsv | LAB_ROOT/MM_review_analysis/figure_restoration_2026-09-16/Figure2/source_data/upset_intersections_a14.tsv; package Figure2/ (identical) | 7f0405f42961e168297824f1118e573453a6edc3fe3cb464792c063ac79c20b5 | copy |
| Fig2/metrics_per_acquisition.tsv | LAB_ROOT/MM_review_analysis/figure_restoration_2026-09-16/Figure2/supp/source_data/metrics_per_acquisition.tsv (= LAB_ROOT/v4_simplified_grouping_2026-09-11/source_data/Figure2/retained/metrics_per_acquisition.tsv); package Figure2/ (identical) | db3499ea77963196082fab113695260bd6d3c4642c71c5811b4ae639186462e9 | copy. Its manifest names LAB_ROOT/canonical/publication_data/analyses/fig2_uniprot_2026-09-09/metrics_per_acquisition.tsv (sha256 1826576e...) as origin; checked here: same 75 rows and identical values in all 23 shared columns (including db_records); the retained table adds three grouping columns (mean_assigned_peptides_per_group, quantified_groups, assigned_coverage_pct) |
| Fig2/arm_CAMPI_FL_UNIPROT__reaggregated_per_acquisition.tsv | LAB_ROOT/simplified_grouping_results/manuscript_cohorts_2026-09-10_2200/cohorts/CAMPI_FL_UNIPROT/reaggregated/per_acquisition.tsv; package Figure2/arms/ (identical) | 22b986eb35f3086641f23c70d6151effe43b29386beadd9adede3853b9205fac | copy, renamed |
| Fig2/arm_CAMPI_FL_UNIPROT__direct_acquisition_summary.tsv | .../cohorts/CAMPI_FL_UNIPROT/direct/acquisition_summary.tsv; package (identical) | 9180536d0a2e4989f54c6fe8f7c024b9b60cfbb90a0c0d678f921598b0f60d12 | copy, renamed |
| Fig2/arm_CAMPI_REF__reaggregated_per_acquisition.tsv | .../cohorts/CAMPI_REF/reaggregated/per_acquisition.tsv | f0dc2fb2ae86939beba26d53fa68ce74d746203eb4198c077195ecefbe89643c | copy, renamed |
| Fig2/arm_CAMPI_REF__direct_acquisition_summary.tsv | .../cohorts/CAMPI_REF/direct/acquisition_summary.tsv | 9b797660a6b296cf1d464e5413bdf9e51ec19daa507f5bc36f9d7c55f1ad8b32 | copy, renamed |
| Fig2/arm_CAMPI_MG__reaggregated_per_acquisition.tsv | .../cohorts/CAMPI_MG/reaggregated/per_acquisition.tsv | 0bf2d7e592b02eb3baeab8b0b3eb2a38f97ac27572aa84dd12d844dcf2ab2194 | copy, renamed |
| Fig2/arm_CAMPI_MG__direct_acquisition_summary.tsv | .../cohorts/CAMPI_MG/direct/acquisition_summary.tsv | 02b6b172fb2123384e50518df107c815c1b310a9b0f5f2a88b0feeb7bb20aba9 | copy, renamed |
| Fig2/arm_CAMPI_FL_REF__reaggregated_per_acquisition.tsv | .../cohorts/CAMPI_FL_REF/reaggregated/per_acquisition.tsv | 130d4b89fa15644a2fcbf75dd30cbc4b39dc8cf76e8abb0bf66a9481ab25be39 | copy, renamed |
| Fig2/arm_CAMPI_FL_REF__direct_acquisition_summary.tsv | .../cohorts/CAMPI_FL_REF/direct/acquisition_summary.tsv | 9ddf09061b1f092444d8f6eebf3ae389b85b8516e1c8041beab44f154d4081f4 | copy, renamed |
| Fig2/arm_CAMPI_FL_REF_MG__reaggregated_per_acquisition.tsv | .../cohorts/CAMPI_FL_REF_MG/reaggregated/per_acquisition.tsv | 7e83bf1ab6548b4e2f5ad156ce94e3ca910270ce89739864a4371fd7081e5adb | copy, renamed |
| Fig2/arm_CAMPI_FL_REF_MG__direct_acquisition_summary.tsv | .../cohorts/CAMPI_FL_REF_MG/direct/acquisition_summary.tsv | e685e28241ad22d097525e3d4d07d50ee410479187f10ee66a287e387807dabd | copy, renamed |
| SFig03/CAMPI_membership_<sample>.tsv (15 files, stored by unit SFig03) | LAB_ROOT/MM_review_analysis/A14_4_3_campi_gap/results/<sample>/peptide_membership.tsv; package Figure2/peptide_membership/<sample>.tsv (identical) | per file in `panels.MEMBERSHIP_SHA256` (equal to the producer's VALIDATION.json sources) | read in place, sha256 checked |
| Fig1/acquisition_sizes_726.tsv (stored by unit Fig1) | see Fig1 PROVENANCE; CAMPI rows equal package Review_v16/Database_sizes/database_sizes_710.tsv (reservoir 375,333,654) | fd9ade3bb4ebd2ee6c129f69819b9697cf20abd7c10c77a12644bef25396e386 (repository copy) | read in place, sha256 checked |
| Fig2/reservoir_protein_presence_by_proteome.tsv | LAB_ROOT/MM_review_analysis/R8_open_requests_2026-09-17/B_blind_reservoir_overlap/results/full/protein_presence_by_proteome.tsv (overlap.py, Slurm 5954792); not in the package; the Codex containment task proposes it as Figure2/reservoir_containment/protein_presence_by_proteome.tsv | 032c84a3383dbd16ecccb95d586c3fc3b2f131d95c4f1757d4ab0983005c68a8 | copy, renamed |
| Fig2/campi_gain_availability.tsv | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/campi_variability_2026-09-30/mechanism_v2/gain_availability.tsv (Slurm 6031428); package CAMPI_variability/mechanism_v2/ (identical) | c9bc2d53c7af91c1346f837e529f95af300a50b3b5bd67bf13a5fb62c76fd74e | copy, renamed |
| Fig2/campi_new_quantified_peptides.tsv | .../campi_variability_2026-09-30/mechanism_v2/new_quantified_peptides.tsv (Slurm 6031428); package (identical) | 609937efda054880dce03eecb1cd2a61522a7850b12fa18c0bb0cbff7e063801 | copy, renamed |
| Fig2/campi_depth_standardization.tsv | .../campi_variability_2026-09-30/summary/depth_standardization.tsv (analyze.py collect, Slurm 6031314); package CAMPI_variability/summary/ (identical) | 7e2073608b4317cb6b0b642d755deebbd5217716029b6a697a8e507a96db03f5 | copy, renamed |
| Fig2/campi_FL_5_frequencies_<sample>.tsv (15) | .../campi_variability_2026-09-30/results/<sample>/FL_5_frequencies.tsv (Slurm 6031312, 6031313); package CAMPI_variability/results/ (identical) | identical to origin, checked per file 2026-10-06 | copy, renamed |
| Fig2/campi_discrepant_prediction_pooled_summary.tsv | LAB_ROOT/MM_review_analysis/A03_4_2_near_misses/results/pooled_summary.tsv (report_v2.py; Slurm 5920659-5920662, 5920711); package Results_tables/CAMPI_discrepant_prediction_pooled_summary.tsv (identical) | 42e3fd63a3562062981d16f45d749833913514f3191d0b09f9af867e1472e2d7 | copy, renamed |
| Fig2/alphanovo_accuracy_pooled.tsv | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/alphanovo_accuracy/results/pooled.tsv (accuracy.py, Slurm 6064877); not in the package | 6ad34ac899e76f30e1f01f5f9105ada2b85bd6996afd6dd4f14f50bcb48397b3 | copy, renamed |
| Fig2/alphanovo_accuracy_per_acquisition.tsv | .../final/alphanovo_accuracy/results/per_acquisition.tsv (Slurm 6064877); not in the package | 82bc0030c9cafaa8b68ea79d00037dc2f44e62374f1801cbb361eeb95ad9a18c | copy, renamed |

Comparison with `SOURCE_DATA_MAP_reanchored.tsv`: map rows RES-25a/b, RES-26a-c, RES-27a-e, RES-31,
RES-32, RES-32b, RES-33, RES-36, RES-37 and RES-38 point to the tables above and their values agree
with the recomputed values in `NUMBERS.tsv`. Rows RES-CM1-CM15 point to `member_summary.tsv`, whose
80 pooled counts the notebook recounts from the membership tables (`recount_membership`); the counts
used in the text (12, 7, 36, 86, 201, 290, 433, 229, 322, 500) agree. Rows RES-28, RES-29, RES-30, RES-34, RES-35
(LOCATED only) are now recomputed and agree (98.4%, 99.9%, 3.38/1.98, 84.8%, 20.3%). The map has no row
for 91.3-100% (line 75; RES-24 not in the archive), the AlphaNovo accuracy numbers (line 81), "six members"
and "3- to 17-fold" (line 83); this unit supplies their tables and code. Several map rows carry
"anchor not in text" flags; the v28d anchors used here are the line numbers in `NUMBERS.tsv`.
