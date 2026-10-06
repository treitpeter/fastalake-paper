# Provenance: Abstract and Discussion statistics

This unit has no figure panel. The table below lists the original producers of the
statistics instead.

## Panels

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| text: cross-fitted gain medians and intervals (lines 35, 183) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration4/crossfit/run/analyze_final.py (sha256 835e0c93c04c41d81a60dc6fbbe1a61ffe25055b3a918af175d9cf1abdd0b2da) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration4/crossfit/run/SUMMARY.tsv; LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration4/crossfit/run/SUMMARY_CAPSCAN_BY_REGION.tsv | abd33b20058b9cec3374d5b63f2afc9fd8718eeccf79e92421e17e54b070d85a; 64b4249c4bdde6cb47aa057bae06dc2277d9529cd8bacc1678d28c0c1a10a5cc | final step 6054853 (slurm-6054853.out beside the outputs) | Python 3.12.4, NumPy 2.2.6, pandas 2.3.3 (interpreter named in final.sbatch) |
| text: same, 347-acquisition CapScan pool (v28f roster) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/public347/crossfit_capscan347/run/analyze_final.py (paths changed only; sha256 6fd5da9be6dae98d49e753bcfe5a24cf6413b6040d8109076411d7ba631812b1) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/public347/crossfit_capscan347/run/SUMMARY.tsv; .../SUMMARY_CAPSCAN_BY_REGION.tsv | 4449fbe156c1abc2177db2e9c3d7d3bac378a5802a948b9037b150508e87fe0f; 4c14a59eb2a002232e767ca88ea76303592428025d18519509ec2495fe4cba5c | 6065360 | as above |
| text: 726 databases, human reference (lines 35, 181) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/fig1bd_n15/code/sizes.py, assemble.py; LAB_ROOT/MM_review_analysis/A01_4_1_database_sizes/frozen/audit_v2.py | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/fig1bd_n15/acquisition_sizes_726.tsv; LAB_ROOT/MM_review_analysis/A01_4_1_database_sizes/results/human_reference.json | fd9ade3bb4ebd2ee6c129f69819b9697cf20abd7c10c77a12644bef25396e386; b705dfc4aba8718c83939a812e3de6614ce0a72eb9260c3e0c7b6ff632cf1079 | see figures/Fig1/PROVENANCE.md | see figures/Fig1/PROVENANCE.md |
| text: reservoirs and computing time (line 193) | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/supplements/revision_2026-09-24/build_revision.py, Computing block lines 201-228, rerun for 347 acquisitions | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/reruns_347/A7_other/resource_accounting/fixed347/source_data/compute_stage_summary.tsv; .../compute_task_records.tsv | 37c1d94d821d5ba9defb5ff5c8495905ce93209a136952dec8b600b74816498f; ea24e72fe85b7e4c0755a849d5822a5347d9d3685e73f765a641072c2d7551b8 | see figures/SFig15 | see figures/SFig15 |
| text: eight CAMPI members (line 35) | LAB_ROOT/MM_review_analysis/A14_4_3_campi_gap/frozen/aggregate.py | LAB_ROOT/MM_review_analysis/A14_4_3_campi_gap/results/member_summary.tsv | bdf6d5f0ded8bc8455f95bd1ac8c4ab52af38294ac7030959a269715873591b8 | see figures/Fig2/PROVENANCE.md | see figures/Fig2/PROVENANCE.md |

Reproduction of the original values (Slurm, repository environment Python 3.12.4, NumPy
2.2.6, pandas 2.3.3): `stats.stat_crossfit_gain` reproduces every 9-50-residue
`gain_cf_pooled` median and both interval bounds in both original SUMMARY files for all
nine strata, to within 7.2e-15 percentage points (floating-point round trip through the
stored tables). The stool interval of the region stream (-8.36 to -2.53%) and of the main
stream (-8.36 to -2.26%) are both reproduced; Methods line 777 states the difference.

v28f roster check. The iteration-4 run (CapScan evidence pooled over 331 acquisitions)
and the 347-acquisition rerun differ in 3 of 578 acquisitions by one peptide each
(crossfit index 325 MG 293 vs 292, 368 FL 3,751 vs 3,750, 387 MG 3,001 vs 3,000; all
CapScan intestinal fluid), which Methods line 777 attributes to the unstable discriminant
solve. Every stratum median and interval at 9-50 residues is identical between the two
runs. `source_data/Fig4/crossfit_per_acquisition.tsv` holds the iteration-4 values, which
the v28f Figure 4b draws (`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/fig4bc_crossfit`), so the Abstract and Figure 4b
are computed from the same table.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| none | none | none | none | NA | this unit has no image; the census covers figure images only |

Codex number recheck (CODEX_TASKS_2026-10-05/out_recheck/REPORT.md), rows on this unit's
numbers:

| Codex row | value | Codex result | this unit | agree? |
|---|---|---|---|---|
| N0072 | PREDICT -3.5 | MATCH, -3.5384714440068077 | -3.5384714440068077 | yes |
| N0073 | CapScan stool -5.5 | MATCH, -5.50576989948352 | -5.50576989948352 | yes |
| N0075 | 21-36 | MATCH, 21.127414023081826 to 35.82251701069244 | the same | yes |
| N0205 | PREDICT interval -4.8 to -1.4 | MATCH | -4.795255221149597 to -1.3745913650359431 | yes |
| N0206 | CapScan stool interval -8.4 to -2.3 | MISMATCH, -8.359531722506663 to -2.4656860269229974 (independent draw order) | -8.359531722506663 to -2.2614498630473157, exact with the original stratum order of the main stream | the bound depends on stream and stratum order; the original order reproduces the printed -2.3. Interval not printed in this unit (Results line 123) |

Codex reports 97 PREDICT bootstrap clusters as 93 patients plus four acquisitions without
a patient identifier; the Fig4 table has the same four (`PREDICT_NOID_1` to `_4`).

Codex statistics register: rows S043, S044 (the two stool medians) and S046-S049 (the four
regions) match this unit's STATS.tsv for estimand, n, draws and seed; the software and
generator they list as missing are filled in STATS.tsv from the code.

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| source_data/Fig4/crossfit_per_acquisition.tsv (owned by Fig4) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration4/crossfit/run/PER_ACQUISITION.tsv | f7616394a3ffd829c92c2256244f5ce174116f7596b579f5db06b7b56a95a6e3 | figures/Fig4/prepare_source_data.py (pseudonyms, column subset). No package sheet: the Source Data map rows ABS-02, ABS-02b, ABS-03, ABS-03b still point to the superseded same-spectrum values in Review_v16/Confidence/ordinary_all_pairs.tsv (9-14%, 47-64%; classified GONE) |
| source_data/Fig1/acquisition_sizes_726.tsv (owned by Fig1) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/fig1bd_n15/acquisition_sizes_726.tsv | fd9ade3bb4ebd2ee6c129f69819b9697cf20abd7c10c77a12644bef25396e386 | byte copy. Package sheet Review_v16/Database_sizes/database_sizes_710.tsv (map row ABS-01, 710 FASTAs, maximum 16,923) is the superseded 710-FASTA roster; the maximum agrees (16,923) |
| source_data/Fig1/human_reference.json (owned by Fig1) | LAB_ROOT/MM_review_analysis/A01_4_1_database_sizes/results/human_reference.json | b705dfc4aba8718c83939a812e3de6614ce0a72eb9260c3e0c7b6ff632cf1079 | figures/Fig1/prepare_source_data.py (path reduced to file name) |
| source_data/SFig15/compute_stage_summary.tsv (owned by SFig15) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/reruns_347/A7_other/resource_accounting/fixed347/source_data/compute_stage_summary.tsv | 37c1d94d821d5ba9defb5ff5c8495905ce93209a136952dec8b600b74816498f | byte copy. Package Supplementary_Figure_15/compute_stage_summary.tsv still has 331 CapScan acquisitions |
| source_data/SFig15/compute_task_records.tsv (owned by SFig15) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/reruns_347/A7_other/resource_accounting/fixed347/source_data/compute_task_records.tsv | ea24e72fe85b7e4c0755a849d5822a5347d9d3685e73f765a641072c2d7551b8 | figures/SFig15/prepare_source_data.py (file names replaced by pseudonymous acquisition ids) |
| source_data/Fig2/member_summary.tsv (owned by Fig2) | LAB_ROOT/MM_review_analysis/A14_4_3_campi_gap/results/member_summary.tsv; package Figure2/member_summary.tsv | bdf6d5f0ded8bc8455f95bd1ac8c4ab52af38294ac7030959a269715873591b8 | byte copy |
