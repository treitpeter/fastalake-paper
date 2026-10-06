# Provenance: Supplementary Figure 9

## Panels

v30 (current). `FIN` is `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`.

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a, b, v30 house style | `FIN/v30/figures_style/SFig09/build_sfig09_style.py` (de9bc16db3ac352a016275157e62089ab8732e1ad154df7059e6974a4acc65f2) | `FIN/v30/figures_style/SFig09/SupplementaryFigure9.png`; `.pdf` | 360f1f3e69ad7299cdcca3b26ea2d0744d0ff2394ad1b07604a4d3a6c0550fdf; dc579ef30da325e72260ca105e620779feea34966ffb329ca0b62105021a40f5 | srun via `FIN/v30/figures_style/_work_sf07_09/run.sh` (job id not recorded) | `fastalake-paper_v28_work/env` python; fl_style with Söhne |
| a, b at 9-50 with lowered panel-a y minimum (v29) | `FIN/len950/SFig09/build/build_ylim.py` (6e5f73f2f882fa322a87a576e5e4019e325552d24bb9ff05567d409cd569f289), table from `FIN/len950/SFig09/build/build.py` (b2fadeb401a648b139f7c1e239a71df346fff5d9cafc4cb8c38622e26abdb382) | `FIN/len950/SFig09/build/sfig09_ylim_9.png` = `FIN/v29/figures/SFig09/SupplementaryFigure9.png` (`expected/SFig09_v29_reference.png`) | 8f9050ddd3f7ed63005c85f4cc630596f641a2791c714992d018ef74313c77d1 | not recorded | as above |

The repository render is IDENTICAL to the v30 PNG and PDF at 300 dpi, and its producer
reference render equals the v29 PNG pixel for pixel (`values_check.json`). v28d producer
(superseded; code kept as `panels.py:figure_ab`, called inside `figure_v29`):

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a, b (one figure) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/four_author_integration_2026-10-02/actual_FAR_261002/build_JAG121_proposal.py (sha256 24ced16fe55d1e4f9d98714fc502f639c2adaeff76c3948de8243491f8cdceb8), lines 41-62 | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/four_author_integration_2026-10-02/actual_FAR_261002/JAG121_label_proposal/Supplementary_Figure_09_proposed.png (and .pdf, .svg) | png 8f374e3683ed5464b0df038007d8ff7b134e93557c050fd8541b47d51e3e45dd; pdf eb9e3548bd4886ee3edecfbf9be1a00d688a2cafbdea94ecceca5e94ebc79a4d | none recorded (run interactively, 2026-10-02) | not recorded by the producer; the repository environment (CPython 3.12.4, matplotlib 3.10.7, pandas 2.3.3, numpy 2.2.6) reproduces both files pixel for pixel |

Regeneration in this repository: `panels.py:figure_ab`, notebook run on Slurm
(p.hpcl92/p.hpcl8/p.hpcl93, 4 cores, 4 GB, 3 s) on 2026-10-06 with the Söhne Buch TrueType
file. Receipts (v28d check, since removed from `checks/`): `panel_ab_png.json` IDENTICAL (0 of 2,365,200 pixels differ),
`panel_ab_pdf.json` IDENTICAL (PDFs rasterised at 300 dpi), `SFig09_v28f.json` IDENTICAL.
Regenerated files: png sha256 a39b68ddc09bf95f63cd2b30e7326b0f9aae13b95ebc0d704192435ec45f0d2e,
pdf sha256 5bc12bf67d05b34ab04a37d6bda4a928ae5bd69d3d061a7fac7f124dfc10b4f8 (file bytes differ from the
producer's; pixels are identical).

The producer's three input tables were written by
LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/v15_MM_coauthor_2026-09-25/build_figures_v15.py
(sha256 e608bba2d87801bf06ad7b7bc406ede9bc76d4d90a3f1cd9cf2fb8f2ab8e6729, lines 50-54) from
the 282-row MM6 join. `panels.py:derive_tables` ports those lines; the notebook re-derives
the three tables and finds them equal to the shipped copies (row sets, tertile labels and
bounds, n identical; gains within 4.6e-13).

## Codex census reconciliation
| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| Supplementary Figure 9 (a, b) | Supplementary Figure 9, word/media/image17.png, sha256 8f374e36... | BYTE_MATCH; WRITE_SITE_HASH_RECEIPT; producing script build_JAG121_proposal.py:58; copy apply_JAG121.py:37; receipt PROPOSAL.json | build_JAG121_proposal.py renders (line 58 is the savefig); apply_JAG121.py copies unchanged (line 36) and embeds (line 61) | yes | Codex did not rerun the producer; this unit reran the ported code and obtained pixel-identical PNG and PDF. No rebuilt candidate is listed for this figure in REBUILT_CANDIDATES.tsv, consistent with no v28f rebuild. |

Codex number recheck (out_recheck/REPORT.md, rows N0118, N0119, N0211) reports MATCH for
80.2/68.6/50.4, 76/75/75, 250 and 80/69/50; this unit agrees. The Codex statistics register
lists no row for Supplementary Figure 9 or for the line-129 descriptive values.

## Source data
| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| capscan_acquisitions_9_50.tsv | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/len950/SFig09/capscan_acquisitions_9_50.tsv` (written by `len950/SFig09/build/build.py`; no package sheet yet) | 5bb9ec898849e853e93ef0d01da9b7eaebab64db2b1a075b89f8d9d44fe5ad39 | copied unchanged; rebuilt from `capscan_assembly_dna_yield_gain_by_acquisition.tsv` and `a08_benchmark282_old.tsv` by `panels.py:acquisitions_9_50` and equal (notebook) |
| capsule_assembly_gain_tertiles_9_50.tsv | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/len950/SFig09/build/tertiles_9.tsv` | 761b3e620045ad96c4cdcaefb49ca4eb3bc6b8eeefb48c3e4e4b2719b8f2638f | copied unchanged; used as the expected panel b values in `render_v30` |
| capsule_assembly_gain_per_specimen.tsv | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/source_data/current/v15/capsule_assembly_gain_per_specimen.tsv; package sheet current/v15/capsule_assembly_gain_per_specimen.tsv (identical bytes) | 800c0e370c61f8557539c1ec149d981965c49c893797910301108ca49cc144e9 | copied; re-derived by panels.py:derive_tables and equal |
| capsule_assembly_gain_acquisitions.tsv | same folder; package sheet current/v15/capsule_assembly_gain_acquisitions.tsv (identical bytes) | 89956d96f81cc50656da184c078c6f14b50e44381b4c1fe3131eb7eca4ddc1c9 | copied; re-derived and equal |
| capsule_assembly_gain_tertiles.tsv | same folder; package sheet current/v15/capsule_assembly_gain_tertiles.tsv (identical bytes) | 9d9b3e2bf774245e74406e5111e69ad02654eebe4eea746641c48b4c71144796 | copied; re-derived and equal |
| capscan_assembly_dna_yield_gain_by_acquisition.tsv | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/main_methods_pass_2026-09-24/review_decisions/MM6_capsule_metagenome_2026-09-25/capscan_assembly_dna_yield_gain_by_acquisition.tsv; package sheet Figure4/capscan_assembly_dna_yield/capscan_assembly_dna_yield_gain_by_acquisition.tsv (identical bytes) | 351bdf31dc30f108e11796f15614cb05299ba8d37a8991d6ca4f10f7fa1a9c4d | copied unchanged |
| CAPSCAN_sanitize_report.tsv | LAB_ROOT/canonical/publication_data/sources/assemblies_prefixed/CAPSCAN_sanitize_report.tsv (not in the package) | 0120dd064a2b20ab39c59177519cb637aeb701bbcd392ca691260b630a9f7205 | copied unchanged |
| a08_benchmark282_old.tsv | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/public347/collect/a08_benchmark282_old.tsv (not in the package) | 4c732679b085a01245d39f324d0884d2eff1986f73914bdd9f37c7ddec2eee44 | copied unchanged |
| a08_benchmark282_new.tsv | same folder, a08_benchmark282_new.tsv (not in the package) | fa6eaf0fddd45be5a5dc989c5c64397421a82713beb0fd084a7dc9aea1700265 | copied unchanged |
| capscan_roster_347.tsv | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/roster/fixed_roster_347.tsv (not in the package) | 473ff6afcf7348163f1fbaa8d624a8078eed9348a15019952363bc86e70f7aa5 | copied unchanged (renamed) |

SOURCE_DATA_MAP_reanchored.tsv rows SF9-a (250), SF9-b (tertile medians) and RES-74
(80.1659) point to the current/v15 sheets above; the values recomputed here agree. The
package has no Supplementary_Figure_09 folder; the figure's data sit under current/v15/ and
Figure4/capscan_assembly_dna_yield/.
