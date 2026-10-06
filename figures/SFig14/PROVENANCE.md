# Provenance: Supplementary Figure 14

## Panels
| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a-d (one figure) | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/supplements/revision_2026-09-24/build_revision.py (sha256 e636bd85164896f53ec1d13ef671754ff800597d2d955a3eb01f5151dced1bbe), lines 278-291 and save() lines 34-52 | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/supplements/revision_2026-09-24/delivery/figures/HISTIDINE_MEASURED_SPECTRA.{pdf,png,svg} | pdf 4b027ebdadb153ca2ce3d297d7fdfd5843e4335ded4db21132032f5de31e7a50; png cf3013f9404667db0d069770d69d658edb7e6d2dcd212ae5533a4e27362507c3 | none recorded (run 2026-09-24 14:33, not under Slurm) | CPython 3.10 (user site-packages, from the build log); matplotlib 3.10.8 (PDF metadata); NumPy and pandas versions not recorded |
| composite (title removed) | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/prepare_supplement_figures.py (sha256 5c52b4b89037916a4503e3b77e8454522186a83c0fdff48cced5f1776141ce4d); title span from inspect_package_inputs.py line 16 (sha256 4f20e9d1e468986ff74626a429aebbe4248fef8f0c6a5f82108a0e64d58d3349) | LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/johannes_main_text_review_2026-09-24/codex_review_2026-09-25/release_candidate/supplementary_figures/Supplementary_Figure_10.{pdf,png} | pdf 1239e2c7b56f37bdf0765d5d2aecdb7212ab8db15b93ea5df1af6420a2879ddd; png 4b3e84a7d0791c0841e7866a490c1f8488b17e5ac123f9e1cfe0e0ce90051f37 (= v28d word/media/image22.png) | none recorded | PyMuPDF version not recorded |
| v30 composite, house style (all panels) | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig14/build_sf14_style.py` (sha256 66cd5621...) with `sf_common.py` (sha256 9a1c9929...), run by `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/run_sf456.sh` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig14/SupplementaryFigure14.pdf` and `.png` | pdf f5aa8c121daf384713fd195e1546f682f3a0aff27001304a98ced6493feb25e2; png 0ce9a0c46aa000236b3dbe5fa60e93a7fcc6024915e05d7f4320d5b0c7f14b14 | not recorded (build.log has no job id) | the repository environment `fastalake-paper_v28_work/env` (run_sf456.sh): CPython 3.12.4, matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3, PyMuPDF 1.27.1 |

v30 result: `panels.py:figure_v30` in the repository environment gives a PNG byte-identical to the v30 image (`checks/SFig14_v30.json`), with every plotted value identical to the producer drawing. The v30 build also wrote `values_check.json` and `pdf_audit.json` (fonts Söhne Buch and Halbfett only, text ink #4e4f51 or white, clip check passed); they are not copied here.

Regeneration in this repository: notebook executed headless with `paper.py build SFig14` on
Slurm jobs 6069090, 6069111 and 6069309 (final run after the fix round; p.hpcl92/p.hpcl8/p.hpcl93, 4 cores, 4 GB
requested, peak RSS 0.43 GB measured in 6069111, 4 s wall) on 2026-10-06, with CPython 3.12.4, matplotlib 3.10.7, NumPy 2.2.6, pandas 2.3.3,
PyMuPDF 1.27.1. Receipts in `checks/`: `panel_abcd_png.json` IDENTICAL (180 dpi PNG),
`panel_abcd_pdf.json` IDENTICAL (PDFs at 300 dpi), `SFig14_v28f_png.json` IDENTICAL,
`SFig14_v28f_pdf.json` IDENTICAL. The regenerated composite PNG is byte-identical to the
v28d image (sha256 4b3e84a7...). The regenerated producer-step PNG has sha256
f3723ba3a2fa4eb92ba842c47b4e5d44880b28a8002fdc8d552f69a9a4443e96 in every run. PDF bytes
change from run to run (matplotlib and PyMuPDF write creation metadata), so no PDF sha256 is
recorded; the PDFs are checked by their pixels rasterised at 300 dpi. Before this, the producer had never been rerun; the only
evidence was hash identity of the delivered PNG.

The plotted tables were written by the upstream spectrum extraction
(LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/frozen/search_revision2.py,
sha256 e8593f5b24caccaf7f6affb9126a3ac6f3ec0825b88de7250968cea75db3b7a2, with
fragment_helpers.py d0dca418ddbd95084a920ff5928d06dacd76341cd9530c8f9f33ed8358c4cd4e), Slurm
job 5937598, output LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/.
They were copied by prepare_current_inputs.py (lines 18-21) into the producer's
`source_data/spectra/`. The notebook re-derives processed peaks and matched fragments for all
four cases from the shipped observed peaks with a port of that code and finds them equal to
the shipped tables (exact float values; matched ion counts equal Sage matched_peaks). The
non-specific searches behind the four PSMs ran Sage 0.14.6: the `version` field of each case's
`results.json`, whose sha256 matches the one recorded in the shipped `PROVENANCE.json`.
HIS_001 to HIS_003 come from the A13 non-specific searches
(LAB_ROOT/MM_review_analysis/A13_4_7_predict_unspecific/results/<acquisition>__unspecific/).
HIS_004 (acquisition 20191024_QX7_PeTr_LC02_SA_MicrobPred_44) does not: that mzML was one of
the five truncated PREDICT files, and its PSM comes from the repaired route. The mzML was
reconverted from RAW (LAB_ROOT/MM_review_analysis/recovery_2026-09-15/raw_recovery/20191024_QX7_PeTr_LC02_SA_MicrobPred_44/,
sha256 e677a105bba5e65d06a3e3a3d69c9b2d153a7895a596e71b24cd3c38e24d0119) and searched again by
LAB_ROOT/MM_review_analysis/code/repaired_predict_searches_20260915.py (sha256
b31b2e94fae08ed9a46173917c412e04a8711444edfd9595188682e2fdbba1ef), Slurm array job 5937565
(repaired_search_full.sbatch). Its results.json is
LAB_ROOT/MM_review_analysis/recovery_2026-09-15/search_recovery/20191024_QX7_PeTr_LC02_SA_MicrobPred_44/unspecific/results/20191024_QX7_PeTr_LC02_SA_MicrobPred_44__unspecific/search/results.json,
sha256 7ad4a697b5d4f06720ac3f246ad48dfb30002c69cbb8b878fd9e088d415535e6, version 0.14.6. The
A13 results.json for the same acquisition (truncated input) has a different sha256
(590e7e9c0bbeb5a00966cb4c13a2e9b8c7e6ab488de577a8af6d2636c0673115) and is not used.

## Codex census reconciliation
| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| Supplementary Figure 14 (a-d) | Supplementary Figure 14, word/media/image22.png, sha256 4b3e84a7... | BYTE_MATCH; COPY_ONLY; producing script build_figures_v15.py:36 (copy); selected origin coauthor_package/supplementary_figures/Supplementary_Figure_12.png; "Original plotting/rendering producer unresolved" | rendered by build_revision.py (lines 278-291, save()), title removed and rasterised by prepare_supplement_figures.py, then copied unchanged by build_figures_v15.py line 36 | partly | Byte identity and the copy step agree. This unit resolves the producer Codex left open: rerunning the ported render and title removal gives a byte-identical PNG. Codex also lists release_candidate/Supplementary_Figure_10.png among byte-identical files, which is the title-removal output. No rebuilt candidate for this figure in REBUILT_CANDIDATES.tsv, consistent with no v28f rebuild. |

The Codex statistics register (STATISTICS_REGISTER.tsv, MISSING.md) and number recheck
(out_recheck/REPORT.md) have no row for Supplementary Figure 14; the figure carries no test.

## Source data
| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| PREDICT_spectrum_case_review.tsv | package Supplementary_Figure_14/PREDICT_spectrum_case_review.tsv (identical bytes); producer input LAB_ROOT/00_PACKAGE_2026-09-17/cleanup_2026-09-22/v10_mm_philip_comment_audit_2026-09-22/nine_point_followup/manuscript_final_pass_2026-09-23/mm_v11_revision_2026-09-23/reappraisal_2026-09-24/supplements/revision_2026-09-24/source_data/original_S22/26d1626714ba_independent_case_review.tsv (sha256 26d1626714ba453dc783b0ba9d091075c251dc79777d6f12207f140b518f8517, CRLF line ends, same values). That file is byte-identical to LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/predict_spectrum_review/full/independent_case_review.tsv, written by frozen/audit_predict_spectrum_cases.py there (sha256 ae0edcea5f876f4394d50ab36a775769fad0863759fa99712d8a9176bc1d4320; lines 158 and 160 define accepted_in_semi and N_flanking_H_inferred) in Slurm job 5940796 (full; pilot job 5940795; jobs/full.sbatch, logs/full_5940796.log), with MANIFEST.json (sha256 45bcd329d2f8cb8964919f20c93d885e6fe12c4498181bb9ba54e1ad5f7e08a1) and full/COMPLETE.json (sha256 029d9eecf5c29d88acb7a40e2a7940fa198a8b1f907725c3f009e14c5991b53d; status PASS, 4 cases, 306 acquisitions, scientific_review PENDING_VISUAL_AND_WRITTEN_SYNTHESIS). The script was adapted from the CapScan case audit (MANIFEST source_worker recovery_2026-09-15/spectrum_case_review/frozen/audit_capscan_spectrum_cases_20260915.py, sha256 9b14e329ba3d68535596cb41ca367a0977a79dfa0036421d0357fd2f89b491a5) and reads A25 results_revision2/PREDICT | cc7ffc04213fa85ce2ab1f7a55fd0abab9e5f9c895b6d9ade4e0bc6045fd5ff8 | copied; matched ions, backbone coverage, alternative counts, accepted_in_semi and N_flanking_H_inferred re-derived by panels.py:stat_case_review (all equal) |
| spectra/HIS_001/processed_peaks.tsv | package Supplementary_Figure_14/spectra/HIS_001/ = A25 results_revision2/PREDICT/HIS_001/ (identical bytes) | 524424e9b8c8eea1921768c4717476b7cd171700b4504465ae2d25a1271de64c | copied; re-derived by panels.py:stat_rederive_spectrum, equal |
| spectra/HIS_001/matched_fragments.tsv | same | ad8564088d3e603c0f76428c8d39d0067909e1c5bc2d18f631ed4fbb3e0c3b1a | copied; re-derived, equal |
| spectra/HIS_001/observed_peaks.tsv.gz | same | 24587300e52c1a504dd5251b77966a0e86b1f441f87869825614f09accddce35 | copied unchanged |
| spectra/HIS_001/bounded_sequence_alternatives.tsv | same | faa42f67f6564d18363a3615c4a1e6ff10658c5ea233a87885fc7b19b877cebf | copied unchanged |
| spectra/HIS_001/PROVENANCE.json | package copy (paths reduced to file names); A25 original differs only in the three path fields | c63c8ca7629ad83294ece0416c865628a255eb621c174c29dd51414ecc3033e8 | copied unchanged from the package |
| spectra/HIS_002/* | as HIS_001 | processed 062827e3579abdba1869585faf0ff08e258deb57ca04ba1cab3d54ae27b5e1c7; fragments e52dacaad2807aebd5ef9526eae7cfa34757727384883e763633b7818d32cc4c; observed 38cf0d2b2a52473898e0fc01cf70ecde422e91e80b5e225e0a300e6928158350; alternatives 3780d3648749517cff8f7d750260f2dacd7abc09330adeda090abc09107da7b7; PROVENANCE.json 31894c05387d9fad6332002d4b5ab57cccd589861b3cb045cd852276269a1469 | as HIS_001 |
| spectra/HIS_003/* | as HIS_001 | processed 2062af57461fe3768b9e7c15e58a2d0a1a427ee480ec558b96bd57869d247342; fragments ce847bb980a96fbcc111b7a9513e6931dfdaa381e036abb8771694e323324e08; observed f9ab97b74b14790829b35467627392520d8a5dc226da2c45cbc923dd63914e3d; alternatives cf0fdf3f4be75168b80e391e4db582e21b848587c496c3b849d3e7b462c8d70d; PROVENANCE.json 78cf7eda7012b30900aa66af1fe3dc86c9a7a8490cf569d774895bf7e613d5b4 | as HIS_001 |
| spectra/HIS_004/* | as HIS_001 | processed 0bb92c37371cd69c80e38c678999b84a1df3e3e0b2b884d8e0553cc1ecd12585; fragments d68208d35f0229fb18c177bfc6bd189173d3dbdf429ca015cd2a4ba20a2d4649; observed 4b2741259600b9f5813a629293619e035010b63d172559626a85f3beeacc63df; alternatives 8e9ed927e00b275fb618af2bc1769dcbab787318bbab9f1f09e8911de6051c31; PROVENANCE.json 8725aee59a9f532eef532e33737545026cac26648d0fc71ebdce31d197df29c0 | as HIS_001 |
| spectra/selected_cases.tsv | package Supplementary_Figure_14/spectra/ = A25 results_revision2/PREDICT/ (identical bytes) | fb903709852b39041f92a55beeac3cf899870d27688b78e2baefea275318f696 | copied; rule checked by panels.py:stat_selection |
| spectra/all_strata.tsv | same | 8d1f3263b5219deb52b06e07ecb8dace03091bd9c92fe2847795db290350e328 | copied |
| spectra/complete_candidate_counts.tsv | same | 0280c6d595144d9125385c717e3c35b096fa632fef42adf80387e8fd5497c706 | copied |
| spectra/eligible_case_selection_psms.tsv.gz | same | 0e9870476d7bbd556bd269b0506945b93c920adb80e8b6d047364ff5a0df0a1d | copied |
| spectra/same_scan_all_search_calls.tsv | LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/same_scan_all_search_calls.tsv (not in the package; written by search_revision2.py line 138) | c41eaa2bea0d0297292dca51688b535ead9acf8e46cb0e9e954c8b8e333618b5 | copied; input of accepted_in_semi in panels.py:stat_case_review |
| spectra/spectrum_audit.tsv | LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A25/results_revision2/PREDICT/spectrum_audit.tsv (not in the package) | 158abf414b43037e1560972bb97e31241aff2f85e44480b557327dd678570963 | copied |

SOURCE_DATA_MAP_reanchored.tsv row SF14 (panels a-d, "four spectra", recomputed 4, VERIFIED)
lists the case table and the five per-case files; all are here with identical bytes. The map
does not list the selection tables, which the package also ships; they are included here.
