# Provenance: methods

## Panels

This unit has no figure panels.

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| none | | | | | |

## Statistics: original producers

| statistic | producer (LAB_ROOT/...) | output (LAB_ROOT/...) | sha256 of output | Slurm job id(s) | interpreter |
|---|---|---|---|---|---|
| line 215 reservoirs | canonical/publication/chain/00_lake.sbatch (lake_builder, 64 threads) | MM_review_analysis/A01_4_1_database_sizes/inputs/{HSTOOL,PREDICT,CAPSCAN}_reservoir_stats.txt; canonical/publication_data/arms/uniprotkb_2026_03_noclust_2026-09-09/cohorts/CAMPI_MOCK/0_lake/CAMPI_MOCK_lake_stats.txt | 4a810dca...abbe75, 2fc2bfbb...d555, abc4d914...a7a, d6f8d57e...f7fd8 | CAMPI 5895248 | Rust lake_builder, hash version 5.3 |
| line 227 comparator FASTAs | supplied by the CAMPI study; MG duplicate-accession records removed before freezing (METAGENOME_COVERAGE_VERIFIED.json: 19,551 to 19,435, sequence sets identical) | canonical/publication_data/arms/uniprotkb_2026_03_noclust_campi_ref_append_2026-09-09/inputs/campi_ref_targets.faa; .../uniprotkb_2026_03_noclust_campi_ref_metag_append_2026-09-09/inputs/campi_mo_targets_fixed.faa | 9d4def2d...fce9; 364c6415...c8d0 | none | none |
| lines 233, 303 reconversion | review/.../final/completeness_653/code/sage_input_integrity.py; review/.../final/second_engine_all/truncation_check/check.sh | review/.../final/completeness_653/out/search_coverage.tsv; review/.../final/second_engine_all/truncation_check/predict_mzml_endings.tsv | 1b75c166...0b5d; 171632df...5755 | 6048085 (old), completeness_653 jobs; 6068722 | not recorded |
| line 267 split bootstrap | review/.../claude_review_2026-10-03/iteration3/disjoint/analyze.py (sha256 1b436e1e...b843) | review/.../iteration3/disjoint/SUMMARY.tsv | 71c6dd9c...bc7a | see iteration3/disjoint/SUBMISSIONS.jsonl | not recorded |
| line 267 cross-fit roster | review/.../iteration4/crossfit/run/analyze_final.py; rerun final/public347/crossfit_capscan347/run/analyze_final.py | review/.../iteration4/crossfit/run/PER_ACQUISITION.tsv; final/public347/crossfit_capscan347/run/PER_ACQUISITION.tsv | f7616394... (prefix, figures/Fig4/prepare_source_data.py); 52f4fa61...e5c3 | 6054800, 6054853; final 6065360 | worker: anaconda 3/2023.03 (Python 3.10.9, NumPy 1.26.4, read by prepare_source_data.py) |
| line 269 unions | sample_specific_merge_all_raw/2026-09-11/code/workflow.py build_union (sha256 8dc29123...6f4f); final/reruns_347/A7_other/union/scripts/run_union.py | sample_specific_merge_all_raw/2026-09-11/cohorts/{HSTOOL,PREDICT,CAPSCAN}_noclust/union_summary.json; final/reruns_347/A7_other/union/fixed347/union_summary.json | af1e3408...5b, 66202a1e...2f, 0162239d...b7, 887144c5...8f | 6064861, 6064862 (347 rerun) | not recorded |
| line 275 antibiotic | review/.../manuscript_completion_2026-10-01/clinical_refit.py (sha256 b4960505...5072) derive() | review/.../manuscript_completion_2026-10-01/clinical_refit/clinical_paired_statistics.tsv | 9518b011...c3bc7 | 6042666 | Python 3.12.4, NumPy 2.2.6, pandas 2.3.3, SciPy 1.16.3 |
| line 299 compute | Slurm accounting (sacct) of jobs 5895248, 5895250, CAMPI chain arrays and array 5450700 | as shipped in source_data/SFig15 (common_resource_accounting) | see figures/SFig15/PROVENANCE.md | 5895248, 5895250, 5450700 | none |

Full input and output sha256 values of the extraction are in the work directory receipt
(`prep_receipt.json`, job `methods_prep`) and repeated in `source_data/methods/README.md`.

## Codex census reconciliation

The Codex figure hash census covers embedded figure images only. This unit has no images.

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| none | none | none | none | NA | text statistics only |

Other Codex evidence compared (results, not code):

| Codex item | Codex result | this unit | agree? |
|---|---|---|---|
| STATISTICS_REGISTER S034-S039 (split bootstrap) | 10,000 draws; seed and software missing | 10,000 draws, seed 20261003, percentile, participant clusters (from analyze.py); all 40 summary rows reproduced | yes; seed filled |
| STATISTICS_REGISTER S055-S056 (antibiotic) | 200,000 sign flips, 50,000 bootstrap, seed 20260907; statistic not reported | same; statistic = mean difference with 67 and 250 exceedances of 200,000; sign-flip seeds 20260907 and 20260908 | yes; statistic and second seed filled |
| out_recheck REPORT N0242 | 578, 296, 282 acquisitions MATCH | 578 (296 PREDICT, 282 CapScan) | yes |
| out_recheck REPORT N0121 | solver perturbation counts 25 and 27 of 578 (text 28) | not in this unit (Supplementary Figure 6 legend) | NA |
| inventory gap: CapScan union | 452,255 is the 331 roster; 460,977 on 347 | 452,255 (331) and 460,977 (347) recounted independently | yes |

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| reservoir_sources.tsv | the four lake-builder statistics files above | as above | `prepare_source_data.py:reservoir_sources` (parse SOURCES BY PRIORITY). No package sheet. Values agree with the reservoir_entries of `source_data/Fig1/acquisition_sizes_726.tsv`. |
| campi_comparator_records.tsv | the two frozen comparator FASTAs above | 9d4def2d...fce9; 364c6415...c8d0 | `prepare_source_data.py:comparator_records` (one row per record; sequence replaced by its SHA-256). No package sheet. |
| union_contributions.tsv | 695 selected FASTAs listed in sample_specific_merge_all_raw/2026-09-11/manifest.json (sha256 9a9839d3...e159) plus 16 participant-15 FASTAs (711 files; 1,042 table rows because the 331 CapScan FASTAs appear in both CapScan rosters) (participant15_added.tsv, sha256 688dea12...559c) | per FASTA in the table | `prepare_source_data.py:union_contributions` (independent recount). Equals input_records and unique_sequences of all four union_summary.json files. Package: Figure3/common_fasta_union_summary.json (healthy stool 153,267) agrees. |
| crossfit_lda_fallback_347.tsv | review/.../final/public347/crossfit_capscan347/run/PER_ACQUISITION.tsv | 52f4fa61...e5c3 | `prepare_source_data.py:crossfit_fallback` (columns i, cohort, lda_fallback). Package sheet Supplementary_Figure_06/metagenome_crossfit/lda_fallback_acquisitions.tsv lists the same three roster indices (168, 185, 496). |
| split_spectra_mg_per_acquisition.tsv | package sheet Supplementary_Figure_06/metagenome_same_spectra/per_acquisition.tsv | 7ba4e68d...93d5 | `prepare_source_data.py:split_mg` (drops `sample`; blanks PREDICT `group`, a specimen-visit code). All other values identical to the sheet. |
| split_spectra_mg_summary_original.tsv | package sheet Supplementary_Figure_06/metagenome_same_spectra/summary.tsv | 71c6dd9c...bc7a | byte-identical copy |
| predict_reconversion.tsv | review/.../final/completeness_653/out/search_coverage.tsv; review/.../final/second_engine_all/truncation_check/predict_mzml_endings.tsv | 1b75c166...0b5d; 171632df...5755 | `prepare_source_data.py:predict_reconversion` (PREDICT primary_public_current rows; path columns dropped; merged with original end-of-file check). No package sheet. |
| compute_nodes.tsv | `scontrol show node` on 2026-10-06 | NA | `prepare_source_data.py:compute_nodes` |
| crossfit_interpreter.tsv | interpreter named in review/.../iteration4/crossfit/design/run.sbatch | NA | `prepare_source_data.py:crossfit_interpreter` |
| (shared) source_data/Fig4/crossfit_per_acquisition.tsv, clinical_paired_statistics.tsv | see figures/Fig4/PROVENANCE.md | 5bcdbe88...b094; 9518b011...c3bc7 | read only |
| (shared) source_data/SFig15/compute_task_records.tsv, compute_stage_summary.tsv, campi_alphanovo_inference_sacct.tsv | see figures/SFig15/PROVENANCE.md | fb248054...c2a; 37c1d94d...98f; d007807d...029 | read only |
| (shared) source_data/Fig1/acquisition_sizes_726.tsv | see figures/Fig1/PROVENANCE.md | fd9ade3b...e386 | read only (cross-check) |
| controlled paired_antibiotic_profiles.tsv | review/.../manuscript_completion_2026-10-01/clinical_refit/paired_antibiotic_profiles.tsv | b99dd07a5cf75150f7cb6520eedcc816cefcd6e562dea9b457c309dfa95f27ec | copied unchanged to the controlled directory, never into the repository |

`review/...` abbreviates `review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03` (and its parent `manuscript_completion_2026-10-01` where stated).
