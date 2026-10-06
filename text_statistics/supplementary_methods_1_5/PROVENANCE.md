# Provenance: supplementary_methods_1_5

`AF` below stands for `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup`, and `CR` for `AF/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03`. Every origin file and its full sha256 are listed in `prepare_receipt.json` (98 files, written by `prepare_source_data.py`).

## Panels

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| none | text-only unit | | | | |

## Statistics

| statistic | producer script | producer output | sha256 | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| AlphaNovo permuted control (553) | CR/final/alphanovo_accuracy/accuracy.py | CR/final/alphanovo_accuracy/results/per_acquisition.tsv | 82bc0030c9cafaa8b68ea79d00037dc2f44e62374f1801cbb361eeb95ad9a18c | 6064877 | CPython 3.10.9 (anaconda 2023.03) |
| bag-aware lookup (561, 579) | LAB_ROOT/MM_review_analysis/A04_4_2_expansion/alphanovo_bags_2026-09-23/bounded_sweep_2026-09-25/run_bounded.py; AF/bags_prominence_review_2026-09-30/audit_and_plot.py | Source Data package Bag_aware_lookup/bounded48_per_acquisition.tsv (copied by SFig02) | dc48dee627b3687a0ee940cbfdf93e205ac3ab33059a2186b7740a75f32507c9 | see figures/SFig02/PROVENANCE.md | see SFig02 |
| 58-acquisition study space (567) | LAB_ROOT/canonical/publication chain stage 02_evidence | Source Data package Supplementary_Figure_02/reference58/stage1_stats.json (copied by Fig1) | 2be7eafbd02937b6df5ac05e1845c1a4d1d64203bfc62bd77ebb23739a21fcd8 | 5872582-5872587 (per Fig1) | see Fig1 |
| 97-acquisition study space (567) | LAB_ROOT/canonical/publication chain stage 02_evidence (production run of 2026-08-02) | LAB_ROOT/canonical/publication_data/cohorts/HSTOOL/1_evidence/evidence_lake.stats.json | d96e655e269271c2b0c92bf6d6467db4069a49fe09e170fcaa244e82f4f70d16 | not recorded | not recorded |
| parser audit (573) | CR/final/reruns_347/A7_other/parser_audit/fixed726/collect.py and the copied Rust binary m07_parser_audit | CR/final/reruns_347/A7_other/parser_audit/fixed726/ACQUISITION_SUMMARY.tsv | 28c39e5eb91f4472730e2777fcc12b96b8a3641a482f13dc193698b162443248 | 6064887, 6064888 (original 710-file audit 6000982) | Rust binary as frozen in M07 |
| CAMPI target FASTAs (609) | supplied by the CAMPI study | LAB_ROOT/canonical/publication_data/arms/campi_db/campi_ref_targets.faa; campi_mo_targets.faa | 9d4def2d1de49d61e88f163e8c456da891b407d77e8959f58da9dbead6d9fce9; b64f27cbc9d5a7a54670252a7021895878d28d929dcf626fa59638c61e663bd0 | none | none |
| member presence (611) | LAB_ROOT/MM_review_analysis/R8_open_requests_2026-09-17/B_blind_reservoir_overlap/overlap.py | (copied by Fig2) | 032c84a3383dbd16ecccb95d586c3fc3b2f131d95c4f1757d4ab0983005c68a8 | 5954792 | see Fig2 |
| input completeness (629) | CR/final/completeness_653/code/sage_input_integrity.py, run653.py | CR/final/completeness_653/out/search_coverage.tsv | 1b75c1666a0795d5660af57ae01f2547366ae57db35f1efa49b26ca6be3dbf74 | 6065843 (original check 6048085) | see methods unit |
| partial-origin human share (651) | CR/iteration1/verify_proposals_methods/partial_human.py; CR/final/reruns_347/A7_other/partial_origin_share/fixed347/partial_human.py | logs baseline331.out and fixed347.out of the rerun (values 0.120, 0.604 to 0.613, 1.964) | taxonomy_rows.json and taxonomy_evidence_signal.npz per cohort, see prepare_receipt.json | 6064869, 6064870 | CPython 3.10 (miniconda), numpy, scipy |
| seed-ortholog check (657) | CR/iteration1/fix_eggnog_native/scripts/prep_native.py, verify_shard.py, summarize_verification.py | CR/iteration1/fix_eggnog_native/verify/native_row_verification.tsv.gz; VERIFICATION_SUMMARY.json | a64fb0da67c402288a531aa785d5c6584e641543e759f27a9cd29382b9e03c4c; 4043d52873190f2374f9afebf2b7ae6b85eb4664dedbc4c8d549d51ee9dfae12 | 6048172, 6048175, 6048176, 6048298 | DIAMOND 2.0.11 (per text) |
| MGnify coverage 88-95% (657) | AF/github_workflow_reconciliation_2026-09-30/resumed_audit/annotation_repair/evidence_repaired/repair_projection.py, collect.py | .../evidence_repaired/full/<arm>/COMPLETE.json (8 arms); CR/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig5/ra/annotation_repair/evidence_repaired/full/CAPSCAN_noclust/COMPLETE.json | in prepare_receipt.json | 6034799, 6034800, 6034819 and the CAMPI arm jobs in full_collection/COMPLETE.json | not recorded |
| CAMPI supplementary annotation (657) | AF/github_workflow_reconciliation_2026-09-30/resumed_audit/annotation_repair/prepare_supplement.py, run_supplement.py | .../supplemental_native_extract/native_sequences.tsv.gz, native_functions.tsv.gz | e4ea31a71071de91899ca745e1c232c62893de1cd13bca28cbd40f5b25b7d640; e5a5006641691f3df75a62dc8631186a8533d22df8eaee9b65f73dd73a605cda | 6034650, 6034651 | eggNOG-mapper 2.1.12 |
| allocation sensitivity (657) | CR/iteration1/proposals_methods/alloc_sensitivity/alloc_sensitivity.py (sha256 f62592df0f5200a6c7fa0d048a68e3ab110e11e1153123d4056d519a2855c490); CR/final/reruns_347/A3_allocation_sensitivity/code/alloc_capscan347.py, collect.py | .../alloc_sensitivity/out/*.tsv; .../A3_allocation_sensitivity/run347/*.tsv, allocated_fraction_medians_347.tsv, capscan_own_stool_331_vs_347.tsv | see Source data table | 6048205, 6048208 (HSTOOL, PREDICT, CapScan 331); 6064914, 6064915, 6064920, 6064925, 6064926 (CapScan 347) | CPython 3.12.4, numpy 2.2.6, scipy 1.16.3, pandas 2.3.3 (miniconda) |

Recomputation in this repository: `prepare_source_data.py` (Slurm 6069195, 6069200, 6069240; peak RSS 2.0 GB, 40 s) and `stats.py --write` (Slurm 6069232), CPython 3.12.4, numpy 2.2.6, scipy 1.16.3, pandas 2.3.3. Every recomputed test and summary agrees with the producer's own output within 2.3e-16 (`stats.py:check_against_producers`).

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| none | none | | | | The Codex census covers the 23 embedded figure images; none belongs to Supplementary Methods 1-5. The Codex number recheck (out_recheck/REPORT.md) has no row for a number of this unit. The Codex statistics register rows S107-S109 (line 657) are filled from the code in STATS.tsv; see OPEN_ITEMS.md. |

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| parser_audit_acquisitions_726.tsv | CR/final/reruns_347/A7_other/parser_audit/fixed726/ACQUISITION_SUMMARY.tsv. Package sheet Methods_records/M07_parser_audit/ACQUISITION_SUMMARY.tsv holds the 710-file version: its rows equal rows 0-709 here in every column except run time. | 28c39e5eb91f4472730e2777fcc12b96b8a3641a482f13dc193698b162443248 | byte copy (`prepare_source_data.py:parser_audit`) |
| search_input_completeness.tsv | CR/final/completeness_653/out/search_coverage.tsv; not in the package | 1b75c1666a0795d5660af57ae01f2547366ae57db35f1efa49b26ca6be3dbf74 | rows of lineages primary_public_current and healthy_stool_58; path columns mzml_cfg and dir removed (`completeness`) |
| hstool97_stage1_stats.tsv | LAB_ROOT/canonical/publication_data/cohorts/HSTOOL/1_evidence/evidence_lake.stats.json; not in the package | d96e655e269271c2b0c92bf6d6467db4069a49fe09e170fcaa244e82f4f70d16 | key-value table; `database` reduced to the file name (`hstool97`) |
| hstool97_pooled_prediction_files.tsv | directory listing of LAB_ROOT/canonical/publication_data/cohorts/HSTOOL/1_evidence/_pooled_csvs | listing | file names and sizes (`hstool97`) |
| origin_rows_HSTOOL.tsv, origin_rows_PREDICT.tsv, origin_rows_CAPSCAN331.tsv | AF/github_workflow_reconciliation_2026-09-30/resumed_audit/annotation_repair/annotated_origin_consensus_2026-10-01/full/<arm>/taxonomy_rows.json and taxonomy_evidence_signal.npz; not in the package | prepare_receipt.json | one row per origin class, signal summed over acquisitions; asserts the producer's printed 0.120, 1.964, 0.604 (`origin_rows`) |
| origin_rows_CAPSCAN347.tsv | CR/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig5/ra/annotation_repair/annotated_origin_consensus_2026-10-01/full/CAPSCAN_noclust/ | prepare_receipt.json | as above; asserts 0.613 |
| mgnify_seed_ortholog_check_by_genome.tsv | CR/iteration1/fix_eggnog_native/verify/native_row_verification.tsv.gz = package sheet Annotation_verification/native_row_verification.tsv.gz (same sha256) | a64fb0da67c402288a531aa785d5c6584e641543e759f27a9cd29382b9e03c4c | per-genome counts; totals asserted equal to VERIFICATION_SUMMARY.json = package Annotation_verification/verification_summary.json (`mgnify`) |
| annotation_mgnify_coverage_by_arm.tsv | evidence_repaired/full/<arm>/COMPLETE.json (8 arms) and the C2 347-roster CapScan COMPLETE.json; not in the package | prepare_receipt.json | counts copied (`annotation_coverage`) |
| campi_supplemental_annotation_queries.tsv | supplemental_native_extract/native_sequences.tsv.gz and native_functions.tsv.gz; not in the package | e4ea31a7...; e5a50066... | accession, I/L sequence hash, flags; source paths and original headers dropped (`campi_supplement`) |
| campi_target_fasta_records.tsv | campi_db/campi_ref_targets.faa and campi_mo_targets.faa; not in the package. text_statistics/methods uses the frozen append copies (campi_comparator_records.tsv), where the MG file has 19,435 records | 9d4def2d...; b64f27cb... | one row per record, sequence replaced by its sha256 (`campi_targets`) |
| allocation_coverage_HSTOOL.tsv, allocation_coverage_PREDICT.tsv | CR/iteration1/proposals_methods/alloc_sensitivity/out/profiles/<arm>/coverage.tsv; not in the package | 2441b850...; 70085379... | byte copy (`alloc_coverage`) |
| allocation_coverage_CAPSCAN347.tsv | CR/final/reruns_347/A3_allocation_sensitivity/run347/profiles/CAPSCAN_noclust/coverage.tsv | 2c04660ce6631fe6... | byte copy |
| allocation_fraction_medians_producer.tsv | CR/final/reruns_347/A3_allocation_sensitivity/allocated_fraction_medians_347.tsv | e18fa519328dfa49... | byte copy (comparison only) |
| allocation_capscan_own_stool_producer.tsv | CR/final/reruns_347/A3_allocation_sensitivity/capscan_own_stool_331_vs_347.tsv | f220676ce7c43423... | byte copy (comparison only) |
| allocation_predict_crossomics_producer.tsv | CR/iteration1/proposals_methods/alloc_sensitivity/out/predict_crossomics_protein_rows.tsv | a26ed8de60f2923d... | byte copy (comparison only) |
| allocation_hstool_ko_precision_producer.tsv | CR/iteration1/proposals_methods/alloc_sensitivity/out/hstool_KO_precision.tsv | aba3ca25544b5003... | byte copy (comparison only) |
| allocation_capscan_own_stool_distances_347.tsv | run347 OG_other_unresolved profiles and C2 roster/design_347.tsv | prepare_receipt.json | producer code of alloc_capscan347.py up to the distance matrix (`own_stool_distances`) |
| allocation_predict_distances.tsv | alloc_sensitivity/out/profiles/PREDICT_noclust, LAB_ROOT/analyses/enterotypes_2026-09-07 (metaG, metaT, metabolites, 77-patient roster), AF/manuscript_completion_2026-10-01/biology_refit/groups/PREDICT (human-assigned groups) | prepare_receipt.json | check_prevalence.py:transform then scipy pdist; patients replaced by pseudonyms U01-U77 in roster order (`predict_distances`) |
| allocation_hstool_ko_pairs.tsv | alloc_sensitivity/out/profiles/HSTOOL_noclust, AF/manuscript_completion_2026-10-01/hstool_precision_refit/design.tsv | prepare_receipt.json | hstool_precision.py:pairs; medians asserted equal to the producer table (`hstool_pairs`) |

Shared tables read from other units (not duplicated): source_data/Fig2/alphanovo_accuracy_per_acquisition.tsv (equal to the producer output, sha256 82bc0030...), source_data/Fig2/reservoir_protein_presence_by_proteome.tsv, source_data/SFig02/bags_bounded48_per_acquisition.tsv (equal to package Bag_aware_lookup/bounded48_per_acquisition.tsv), source_data/SFig02/bags_bounded48_orderings.tsv, source_data/Fig1/stage1_stats.json.
