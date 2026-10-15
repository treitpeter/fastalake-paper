# Open items: supplementary_methods_1_5

## Manuscript sentences that change for v28f (347 CapScan acquisitions, 15 participants, 726 prediction files)

- MISMATCH: line 573 "all 710 retained prediction files (112,835,982 rows)" becomes "all 726 retained prediction files (114,971,345 rows)"; parser results are unchanged (all 726 files identical, still three lookup-only queries) (NUMBERS.tsv rows for line 573).
- MISMATCH: line 629 "the 637 public-sequence CapScan and PREDICT searches" becomes "the 653 ..." (306 PREDICT + 347 CapScan); lowest coverage stays 99.6% (NUMBERS.tsv line 629).
- MISMATCH: line 651 "0.60% (CapScan)" becomes "0.61%" (0.6128% on 347 acquisitions; 0.6037% on 331) (origin_rows_CAPSCAN347.tsv).
- MISMATCH: line 657 "All fourteen CapScan participants" becomes "All fifteen"; P stays at 1 × 10⁻⁵ in every rule and prevalence filter (allocation_capscan_own_stool_distances_347.tsv).

## Text and code disagree

- MISMATCH: line 629 states that the highest scored scan "had to reach at least 95% of the scan range"; the original check (sage_input_integrity.py, flag at coverage < 0.9) used 90%. Every one of the 711 searches passes both thresholds (lowest 99.63%), so the result holds as written; the text should either say 90% or say that all searches exceeded 95% (stats.py:stat_completeness).
- MISMATCH: line 609 "19,435 MG records". The MG FASTA of the standalone MG search (campi_db/campi_mo_targets.faa) holds 19,551 records: each of the 116 contaminants appears twice with the same name and sequence. 19,435 is the number of distinct records, and the record count of the deduplicated copy used for the FL + Reference + MG searches (text_statistics/methods, campi_comparator_records.tsv). The sentence is correct for the appended database; the standalone search carried 116 duplicate records (campi_target_fasta_records.tsv).

## Statistical details the text does not state (filled from the code in STATS.tsv)

- TEXT-GAP: line 553, permuted-prediction control: one permutation per acquisition (Python random.Random, seed 20261005 plus the character-code sum of the acquisition name); pooled 45 exact matches of 801,613 reference identifications (0.0056%); per acquisition 0 to 0.017%.
- TEXT-GAP: line 657, CapScan own-stool test (Codex S107): two-sided permutation test of the mean own-stool advantage (mean Hellinger distance of a participant's four intestinal-region OG profiles to other participants' stool minus to own stool), 99,999 random permutations of stool labels across participants, seed 20260910, P = (1 + count) / 100,000; n = 15 participants (v28f). P = 1 × 10⁻⁵ is the smallest value this test can return. Twelve tests (4 rules x 3 prevalence filters 0, 0.1, 0.2), no multiplicity correction. Acquisitions are averaged within specimen and specimens within participant-region before the test.
- TEXT-GAP: line 657, PREDICT cross-omics sensitivity (Codex S108): Mantel test with Spearman correlation of Euclidean distances between Hellinger-transformed profiles (metabolites log10, features present in at least 70% of patients, median-imputed, z-scored), two-sided, 1,999 permutations of patients, seed 20260907, n = 77 patients; no multiplicity correction in the sensitivity (the production analysis applies Benjamini-Hochberg within families of 28). The 0.16-0.64 range covers the three restricted rules and the DNA, RNA and metabolite layers over KO, OG and pathway levels (27 tests); all 60 rows are in STATS.tsv.
- TEXT-GAP: line 657, healthy-stool KO correlations (Codex S109): median over acquisition pairs (66 technical-replicate pairs, 45 inter-individual pairs) of the Pearson r of log2 KO intensity over KOs positive in both acquisitions; descriptive. The printed ranges 0.89-0.91 and 0.69-0.73 pool the profile with and without human-assigned signal; with one profile alone the ranges are 0.90-0.91 and 0.70-0.73 (all) or 0.89-0.90 and 0.69-0.73 (human-assigned withheld).
- TEXT-GAP: line 657, "allocated a median of 7-61% instead of 62-81%": medians over acquisitions of allocated / total LFQ intensity in the profile with human-assigned signal withheld, per cohort (58 healthy stool, 306 PREDICT, 346 non-empty CapScan), level and rule.
- TEXT-GAP: line 657, "88-95% of compatible proteins": share of compatible proteins whose MGnify annotation was joined by gene identifier with a verified sequence, before the seed-ortholog check, over the five CAMPI arms and the three cohorts (88.3% PREDICT to 94.7% CAMPI reference; CapScan 89.3% on 331 and 89.2% on 347 acquisitions).

## Not traced or not checked here

- UNTRACED: line 567, 26,881,038 is read from the statistics file of the original 97-acquisition production run (evidence_lake.stats.json); the Slurm job id of that run was not found and the lookup was not re-run.
- UNTRACED: line 567, "Re-extraction using only the 58 retained acquisitions produced identical acquisition-specific selected FASTAs for all 58 acquisitions" is not checked by this unit.
- UNTRACED: line 573, "Three lookup sequences absent from selection evidence had no matches in the saved candidate FASTAs": the count of three is recomputed; the absence of matches is not.
- UNTRACED: line 657, "Most failures involved proteins shorter than half the seed ortholog or a seed shared with the adjacent gene" prints no number; recomputed shares are 70.5% and 42.9% of failures against 2.6% and 2.1% of passing rows (stats.py:stat_annotations).
- AUTHOR: line 657 says the sensitivity values were calculated before the seed-ortholog check. The HSTOOL and PREDICT values come from the 2026-10-03 run on the 2026-10-01 annotated-origin consensus; the CapScan 347 rerun used the participant-15 consensus built with the same 2026-10-01 rules. That this consensus also predates the check was not verified here.
- UNTRACED: numbers in lines 493-659 outside this unit's inventory are not checked here and need an owner: 296 and 335 assembly FASTAs (499, 527, 547), 21 and 15 CAMPI acquisitions and the healthy-stool design counts 12, 15, 10, 11, 10 and 39 (501), PREDICT 306, 302, 93, 296, 91, 42, 85, 77 and CapScan 347, 308, fifteen, sixteen, fourteen (501, 503), acquisition parameters (509, 517, 523), 331 RAW files (523), reservoir sizes 233,497,702, 272,919,212, 264,925,283 and 375,333,654 (545-549, 611), AlphaNovo preprocessing constants (553), filtering constants (557), 13,512 observations, 60 acquisitions and two seeds (561, Supplementary Table 6), 215, 193 and 22 fragment pairs (581), 296 and 282 matched acquisitions (617, 619, 633, 643, 647), five truncated PREDICT files (625), and 58, 306, 347, 15 and 346 grouping rosters (643).
- AUTHOR: the Source Data package does not contain the AlphaNovo accuracy table, the partial-origin rows, the input-completeness table, the member-presence table or the allocation-sensitivity tables, and its Methods_records/M07_parser_audit sheets hold the 710-file version. The tables in source_data/supplementary_methods_1_5/ (and Fig2, SFig02) can fill these sheets.

## Controlled inputs

- CONTROLLED: PREDICT_fullstack_primary_samples.tsv (LAB_ROOT/analyses/enterotypes_2026-09-07) | index (acquisition-derived sample id) and idpatient | 5e967b0aa182f2dfd6b2d902d2d1b393dd96caf43ef619575ee5c5b972aa09a7 | prepare_source_data.py:predict_distances (not the notebook) | defines the 77-patient roster and its order; the repository table carries only pseudonyms U01-U77 in that order.

## Zenodo

- ZENODO: native_row_verification.tsv.gz | 249,402,432 bytes | a64fb0da67c402288a531aa785d5c6584e641543e759f27a9cd29382b9e03c4c | supplementary_methods_1_5 | 9,368,708-row MGnify seed-ortholog verification (already in the Source Data package as Annotation_verification/native_row_verification.tsv.gz) | lets a reader recompute mgnify_seed_ortholog_check_by_genome.tsv.
- ZENODO: allocation-sensitivity profile matrices (HSTOOL, PREDICT: iteration1 alloc_sensitivity out/profiles; CapScan 347: A3_allocation_sensitivity run347/profiles; npz plus axes json, up to 42 MB each; hashes in prepare_receipt.json) | supplementary_methods_1_5 | lets a reader re-derive the distance and pair tables with prepare_source_data.py.

## Shared requests

- PIPELINE: the task did not allow writing in pipelines/. Upstream ids created by this unit, all with pipeline_dir null: healthy_stool_97_study_lookup, prediction_parser_audit, origin_annotation_share, mgnify_seed_ortholog_check, annotation_repair_coverage, campi_supplemental_eggnog, functional_allocation_sensitivity. Please check them against ids chosen by other units at the same time. Reused ids: campi_alphanovo_accuracy, bag_mass_lookup_search, healthy_stool_stage1_lookup, campi_reservoir_containment, sage_input_completeness, predict_crossomics_refit, healthy_stool_ko_precision_refit.
- PIPELINE: functional_allocation_sensitivity: needs alloc_sensitivity.py, alloc_capscan347.py, make_configs.py, collect.py with their receipts (jobs 6048205, 6048208, 6064914-6064926).
- PIPELINE: prediction_parser_audit: needs the M07 Rust source of m07_parser_audit (only the binary was located in the rerun folder) and collect.py.

## v30 reference (added 2026-10-06)

- UNTRACED: the v28d quote of NUMBERS.tsv rows 4 (48), 5 (15), 6 (32), 7 (one), 32 (710), 33 (112,835,982), 40 (637), 41 (58), 43 (0.12%), 44 (0.60%), 45 (1.96%), 47 (17,935), 59 (fourteen), 60 (1 × 10⁻⁵) was not located in the v30 text (`visible_v30f.txt`); the wording or the value changed or the sentence was removed. v30_status NOT_IN_V30; the v28d columns are unchanged.

## v33 sync (2026-10-07)

- UNTRACED: v33 Supplementary Methods (paragraph about assembly-local names) prints the 295-acquisition header audit (2,393,928 of 29,163,863 sequences; 2,494,126 of 2,494,236 shared names). The table is `source_data/supplementary_methods_1_5/header_audit.tsv` (`pipelines/metagenome_header_audit`), but no NUMBERS row covers these values.
