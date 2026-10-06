# Open items: supplementary_methods_8_11

## Manuscript sentences that must change for v28f

- MISMATCH: v28d line 855 "All 331 CapScan and 306 PREDICT acquisitions were searched": 347 CapScan non-specific searches completed (331 production plus 16 of participant 15), all with identical settings; should read 347 (NUMBERS.tsv, line 855, v28f_change:347).
- MISMATCH: v28d line 849, evenness estimator. The text says evenness was "calculated for each participant and sampling region after averaging acquisitions within specimen and specimens within participant–region, then averaged the four capsule-region indices". The producer (LAB_ROOT/.../coauthor_comment_check_v25_2026-10-05/compute/C4_dominance_evenness_347, functions `per_acquisition` and `summarise` of `dominance_evenness.py`) computes Pielou J and the leading share per acquisition, then takes the median over all fluid acquisitions of a participant (regions pooled) and the median over stool acquisitions, and tests with an exact Wilcoxon signed-rank test. The printed 0.623 / 0.669 and 22.5% / 12.7% (Results, line 159) are these participant medians. Either the sentence is changed to describe the computed estimator or the estimator is recomputed as described (STATS.tsv evenness rows).
- MISMATCH: v28d lines 863 and 869 "fourteen participants" for paired CapScan comparisons is fifteen on the 347 roster (recomputed by figures/Fig6; here `terminal_fraction__length_9` rows have 15 participants).
- MISMATCH: v28d line 823 "11 of 24 ... three of four non-host protein–metabolite" is 9 of 24 and one of four on the 347 roster (recomputed by figures/SFig11).
- MISMATCH: v28d line 903 "331 CapScan manuscript acquisitions" is 347 (recomputed by figures/SFig15).

## Statistical details the text does not state (filled in STATS.tsv from the code)

- TEXT-GAP: CapScan regional and own-stool bootstrap intervals use seed 20261003 (2,000 participant resamples); the own-stool interval is not described at all in Supplementary Methods 8 (STATS.tsv own-stool rows).
- TEXT-GAP: functional-profile permutations use the same seed as native profiles, 20260907; line 815 states it for native profiles only, lines 827 and 835 do not.
- TEXT-GAP: held-out EDS prediction: label permutations are drawn with seed 20260915 (the first split seed) and the patient bootstrap uses seed 20260918; neither is stated (line 803-805).
- TEXT-GAP: antibiotic sign flips use seed 20260907 plus the endpoint index (tables/STab3/panels.py), not stated at line 795.
- TEXT-GAP: histidine and terminal-fraction tests: sign-flip and bootstrap seeds (A28: 20260914 + n and 20261014 + n; A24: sha256-derived per n) and bootstrap draw counts per endpoint are not stated; exact enumeration for n ≤ 16 is stated.
- TEXT-GAP: individual-function correlations use scipy.stats.spearmanr (two-sided, t approximation, df = n − 2); implementation not stated (line 829).
- TEXT-GAP: no interval is computed for the cross-omics, controlled-distance or individual-function correlations; the text neither states intervals nor their absence.
- TEXT-GAP: the evenness and leading-share tests (Results line 159) use scipy.stats.wilcoxon with method="exact"; statistic (W = 5 and W = 3), implementation and correction family (none) are not stated.
- TEXT-GAP: statistical software versions for the downstream analyses are not stated (Python 3.12.4, numpy 2.2.6, scipy 1.16.3, statsmodels 0.14.5, scikit-learn 1.7.2 for the biology analyses; Python 3.10.9, numpy 1.26.4, scipy 1.12.0 for the histidine and terminal-fraction inference). Line 907 lists search and annotation software only.

## Numbers not recomputed

- UNTRACED: analysis parameters printed in the text and used as constants by the producers but not exercised by repository code, so their NUMBERS.tsv rows are OPEN: 80% leader eligibility (line 783, two rows), 20% and five percentage points (785), 10% training prevalence, ridge grid 0.01-100 (only 0.01, 1, 10 and 100 were ever selected) and five nearest profiles (803; producer `analyze_score.py` `features`, `ALPHAS`, `distance_predict`), 10% native prevalence (813, 819), 70% metabolite finiteness (815, 827, 833), peptide q≤0.01 and top-30% (855), q<0.05 (889), 95% interval level (795, 843).
- UNTRACED: seeds that outputs cannot confirm: 20260907 (lines 795 and 815) and 20260910 (849; P is at its minimum 1/100,000 for any seed).
- UNTRACED: "all 296 available PREDICT assembly FASTAs" (line 791): count of assembly FASTAs pooled into the PREDICT reservoir not located in any shipped table.
- UNTRACED: "twelve ... strata: eight CapScan strata" (line 877): the Source Data strata table (`source_data/SFig14/spectra/all_strata.tsv`) holds only the four PREDICT strata; the eight CapScan strata and their selection denominators were not located.
- UNTRACED: line 909, the public-module parity check ("fasta_lake/study.py reproduced the production healthy-stool grouping exactly, with identical output tables"): the receipt of that comparison was not located. The two counts (59,428 groups under direct peptide assignment, 1,909,973 features) are traced only to the production grouping receipts.
- CONTROLLED: line 793 "nine core antibiotic fields" and the line 795 bootstrap (50,000 resamples, seed 20260907) are exercised by tables/STab3 with its controlled inputs; nothing new requested here.

## Shared requests

- shared request: tools/check_numbers.py does not parse Unicode superscript exponents ("1 × 10⁻⁵"), so the NUMBERS.tsv row for line 847 "1 × 10⁻⁵" (recomputed 1e-05) is reported as MISMATCH; figures/Fig5 has the same false failures (lines 155, 159, 847).
- PIPELINE: predict_eds_heldout_prediction (LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/eds_phenotype/code/analyze_score.py driven by .../verified50/B/nested_eds.py), capscan_metagenome_profiles (LAB_ROOT/analyses/capscan_metag_analysis_2026-09-10) and healthy_stool_study_grouping are not in pipelines/; pipeline_dir is null in UPSTREAM.json.

## v30 reference (added 2026-10-06)

- UNTRACED: the v28d quote of NUMBERS.tsv rows 12 (296), 14 (95%), 15 (50,000), 16 (42), 17 (20260907), 18 (200,000), 19 (two), 29 (eight), 38 (1,999), 39 (20260907), 40 (ten), 41 (six), 47 (77), 48 (70%), 59 (70%), 78 (fourteen), 82 (four), 83 (fourteen), 84 (99,999), 85 (20260910), 88 (331), 89 (306), 90 (5–50), 92 (20), 93 (two), 94 (9–50), 95 (0.01), 96 (30%), 102 (0.05) was not located in the v30 text (`visible_v30f.txt`); the wording or the value changed or the sentence was removed. v30_status NOT_IN_V30; the v28d columns are unchanged.
