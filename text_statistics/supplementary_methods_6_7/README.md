# Supplementary Methods 6-7

Supplementary Methods 6 (validation, error controls and spectrum inspection) and 7 (transfer, pooling and common study search spaces) describe the controls behind the database-selection claims and print about 220 numbers. This unit recomputes each printed result number in v28d lines 661-777, plus the matched rosters of line 633, from `source_data/supplementary_methods_6_7/`. The sections have no figure panels. Their figures and tables (Supplementary Figs. 1-7, Supplementary Tables 1-2) are separate units. Run with `python paper.py build supplementary_methods_6_7`.

What is counted as a printed number: every result, count, n, percentage, interval bound and comparison value. Design parameters typed in the text (seeds, q thresholds, length windows, cluster identities, software versions) are not results. They are the named constants at the top of `stats.py` and are listed under Choices below. Numbers written as words ("eight", "Thirteen", "none") are checked as words.

## Sections

| section | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| sequence controls | selected-FASTA sizes and accepted peptides of authentic and scrambled workflows; scrambled-query match fractions | "499-557 ... 489-536 ... 12,606 ... zero ... 29,663 ... 8.11-8.65% ... 0.72-1.03%" (667) | 11 workflows on 1 acquisition; distinct queries | none (randomised control) | Sequence controls and lookup specificity | scramble_end_to_end.tsv, scramble_lookup_comparison.tsv | scramble_end_to_end_workflow (not in pipelines/) |
| database size | selected FASTAs and reservoirs | "726 primary acquisitions", "20,431", "710 selected FASTAs" (673-675) | 726 acquisitions | none | Selected-database size | acquisition_sizes_726.tsv, human_reference.tsv | cohort_database_sizes |
| selection sensitivity | grid of 19 settings | "58 ... 19 settings ... 1,102" (679) | 1,102 acquisition-setting observations | none | Selection sensitivity | selection_grid_acquisitions.tsv | healthy_stool_selection_sweep |
| CAMPI recovery | acquisitions, members, contaminants, pooled reference peptides, depth draw | "116", "76,127", "38,254" (683-687) | 15 acquisitions; distinct peptides | none | CAMPI community recovery | campi_spectral_inventory.tsv, campi_reference_peptides.tsv, campi_member_summary.tsv, campi_reference_fasta_records.tsv | campi_variability, campi_reference_gap |
| clustering | configurations | "105" (695) | 15 x 7 configurations | none | Candidate clustering | clustering_per_acquisition_identifications.tsv | campi_candidate_clustering |
| fusion control | (E + E0 + 2B + U)/N and its components, two seeds | "2.59-2.79% ... 65-73% ... 0.36-0.45%" (715) | 15 acquisitions x 2 seeds; peptide | none (descriptive estimator; seed range) | CAMPI target-entrapment fusion control | fusion_campi_acquisition_components.tsv | campi_fusion_375m |
| reservoir occurrence | accepted and shuffled control sequences found in the reservoir | "1,517/1,544 (98.25%) ... 350/3,177 (11.02%)" (721) | unique sequences | none | Reservoir occurrence | label_validity_presence.tsv, label_validity_result.tsv | campi_fusion_375m |
| donor entrapment | roster, seeds, ties | "58", "No", "Three", "five" (727-729) | 63 target acquisitions x 3 seeds | none | Donor-selected entrapment | donor_recount.tsv | healthy_stool_donor_entrapment |
| simulation | mean FDP per arm and regime, bootstrap intervals, paired differences | "1.25% ... 1.91% ... 0.71 (0.35-1.12) ... 0.41" (733) | 12 + 12 synthetic mixtures (acquisition) | none; 95% percentile bootstrap intervals | Restricted-reference simulation | simulation_all_outcomes.tsv, simulation_design.tsv | restricted_reference_simulation |
| transfer and pooling | donor assignments; shuffle-only acceptance | "nine ... one", "2 of 376,412", "12 of 411,354" (739-743) | 58 targets; peptide-acquisition observations | none | Healthy-stool transfer and pooling | transfer_pool_comparisons.tsv, hstool_pooling_per_acquisition.tsv | healthy_stool_database_transfer, healthy_stool_pooling_controls |
| common search space | union sizes; features outside own FASTA | "153,267 ... 555,987 ... 452,255 ... 43.2%" (751) | protein sequences; quantified features | none | Common study search space | common_fasta_union_summary.tsv, union_contributions.tsv, representative_presence.tsv | healthy_stool_common_fasta, cohort_selected_fasta_unions, peptide_centric_downstream |
| spectrum split | roster, participants, reconverted inputs, fallbacks, zero yield | "93 ... 302 ... four", "Thirteen ... eight", "0.0473" (755-761) | 710 acquisitions; 93 + 15 participants | none | Spectrum-split selection | split_half_acquisitions_710.tsv, split_half_model_fallbacks.tsv, split_half_zero_yield.tsv | cohort_spectrum_split_selection |
| joint search | per-acquisition confidence gains, permutation control, pair CV and r, leading species | "32.0% ... 32.9% ... 17 ... none", "20.1% to 15.2% ... 70.8% to 66.5%" (765-769) | 58 acquisitions; 66 and 105 pairs | none | Joint healthy-stool search | joint_counts_per_acquisition.tsv, joint_fileq_coverage.tsv, microbial_group_pairs.tsv, technical_agreement_summary.tsv | healthy_stool_joint_search_q, healthy_stool_common_fasta |
| cross-fitting | fold sizes, validation checks, perturbation and fallback sensitivity, gains | "9.9-10.1% ... 33 of 36 ... 0.16 ... −8.4 to −2.5% ... 67% and 96%" (773-777) | 578 acquisitions (296 PREDICT, 282 CapScan) | none; 95% participant cluster bootstrap intervals | 10-fold spectrum cross-fitting | crossfit_*.tsv, confidence_ordinary_all_pairs.tsv, minlen7_benchmark_per_acquisition.tsv | capscan_crossfit_10fold, clinical_public_only_comparison, minlen7_search_sensitivity |
| confidence rosters | matched rosters and CAMPI fusion searches | "296 ... 282 ... 30 ... three" (633) | acquisitions; searches | none | Additional-confidence rosters | confidence_ordinary_all_pairs.tsv, confidence_fusion_all.tsv, confidence_pep01_zero_yield_cases.tsv | clinical_public_only_comparison, campi_fusion_375m |

No line of these sections reports a hypothesis test. STATS.tsv lists the 43 intervals and descriptive estimators with their construction, seed and implementation.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| SCRAMBLE_SEEDS | 20260901, 20260913, 20260914, 20260915, 20260916 | five realisations of each shuffle type (line 665) | all five curves in scramble_lookup_comparison.tsv |
| Q_CUTOFF | 0.01 | peptide q threshold of every endpoint in these sections | 0.001-0.05 in the fusion and confidence tables |
| FULL_LENGTH, NESTED_LENGTH, SHORT_LENGTH | 5-50, 9-50, 5-8 residues | accepted population, lookup-length subset and short peptides (lines 715, 721, 777) | all three reported |
| LOOKUP_TOP_FRACTION | 0.30 | production lookup queries (lines 677, 679) | 5-100% in the selection grid (line 679) |
| SELECTION_GRID_LENGTHS | 6-12 | minimum evidence lengths of the grid (line 679) | selection_grid_acquisitions.tsv |
| CLUSTER_IDENTITIES, CLUSTER_COVERAGE | 0.99-0.50; 0.8 | Linclust settings (line 695) | Supplementary Table 1 |
| SECOND_ENGINE | MSFragger 4.1 + Percolator 3.06.5 | second search workflow (line 699) | Supplementary Fig. 3 unit |
| AUDIT_SAGE | Sage v0.14.7, five candidates | spectrum-level audit (line 703) | Fig. 4 unit |
| FUSION_SEEDS | 20260907, 20260908 | the two fusion-control realisations (from the data; not printed) | both seeds reported separately |
| LABEL_SHUFFLE_SEED | 20260925 | composition-preserving shuffles (line 721) | none |
| DONOR_ENTRAPMENT_SEEDS | 20260910, 20260911, 20260912 | FDRBench seeds (from the data; not printed) | per-seed rows in donor_recount.tsv |
| FDRBENCH_VERSION | 0.0.4 | entrapment generator (line 727) | none |
| SIM_BOOTSTRAP_DRAWS, SIM_BOOTSTRAP_SEED, SIM_CI | 20,000; 2026091051; 2.5% and 97.5% quantiles | as in the original collect.py (line 733; seed not printed) | t intervals in simulation_uncertainty.tsv (not printed) |
| JOINT_PERMUTATION_SEED | 0 | one within-acquisition label permutation in the original extract_fileq.py (not printed) | (D+1)/T variant reported |
| CROSSFIT_FOLDS, CROSSFIT_FOLD_KEY | 10; SHA-256 of "FastaLake-D11-crossfit-20261004:<native scan>" | fold assignment (line 773) | one-fold control (line 775) |
| CROSSFIT_BOOTSTRAP_DRAWS, CROSSFIT_BOOTSTRAP_SEED, CROSSFIT_CI | 10,000; 20261004; 2.5th and 97.5th percentiles | as in the original analyze_final.py (line 777) | region stream restarted with the same seed (line 777) |
| PERTURBATION_MOVE | 0.5 percentage points | threshold for single acquisitions moved by a one-ulp perturbation (line 777) | none |
| CROSSFIT_REGIONS | duodenum, jejunum, ileum, colon, stool | order of the region bootstrap stream | none |

## Known limitations

- Several numbers are read from producer summaries rather than recomputed from lower-level data: the leading-species agreement (technical_agreement_summary.tsv), the healthy-stool union size (common_fasta_union_summary.tsv; independently matched by union_contributions.tsv), the fallback list and zero-yield record of the spectrum split, and the cross-fit validation and one-fold-control records. The searches themselves are not rerun.
- The cross-fit bootstrap reproduces the original intervals exactly only with the original cluster order. The package table carries pseudonyms, so `crossfit_per_acquisition.tsv` adds `cluster_order`, the rank of the original participant identifier (identifiers are not included).
- Spectrum-split numbers are for the 710-acquisition roster (CapScan 331), as in v28d and the v28f Supplementary Fig. 6a-c. The 726-acquisition re-run is in source_data and printed by the notebook, but not used for NUMBERS.tsv (OPEN_ITEMS.md).
- Five numbers disagree with the stored producer records or change with the v28f roster; they carry `v28f_change` status and are explained in OPEN_ITEMS.md.
- `prepare_source_data.py`, which copied the tables from the lab archive, is kept outside the repository because it names lab paths; its transformations are listed in PROVENANCE.md.
