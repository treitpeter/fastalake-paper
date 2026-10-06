# Source data: Supplementary Figure 4

Tables read by `figures/SFig04/SFig04.ipynb`. All come from the minimum-length-7
sensitivity analysis (`minlen7_search_sensitivity`, see `figures/SFig04/UPSTREAM.json`). How
each file was derived from the original outputs is in `figures/SFig04/prepare_source_data.py`
and `figures/SFig04/PROVENANCE.md`. No clinical covariates are included. Acquisitions are
identified by their raw-file names.

Length populations used throughout:

| code | meaning |
|---|---|
| `baseline5` | production search, 5-50 residues (published results) |
| `nested7` | production results filtered to 7-50 residues, original q values |
| `fresh7` | new search with minimum length 7, confidence re-estimated |
| `nested9` | production results filtered to 9-50 residues, original q values |
| `fresh7_nested9` | new 7-50 search filtered to 9-50 residues |

Database arms: `FL` FastaLake (CAMPI, healthy stool); `FL_public` public-sequence FastaLake;
`FL_cohortMG` assembly-inclusive FastaLake; `FL_cohortMG_lfq` the five repaired PREDICT
FL_cohortMG searches rerun with LFQ on; `MG` specimen-matched metagenome; `REF` CAMPI
reference proteomes; `FL_R` FastaLake plus reference; `FL_R_M` FastaLake plus reference
plus metagenome.

## benchmark_matched_per_acquisition.tsv
sha256: 0d0e4cb3340b2ed8db6943a0e97e834a663f1872a6c3ab1695b3d75a835b6753
rows: 2890 (one row per acquisition and length population; 296 PREDICT and 282 CapScan matched pairs, 5 populations)
produced by: minlen7_search_sensitivity

| column | unit | meaning |
|---|---|---|
| cohort | id | PREDICT or CAPSCAN |
| arm | id | FastaLake arm, always FL_public |
| comparator | id | comparator database, always MG |
| sample | id | acquisition (raw-file name) |
| population | id | length population (table above) |
| arm_peptides | count of peptides | accepted distinct peptides, FastaLake search |
| comparator_peptides | count of peptides | accepted distinct peptides, MG search of the same spectra |
| shared | count of peptides | peptides accepted in both searches |
| gain_pct | % | 100 (arm_peptides / comparator_peptides - 1) |
| comparator_identity_recovery_pct | % | 100 shared / comparator_peptides |

## per_acquisition_length_counts.tsv
sha256: 4ba1f5dcd266997502f4e96de7dfc2897e524310df3b43baf37067cb650c058e
rows: 1990 (one row per repeated search: acquisition and database arm, v28d roster with 331 CapScan)
produced by: minlen7_search_sensitivity

| column | unit | meaning |
|---|---|---|
| cohort | id | CAMPI, HSTOOL, PREDICT or CAPSCAN |
| arm | id | database arm |
| sample | id | acquisition |
| baseline5, fresh7, fresh7_nested9, nested7, nested9 | count of peptides | accepted distinct peptides in each length population |
| gained_fresh_nested7 | count of peptides | in fresh7 but not nested7 |
| lost_fresh_nested7 | count of peptides | in nested7 but not fresh7 |
| fresh7_vs_base5 | % | 100 (fresh7 / baseline5 - 1) |
| fresh7_vs_nested7 | % | 100 (fresh7 / nested7 - 1) |
| nested7_vs_base5 | % | 100 (nested7 / baseline5 - 1) |

## per_acquisition_length_counts_participant15.tsv
sha256: 53c749ed5a550d77cae7f97a2d4b183708dea2f31a0ab228ac5aac826e74223d
rows: 16 (one row per search; the 16 CapScan participant-15 FL_public acquisitions added for the 347 roster)
produced by: minlen7_search_sensitivity

Same columns as `per_acquisition_length_counts.tsv`. Panel c appends these rows, giving 347
CapScan FL_public acquisitions.

## campi_member_counts.tsv
sha256: 7033b5a9e91b29446e095a0227bfb89dee83b460fb6bef3480c7ffee63ce0f49
rows: 200 (one row per arm, length population and CAMPI community member)
produced by: minlen7_search_sensitivity

| column | unit | meaning |
|---|---|---|
| arm | id | FL, FL_R, FL_R_M, MG or REF |
| population | id | length population |
| member | id | CAMPI community member |
| pooled_specific | count of peptides | distinct member-specific peptides pooled over the 15 acquisitions |
| acquisitions_ge1 | count of acquisitions | acquisitions with at least one member-specific peptide |
| median_per_acq | count of peptides | median member-specific peptides per acquisition |
| per_acquisition | count of peptides | per-acquisition counts, semicolon-separated, acquisition order S01-S15 |

## campi_fusion_pooled_q01.tsv
sha256: 669a785968dca936fe9ed9acdc45bf205558ca711c845bf77459374e4e3bc073
rows: 10 (one row per control seed, search and length floor, at peptide q <= 0.01)
produced by: minlen7_search_sensitivity

| column | unit | meaning |
|---|---|---|
| seed | id | fusion-control construction seed |
| search | id | baseline (production 5-50 search) or fresh7 (new minimum-7 search) |
| min_residues | residues | length floor applied to the accepted peptides |
| acq | count of acquisitions | acquisitions pooled |
| n | count of peptides | sum of accepted canonical peptides |
| num | count of peptides | sum of paired numerators |
| unp | count of peptides | sum of unpaired accepted peptides (worst case) |
| pooled_diag_pct | % | 100 (num + unp) / n, drawn in panel e |
| paired_only_pct | % | 100 num / n |

## campi_fusion_per_acquisition_q01.tsv
sha256: 749eb922b7f68756c44aaed57905904e2fd926de573577873c0e2150b11e3802
rows: 150 (one row per seed, acquisition, search and length floor, at q <= 0.01)
produced by: minlen7_search_sensitivity

Per-acquisition rows whose sums give `campi_fusion_pooled_q01.tsv`
(`panels.stat_e_pooled_diagnostic`). Columns used: `seed`, `sample` (id), `search` (id),
`q_threshold` (dimensionless), `min_residues` (residues), `n_canonical`,
`paired_numerator`, `n_unpaired_accepted` (count of peptides). The other columns are
kept as written by the fusion-control run: `n_target_only`, `n_entrapment_only`,
`n_shared`, `n_paired_entrapment` (count of peptides);
`paired_upper_diagnostic_with_unpaired_worst_case`,
`combined_upper_diagnostic_with_unpaired_worst_case`, `paired_subset_fraction`
(fraction); `whole_procedure_fdr_control_established` (boolean, always False).

## healthy_strategies_per_acquisition.tsv
sha256: 8f8e7b225b9c9652c382e2a3bf875079a230d0f685f06c281be4f3f069aa598a
rows: 232 (one row per acquisition, length population and strategy; 58 acquisitions)
produced by: minlen7_search_sensitivity

| column | unit | meaning |
|---|---|---|
| sample | id | healthy-stool acquisition |
| population | id | baseline5 or fresh7 |
| strategy | id | common (common FASTA) or joint (joint search with study-level confidence) |
| own | count of peptides | accepted peptides with the acquisition's own database |
| test | count of peptides | accepted peptides with the strategy |
| count_change_pct | % | 100 (test / own - 1) |
| own_identity_recovery_pct | % | percentage of own-database peptides also accepted with the strategy |

## grouping_study_groups.tsv
sha256: 8b9f580b4b783f8987c5887065b345b705af0a9817ef734a84cff9c73094a9bc
rows: 8 (one row per study-level grouping task)
produced by: minlen7_search_sensitivity

| column | unit | meaning |
|---|---|---|
| task | id | cohort and database of the grouping |
| roster | id | M09 (original run) or 347 (CapScan regrouped on the 347-acquisition roster) |
| arms_used | id | database arm and number of acquisitions |
| acquisitions_baseline, acquisitions_fresh7 | count of acquisitions | acquisitions grouped |
| baseline | count of protein groups | inferred study groups, production 5-50 searches |
| fresh7 | count of protein groups | inferred study groups, new minimum-7 searches |

## search_inventory.tsv
sha256: 5b8811d1f533088cf0f59875da8373abaf8fce113000a2d51629e5ebb91e2bd9
rows: 2022 (one row per repeated search)
produced by: minlen7_search_sensitivity

| column | unit | meaning |
|---|---|---|
| run | id | M09_331 (1,990 original searches), A2_347 (16 CapScan FL_cohortMG participant-15 searches), public347 (16 CapScan FL_public participant-15 searches) |
| task_index | id | task index within the run |
| cohort | id | cohort |
| arm | id | database arm |
| sample | id | acquisition |
| lfq | boolean | label-free quantification enabled |
| sage_version | id | search engine version |
| sage_exit | dimensionless | exit code |
| receipt_pass | boolean | the run's receipt checks passed (settings equal to production except minimum length, no accepted peptide shorter than 7) |
