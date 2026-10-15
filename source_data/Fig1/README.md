# Source data: Figure 1

Tables read by `figures/Fig1/Fig1.ipynb`. All tab-separated files have one header row.
Origins, sha256 and the comparison with the Source Data package are in
`figures/Fig1/PROVENANCE.md`; `figures/Fig1/prepare_source_data.py` documents the copy.
No table contains clinical covariates. Acquisitions are named by their raw-file names.

Cohort labels: `CAMPI` (defined eight-member community, 15 acquisitions), `HSTOOL`
(healthy-stool variability series, 58), `PREDICT` (306), `CapScan` or `CAPSCAN` (347;
15 participants). Healthy stool is the "healthy stool" of the text.

## acquisition_sizes_726.tsv

One row per acquisition-specific FASTA (726 rows). Panel b, Results line 55.

| column | meaning |
|---|---|
| cohort | cohort label |
| sample | acquisition (raw-file name) |
| reservoir_entries | entries in the full study reservoir the FASTA was selected from |
| below_human_reference | `True` if `entries` < 20,431 |
| path | FASTA file name (always `parsimony.fasta`) |
| bytes, mtime_ns, sha256 | size, modification time (ns) and sha256 of that FASTA |
| entries | number of FASTA entries (plotted) |
| residues | total residues |
| unique_exact_sequences | distinct sequences |
| duplicate_accession_count | accessions occurring more than once |
| sequence_set_sha256 | sha256 of the sorted sequence set |

## human_reference.json

Reviewed human Swiss-Prot canonical proteome (archive 2026_03, isoforms excluded). Keys:
`entries` (20,431, the dashed line of panel b), `residues`, `unique_exact_sequences`,
`sha256` and `bytes` of the FASTA, `sequence_set_sha256`, header audits
(`isoform_accessions`, `nonhuman_headers`, `nonreviewed_headers`,
`duplicate_accession_count`), `path` (file name), `_origin_sha256` (record it was reduced from).

## length_by_arm.tsv

One row per arm and peptide length 9 to 50 (6 arms x 42 lengths). Worked acquisition
`20181231_QX7_PhGe_LC95_SA_WorkflowReplicates_2`. Panel c, Results line 59, panel a (24,560; 36,588).

| column | meaning |
|---|---|
| arm | `REAL` (authentic predictions) or `SIM_<seed>` (all-position scramble) |
| seed | scramble seed; empty for `REAL` |
| length | query length (residues) |
| distinct_queries | distinct I/L-normalised queries of that length (denominator) |
| matching_queries | of those, queries contained in at least one reservoir protein |
| fraction_matched | matching_queries / distinct_queries (plotted, x 100) |
| retained_rows, matching_rows, row_weighted_fraction | the same counted over prediction rows |
| within_protein_windows | number of length-L windows in the 233,497,702-entry reservoir |
| uniform20_reference | match probability for uniform 20-residue sequences without I/L equivalence |
| uniform20_il_mean | mean probability under uniform residues with I/L equivalence, over the arm's queries |
| measured_composition_mean | mean probability under the reservoir's residue composition |

## five_seed_summary.tsv

One row per length 9 to 50. `seeds_with_queries`; `scrambled_min`, `scrambled_median`,
`scrambled_max`: minimum, median and maximum of `fraction_matched` over the five seeds
(median plotted in panel c); `uniform20_reference`, `uniform20_il_median`,
`composition_median`: reference curves (not plotted in Figure 1).

## acquisition_support_summary_726.tsv

Three rows per acquisition (populations), 2,178 rows. Panel d, Results line 63.

| column | meaning |
|---|---|
| cohort, sample | as above (`CAPSCAN` upper case here) |
| population | `all_selected` (all proteins in the FASTA; plotted), `identified_accessions`, `native_protein_q01` |
| n | proteins in the population |
| one_top_prediction | proteins containing exactly one distinct top-30% production query |
| one_parsimony_prediction | proteins containing exactly one distinct eligible (parsimony) prediction |
| with_FL_corroborated_prediction | proteins containing an eligible prediction also accepted by the final search |
| with_REF_correct_top, with_REF_correct_all | CAMPI only: proteins with a prediction matching the reference-search peptide of the same spectrum |
| *_median columns | per-protein medians of prediction sequences, spectra and sharing counts |
| *_fraction columns | the counts above divided by `n` (`one_top_prediction_fraction` plotted, x 100) |

## worked_acquisition_prediction_support.json

A09 record of the worked acquisition, with absolute paths reduced to file names. Keys used:
`accepted_peptides` (29,663 distinct accepted peptides at peptide q <= 0.01) and
`accepted_not_predicted_any` (4,258 absent from all archived predictions). Also
`filtering` (prediction rows, eligible rows, retained top-30% rows, 24,560 distinct
queries, score threshold), `populations` (as in the support table), `inputs` and
`outputs` (file names, sizes, sha256).

## donor_recount.tsv

One row per donor-selected search, length range and q threshold (189 searches x 2 x 5 = 1,890 rows).

| column | meaning |
|---|---|
| index | search index |
| target | acquisition whose spectra were searched |
| donor | acquisition the transferred database was selected from (= target for `T0_self`) |
| tier | `T1_within_group` (donor from another acquisition of the same design; panel e) or `T0_self` |
| seed | entrapment seed (20260910, 20260911, 20260912) |
| group | healthy-stool experimental design |
| min_length | 5 (5-50 residues) or 9 (9-50) |
| q_threshold | peptide q threshold (0.001 to 0.05) |
| n_target (T), n_entrapment (E) | accepted target and entrapment peptides |
| partner_undiscovered (E0) | accepted entrapments whose paired target was not accepted |
| entrapment_ranks_higher (B) | entrapments scoring above their accepted target |
| exact_score_ties | entrapment-target pairs with equal scores |
| paired_numerator | E + E0 + 2B |
| lower_fraction | E / (T + E) |
| combined_estimate | 2E / (T + E) |
| paired_estimate | (E + E0 + 2B) / (T + E) (FDRBench paired estimator) |
| paired_ties_upper | (E + E0 + 2B + 2 ties) / (T + E) |

## donor_seed_means.tsv

`donor_recount.tsv` averaged over the three seeds per target, donor, tier, group,
min_length and q_threshold (630 rows; columns `n_target`, `n_entrapment`,
`paired_estimate`, `lower_fraction`, `combined_estimate` are seed means). Panel e plots
`paired_estimate` x 100 for `tier == T1_within_group`, `q_threshold == 0.01` (58 targets per length).

## scramble_end_to_end.tsv

One row per workflow, gate, q and length range (11 workflows x 18 = 198 rows). Panel a,
Results lines 53, 59, 61.

| column | meaning |
|---|---|
| index | workflow index (0 authentic, 1 to 10 scrambled) |
| arm | `authentic`, `all_position` (all residues permuted) or `terminal` (C-terminal residue fixed) |
| seed | scramble seed (0 for authentic) |
| selected_proteins | candidate proteins retrieved by the acquisition's queries |
| parsimony_proteins | proteins in the FASTA after razor parsimony |
| gate | `q_only` (peptide q only) or `pep01` (peptide q and Sage `posterior_error` <= -2, that is log10 PEP, PEP <= 0.01); rank-1 targets only |
| q | peptide q threshold |
| min_length, max_length | accepted-peptide length range |
| n_peptides | accepted distinct target peptides in the final Sage search |

## selection_grid_acquisitions.tsv

One row per selection setting and healthy-stool acquisition (19 settings x 58 = 1,102 rows).
Figure 1 uses the row `arm == tp0.30` of the worked acquisition (24,560 queries; 6,671,447
candidates).

| column | meaning |
|---|---|
| arm | setting: `tp<fraction>` (top fraction at minimum length 9), `ml<n>` (minimum length at top 30%), `ml9_pool97` |
| kind | `top_percent` or `length` |
| min_length, top_percent | the setting |
| sample | acquisition |
| s3_seeds | distinct lookup queries |
| s3_match_rate | fraction of queries with a reservoir match |
| s3_proteins | candidate proteins |
| razor_selected | proteins after razor parsimony |
| razor_unique_peptides, razor_mapped_peptides | predictions unique to / mapped to selected proteins |
| psms_q01, decoy_frac_q01, peptides_q01, proteins_q01 | search results at q <= 0.01 |
| scored_candidates_median | median candidates scored per spectrum |
| status | run status |

## stage1_stats.json

Study-level lookup of the pooled top-30% queries of the 58 healthy-stool acquisitions in
the healthy-stool reservoir. Keys used: `proteins_searched` (233,497,702 reservoir
entries), `total_peptides` (250,197 distinct pooled queries), `proteins_matched`
(18,007,876 study-supported proteins). Also `match_rate`, `normalize_il`,
`peptide_length_range`, `score_filter`, `total_seconds`.

## scramble_study_lookup_stats.json

Study-level lookup of the composition-matched scrambled queries (seed 20260901; 36,588
distinct queries) in the same 233,497,702-entry healthy-stool reservoir, as `stage1_stats.json`
for the authentic queries. Keys used: `proteins_searched`, `proteins_matched` (5,841; the
Fig. 1a table cell Scrambled / Study-supported proteins). Producer output
`FINAL/v30/fig1_direct233M/scr/matches.stats.json` with `database` reduced to the file name.
sha256: fc787f3313ecde21e191c954af2e27f56dfc9140a837760c2674e8ace94577a9

### worked_acquisition_accepted_peptides.tsv

Every accepted peptide of the worked healthy-stool acquisition, flagged by whether it occurs
among the archived AlphaNovo predictions (Results: 4,258 of 29,663, 14.4%; Source Data map
RES-17). Produced by `pipelines/cohort_prediction_support/res17_table.py` from the predictions
CSV (sha256 8826920d...) and the Sage `results.sage.tsv` (sha256 92339fc9...).

| column | unit | meaning |
|---|---|---|
| acquisition | text | raw file name of the worked acquisition |
| peptide_IL | sequence | accepted peptide: target, rank 1, peptide q ≤ 0.01, modifications removed, I written as L |
| length | residues | peptide length |
| best_peptide_q | q-value | lowest Sage peptide q among this peptide's PSMs |
| in_alphanovo_predictions | yes/no | whether the I/L-collapsed sequence occurs among the unmodified AlphaNovo predictions |
