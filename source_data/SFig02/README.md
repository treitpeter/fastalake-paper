# Source data: Supplementary Figure 2

Origins, package comparison and transformations: `figures/SFig02/PROVENANCE.md`.

## selection_length_lookup.tsv
sha256: 63594d955d35bd564a32fe68317274537519376858d3bc0daee1b09cc3b1cc29
rows: 252 (one row per arm, seed and query length 9-50; six arms: authentic and five scrambles)
produced by: database_selection_specificity

| column | unit | meaning |
|---|---|---|
| arm | id | `REAL` (authentic predictions) or `SIM_<seed>` (composition-preserving scramble) |
| seed | id | scramble seed; empty for `REAL` |
| length | residues | query length |
| distinct_queries | count of sequences | distinct I/L-normalised queries of this length (denominator) |
| matching_queries | count of sequences | distinct queries contained in at least one reservoir protein |
| fraction_matched | dimensionless | matching_queries / distinct_queries; empty when the denominator is 0 |
| retained_rows | count of predictions | prediction rows behind the distinct queries |
| matching_rows | count of predictions | prediction rows whose query matched |
| row_weighted_fraction | dimensionless | matching_rows / retained_rows |
| within_protein_windows | count | number of length-L windows within reservoir proteins |
| uniform20_reference | dimensionless | expected match probability for a uniform 20-letter alphabet |
| uniform20_il_mean | dimensionless | mean I/L-aware uniform-alphabet match probability over the arm's queries |
| measured_composition_mean | dimensionless | mean match probability under the measured residue composition |

## lookup_query_outcomes.tsv
sha256: 64182cf62ebb7d9992a0ffa2541d51590ead3a0f780dca2717f1e157e9100eba
rows: 207,503 (one row per distinct query within arm and seed)
produced by: database_selection_specificity

| column | unit | meaning |
|---|---|---|
| arm | id | as above |
| seed | id | as above |
| query | id | I/L-normalised peptide sequence |
| length | residues | query length |
| retained_rows | count of predictions | prediction rows with this query |
| matched | boolean | query contained in at least one reservoir protein |
| match_occurrences | count | occurrences of the query in the reservoir |
| uniform20_il_probability | dimensionless | I/L-aware uniform-alphabet match probability |
| measured_composition_probability | dimensionless | match probability under the measured composition |

## selection_grid_acquisitions.tsv
sha256: 8b472997023e29d4e24b5a9f7a06d9aceda6f45fc4c45391825c075d85aa3425
rows: 1,102 (one row per arm and acquisition; 19 arms x 58 healthy-stool acquisitions)
produced by: healthy_stool_selection_sweep

| column | unit | meaning |
|---|---|---|
| arm | id | `ml<L>` minimum length L at top 30%, `tp<f>` top fraction f at length nine, `ml9_pool97` older 97-pool reference |
| kind | id | `length`, `score` or `reference` |
| min_length | residues | minimum lookup length |
| top_percent | dimensionless | top-scoring fraction retained (0-1) |
| sample | id | acquisition (raw file name) |
| s3_seeds | count of sequences | lookup queries of the acquisition (`total_peptides` of the per-acquisition lookup statistics) |
| s3_match_rate | dimensionless | `match_rate` of the per-acquisition lookup statistics |
| s3_proteins | count of proteins | candidate proteins matched by the acquisition's queries |
| razor_selected | count of proteins | proteins in the selected FASTA |
| razor_unique_peptides | count of peptides | `unique_peptides` of the selection-parsimony statistics |
| razor_mapped_peptides | count of peptides | `total_peptides_mapped` of the selection-parsimony statistics |
| psms_q01 | count of PSMs | target PSMs at spectrum q <= 0.01 |
| decoy_frac_q01 | dimensionless | decoy PSMs / (target + decoy PSMs) at spectrum q <= 0.01 |
| peptides_q01 | count of peptides | distinct target peptides at peptide q <= 0.01 (modifications stripped, I/L distinct) |
| proteins_q01 | count of proteins | distinct target protein accessions at protein q <= 0.01 |
| scored_candidates_median | count | median Sage `scored_candidates` over target PSMs at spectrum q <= 0.01 |
| status | id | run status |

## selection_grid_summary.tsv
sha256: fcec37258d159150f7307c1fad691862c1be1226c75eb0692e2527af41244dfb
rows: 19 (one row per arm)
produced by: healthy_stool_selection_sweep

| column | unit | meaning |
|---|---|---|
| arm, kind, min_length, top_percent | as above | arm definition |
| n | count of acquisitions | acquisitions evaluated |
| n_pooled | count of acquisitions | acquisitions pooled upstream (58, or 97 for the reference arm) |
| peptides_q01_median, _q1, _q3 | count of peptides | median and exclusive quartiles over acquisitions |
| razor_selected_median, _q1, _q3 | count of proteins | as above for selected proteins |
| proteins_q01_median, _q1, _q3 | count of proteins | as above for accepted protein accessions |
| decoy_frac_q01_median, _q1, _q3 | dimensionless | as above for the decoy fraction |

## bags_bounded48_per_acquisition.tsv
sha256: dc48dee627b3687a0ee940cbfdf93e205ac3ab33059a2186b7740a75f32507c9
rows: 336 (one row per acquisition and lookup arm; 48 acquisitions x 7 arms)
produced by: bag_mass_lookup_search

| column | unit | meaning |
|---|---|---|
| cohort | id | CAMPI, HSTOOL (healthy stool) or CAPSCAN |
| sample | id | acquisition |
| arm | id | `exact`, `bags_le2`, `bags_le4`, `bags_le8`, `bags_all`, `mass_all`, `comb_all` |
| peptides | count of peptides | accepted canonical peptides at nominal peptide q <= 0.01 |
| spectra | count | spectrum count reported by the search receipt; not used by this figure |
| gained | count of peptides | peptides accepted in this arm but not in `exact` |
| lost | count of peptides | peptides accepted in `exact` but not in this arm |
| retained | count of peptides | peptides accepted in both |
| net | count of peptides | gained - lost |

## bags_bounded48_orderings.tsv
sha256: 702173f1d41bd92f8f3175358a0e146906668856e017992181846003439e53ec
rows: 819 (one row per acquisition and number of bag orderings)
produced by: bag_mass_lookup_search

| column | unit | meaning |
|---|---|---|
| cohort | id | as above |
| sample | id | acquisition |
| orderings | count | orderings admitted by the prediction's recorded bags, including the original (capped at 64) |
| predictions | count of predictions | lookup-eligible prediction rows with this number of orderings |
| lookup_eligible_predictions | count of predictions | all lookup-eligible prediction rows of the acquisition |

## bags_fragment_audit_summary.tsv
sha256: fadd04d3e48ad5add84a4b63bbdef9b2008de38823986507c1f635b477e92c41
rows: 150 (one row per acquisition, contrast, change and sequence relation; six acquisitions)
produced by: bag_mass_lookup_search

| column | unit | meaning |
|---|---|---|
| cohort | id | CAMPI, CAPSCAN, HSTOOL or PREDICT |
| sample | id | acquisition; the PREDICT acquisition is pseudonymised as `PREDICT_acq_01` |
| contrast | id | arm compared with its comparator, e.g. `bags_vs_exact` |
| change | id | `gained` or `lost` |
| relation | id | sequence relation of the pair, e.g. `composition_preserving_rearrangement` |
| spectrum_events | count of pairs | same-spectrum comparison pairs |
| zero_diagnostic | count of pairs | pairs without a diagnostic annotated fragment ion |
| positive_diagnostic | count of pairs | pairs with at least one diagnostic ion |
| median_diagnostic | count | median number of diagnostic ions per pair |
