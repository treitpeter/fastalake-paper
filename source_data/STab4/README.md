# Source data: Supplementary Table 4

Tables read by `tables/STab4/STab4.ipynb` in addition to the Supplementary Figure 10
tables it shares (`source_data/SFig10/`: patients.tsv, focused_values.tsv,
focused_features.tsv, legacy_85_values.tsv, legacy_features.tsv, legacy_259_values.tsv,
augmented_pair_screen.tsv, all_EDS_associations.tsv, selected_pair_summary.tsv,
separate_adjustments_all_models.tsv; see `source_data/SFig10/README.md`).

They were derived from the published workbook and the Source Data package by
`tables/STab4/prepare_source_data.py` (origins and hashes in `tables/STab4/PROVENANCE.md`).
Patients are keyed by the pseudonym of the Source Data package (PREDICT_Pnnn).
Acquisition and microbiome sample identifiers, clinical group and antibiotic state are
not here (controlled input, `controlled/CONTROLLED_INPUTS.md`).

Read floats with `pandas.read_csv(..., sep="\t", float_precision="round_trip")`; the
values then equal the published cells bit for bit.

## observations_259.tsv
sha256: dfb2eadb66f4adbe4fe08d06b1870452b0f34d0ea439ed4b11171f397b53a83c
rows: 259 (one row per matched proteome-metabolome observation, repeated visits included, in the column order of sheet All20_features_259_observations)
produced by: predict_eds_screen

| column | unit | meaning |
|---|---|---|
| observation_index | id | 1-259, column order of sheet All20_features_259_observations; same index as observation_index in source_data/SFig10/legacy_259_values.tsv |
| patient_id | id | pseudonymised patient identifier (PREDICT_Pnnn); empty for 3 observations that carry no patient identifier in the producer table either |
| Human | fraction | human-assigned share of total LFQ intensity of the acquisition |
| pg_parsimony | count of protein groups | quantified protein groups of the acquisition (protein-group depth) |
| dysbiosis_score | dimensionless | metagenomic Enterotype Dysbiosis Score (0-1) |

## patients_85_observation_index.tsv
sha256: 2bdc49153ce168703b8e1e316472158f1fbfae336c235d4982a040176ee44fd8
rows: 85 (one row per patient, EDS order)
produced by: predict_eds_screen

| column | unit | meaning |
|---|---|---|
| patient_position | id | display position 1-85 (as in source_data/SFig10/patients.tsv) |
| observation_index | id | the observation (1-259) used for the patient: earliest matched proteome-metabolome observation (published column `index` + 1) |

## read_me.tsv
sha256: 615a22f724377b0d20667e91f847dc78aab6e5b3a2ad9161ee8a8a32fe9ba5be
rows: 9 (one row per line of the published Read_me sheet)
produced by: Source Data package builder (specials.py ST3_README, renumbered to Supplementary Table 4)

| column | unit | meaning |
|---|---|---|
| line | id | line number in the sheet |
| text | id | text of the line; empty for the blank line |

## five_pair_all_models.tsv
sha256: 1f6400e71348f3d87d5248f407f7b9f8f650f4334e9a5eddf21938d315ea1028
rows: 320 (5 pairs x 2 subsets x 2 normalizations x 16 covariate combinations)
produced by: predict_eds_pair_adjustments (copied unchanged from the Source Data package, Biology_refit/eds_refit/separate_adjustments/all_models.tsv)

| column | unit | meaning |
|---|---|---|
| id | id | pair name |
| metabolite | id | metabolite feature id |
| protein | id | KO or KEGG map id |
| layer | id | KO or Pathway |
| scope | id | pair_complete or all_pairs_complete (the 84 patients shared by all five pairs) |
| normalization | id | native or annotated_relative |
| model | id | covariates joined by '+', or none |
| n | count of patients | subset size |
| rho | dimensionless | partial rank correlation |
| delta_from_unadjusted | dimensionless | rho minus the unadjusted rho of the same subset and normalization |
| covariate_df | count | covariate columns (rank of the design minus the intercept) |
| nominal_p | dimensionless | nominal two-sided P of the partial correlation (not a discovery test) |
| lower95, upper95 | dimensionless | pointwise 2.5th and 97.5th bootstrap percentiles of rho (pair_complete only) |
| delta_lower95, delta_upper95 | dimensionless | same for the change from the unadjusted rho |
| bootstrap_draws | count | bootstrap draws (1000; empty for all_pairs_complete) |
