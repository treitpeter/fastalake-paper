# Supplementary Table 4

Enterotype Dysbiosis Score (EDS) displays in the PREDICT cohort. The workbook keeps the
earlier 20-feature, 259-observation selection and its 85-patient display (Fig. 4h,
Supplementary Fig. 10b), and gives, for the three post hoc metabolite-enzyme pairs of
Supplementary Fig. 10a, the six displayed measurements, the feature associations, the q
values saved from the complete screen families and all 192 descriptive adjustment
estimates. `STab4.ipynb` rebuilds all eleven sheets from `source_data/`, compares them
cell for cell with the published workbook, and recomputes the statistics the table
reports. The table is unchanged between v28d and v28f.

## Sheets

| sheet | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| Read_me | sheet guide, as published | 481 | NA | none | Table | source_data/STab4/read_me.tsv | Source Data package builder (PROVENANCE.md) |
| All20_features_259_observations | 20 features x 259 observations, raw values | 481, 801 | 259 observations (repeated visits) | none | Table | source_data/SFig10/legacy_259_values.tsv, legacy_features.tsv | predict_eds_screen; selection predict_eds_feature_selection |
| Fixed_feature_key | feature block, label, positive patients, rho with EDS | 481, 801 | 85 patients | Spearman (descriptive) | Table; Statistics: overview features and z scores | source_data/SFig10/legacy_features.tsv | predict_eds_feature_selection, predict_eds_screen |
| Observation_metadata | per observation: pseudonym, human-assigned share, protein-group depth, EDS | 481 | 259 observations | none | Table | source_data/STab4/observations_259.tsv | predict_eds_screen |
| 85_patient_order | per patient in EDS order: observation used, pseudonym, share, depth, EDS | 481, 799 | 85 patients | none | Table | source_data/SFig10/patients.tsv, source_data/STab4/patients_85_observation_index.tsv | predict_eds_screen |
| 85_patient_raw, 85_patient_z | 20 features x 85 patients, raw and within-feature z | 481, 801 | 85 patients | none | Table; Statistics: overview features and z scores | source_data/SFig10/legacy_85_values.tsv | predict_eds_screen |
| Focused_raw | six displayed measurements per patient | 481, 801 | 85 patients | none | Table | source_data/SFig10/focused_values.tsv and the two patient tables | predict_eds_focused_pairs |
| Focused_features | association of each displayed feature with EDS, full-family q | 481 | 85 or 84 patients | Spearman, BH | Statistics: the six displayed features | source_data/SFig10/focused_features.tsv, all_EDS_associations.tsv | predict_eds_screen |
| Focused_pairs | pair n, rho, saved full-family q, joint-adjusted rho and q | 137, 481 | 85, 84, 85 patients | Spearman and partial rank correlation, BH | Statistics: the three pairs | source_data/SFig10/selected_pair_summary.tsv, augmented_pair_screen.tsv | predict_eds_screen, predict_eds_focused_pairs |
| Focused_adjustments | 192 partial rank correlations | 137, 481, 799 | 85, 84, 85 and common 84 patients | point estimates only | Statistics: 192 adjustment estimates | source_data/SFig10/separate_adjustments_all_models.tsv | predict_eds_focused_pairs |
| (not a sheet) five-pair decomposition | counts behind v28d line 799 (five pairs, 16 combinations, 84 patients, 1,000 draws) | 799 | patients | percentile bootstrap (not rerun) | Statistics: five-pair decomposition | source_data/STab4/five_pair_all_models.tsv | predict_eds_pair_adjustments |

The repository workbook (`build/STab4/Supplementary_Table_4.xlsx`) withholds acquisition
and microbiome sample identifiers, clinical group and antibiotic state; sample
identifiers are replaced by observation index (1-259) or patient display position
(1-85, EDS order). With `FASTALAKE_CONTROLLED_DIR` set (see
`controlled/CONTROLLED_INPUTS.md`), the notebook also rebuilds the complete published
layout and compares it with the published workbook itself.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| N_PATIENTS | 85 | earliest matched proteome-metabolome observation per identified patient (v28d line 799) | all 259 matched observations: sheet All20_features_259_observations |
| MIN_POSITIVE_N | 43 | minimum patients with positive finite values for a test (v28d line 799) | none reported |
| Q_CUTOFF | 0.05 | candidate and significance threshold (v28d line 799) | none reported |
| PAIRS | three metabolite-KO pairs (HILIC_FT23003-K01812, HILIC_FT24070-K01667, C18_F2423-K01442) | selected post hoc for annotated reaction links and distinct topics (v28d line 799) | 21 candidate pairs: figures/SFig10 |
| COVARIATES | human-assigned share, log2 protein-group depth, antibiotic state, clinical group | sensitivity covariates (v28d line 799) | all 16 combinations: Focused_adjustments |
| DEPTH_TRANSFORM | log2, then ranked | v28d line 799 | none reported |
| NORMALIZATIONS | native, annotated_relative | v28d line 799 | both in Focused_adjustments |
| SCOPES | pair_complete, common_84 | v28d line 799 | both in Focused_adjustments |
| JOINT_MODEL | host_depth_antibiotic_clinical | joint adjustment (v28d lines 137, 799) | EDS-adjusted model: augmented_pair_screen.tsv |
| FIVE_PAIR_BOOTSTRAP_DRAWS | 1000 | v28d line 799 | none reported |
| FIVE_PAIR_SEED_BASE | 2026093032 (+ pair index) | producer run_analysis.py; not stated in v28d | none |
| FIVE_PAIR_INTERVAL | 2.5th and 97.5th percentiles | pointwise 95% percentile intervals (v28d line 799) | none |

## Known limitations

- The repository workbook is not byte-identical to the published one: it withholds four
  columns and replaces sample identifiers. All other cells are equal (checked cell for
  cell in `checks/workbook_vs_published_redacted.tsv`; with the controlled input every
  cell of the published workbook, `checks/workbook_vs_published_full.tsv`).
- 144 of the 192 adjustment estimates and the three joint-adjusted correlations need
  antibiotic state and clinical group, which are controlled inputs.
- Saved full-family q values are recomputed from the saved screen P values; the screen
  itself (predict_eds_screen) is not rerun here.
- The five-pair bootstrap of v28d line 799 is not part of the workbook and is not rerun.
- The Read_me sheet cross-references figure panels differently from the v28d legends
  (OPEN_ITEMS.md).
- No schematics: there is no ARTWORK_NUMBERS.md.
