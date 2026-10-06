# Supplementary Table 3

PREDICT antibiotic endpoints under alternative handling of rifaximin recorded for hepatic
encephalopathy. The primary exposure rule does not use the field that records rifaximin
for hepatic encephalopathy. This table repeats the paired on-minus-off analysis of
Fig. 4e-g with that field handled three other ways, so that a reader can see that the
two endpoints (human-assigned share of total LFQ intensity and quantified protein groups)
keep their direction and size under each definition.

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| Endpoints sheet | per definition and endpoint: paired patients and acquisitions, state means, mean and median on-minus-off difference, patients higher on, 95% CI, sign-flip P, BH q | 137 (definition B), 479, 793 | 42, 40, 35, 38 patients (definitions A-D) from 306 acquisitions | mean paired difference; percentile patient bootstrap; two-sided Monte Carlo sign-flip test; Benjamini-Hochberg over two endpoints per definition | `Table` (rebuild), `Statistics (controlled inputs)` (recompute) | `endpoints_by_definition.tsv`, `endpoints_producer_output.tsv` | predict_rifaximin_sensitivity (ported in `panels.py`), predict_clinical_refit |
| Composition sheet | mean origin composition of total LFQ intensity by state, per definition | 479 | as above | descriptive | `Table`, `Statistics (controlled inputs)` | `composition_by_definition.tsv` | as above |
| Read_me sheet | sheet descriptions | 479 | NA | NA | `Table` | text constant `panels.README_LINES` | NA |
| Methods counts | acquisitions reclassified, patients affected, paired n per definition | 793 | 306 acquisitions | counts | `Statistics (controlled inputs)` | `definition_counts.tsv` | as above |

Definitions (column `scenario`):

| scenario | rule |
|---|---|
| A_current | primary rule: On if any core or treatment antibiotic flag is positive, Off if all nine core flags are negative, else Unknown; the hepatic-encephalopathy rifaximin field is not used |
| B_HE_as_exposure | that field added to the core flags |
| C_HE_patients_excluded | patients with any acquisition carrying that field removed |
| D_HE_off_rows_excluded | only acquisitions that are Off under A and carry that field set to Unknown |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `SEED` | 20260907 | bootstrap seed of the primary analysis (Methods, v28d line 795); sign flips use `SEED + endpoint index` | none |
| `N_BOOTSTRAP` | 50,000 | as primary analysis (v28d line 795) | none |
| `N_SIGNFLIP_BATCHES` x `SIGNFLIP_BATCH` | 100 x 2,000 = 200,000 | as primary analysis (v28d line 795) | an independent implementation with 10^6 draws and another seed gave P within Monte Carlo error (PROVENANCE.md) |
| `CI_QUANTILES` | 0.025, 0.975 | 95% percentile interval | none |
| `TIE_TOLERANCE` | 1e-12 | null means within 1e-12 of the observed absolute mean count as extreme | none |
| `PSEUDOCOUNT` | 1 | plus-one correction (v28d line 795) | none |
| `ENDPOINTS` | human-assigned share x 100; protein-group count | the two primary endpoints; one BH family per definition (v28d line 893) | none |
| `CORE_FLAGS` | nine core antibiotic fields | primary exposure rule (v28d line 793) | definitions B-D in this table |
| `HE_FIELD` | rifaximin for hepatic encephalopathy | the field the sensitivity varies | definitions B-D |

## Known limitations

- The statistics need patient-level clinical inputs (antibiotic fields, patient identifiers). They are controlled (`controlled/CONTROLLED_INPUTS.md`); without them the notebook rebuilds and checks the published workbook only.
- P values are Monte Carlo estimates from 200,000 sign assignments; their resolution is about 5 x 10^-6.
- The table is a sensitivity analysis of the exposure definition only; acquisition aggregation, endpoints and resampling are those of the primary analysis.
