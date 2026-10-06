# source_data/STab3: data dictionary

Supplementary Table 3, PREDICT antibiotic endpoints under alternative handling of
rifaximin for hepatic encephalopathy (HE). All tables are aggregates over patients; none
contains patient or acquisition identifiers or clinical fields. How each file was derived:
`tables/STab3/prepare_source_data.py` and `tables/STab3/PROVENANCE.md`.

Definitions (`scenario`): `A_current` primary rule (HE field not used);
`B_HE_as_exposure` HE field counted as exposure; `C_HE_patients_excluded` patients with
any HE record removed; `D_HE_off_rows_excluded` Off acquisitions with an HE record set to
Unknown.

## endpoints_by_definition.tsv

The Endpoints sheet of the published workbook (Source Data package
`Supplementary_Tables/Supplementary_Table_3.xlsx`, sheet `Endpoints`), one row per
definition and endpoint (8 rows).

| column | meaning | unit |
|---|---|---|
| scenario | definition id | - |
| definition | definition in words | - |
| metric | endpoint | - |
| n_patients | patients with acquisitions in both known states | patients |
| mean_off, mean_on | mean over patients of the patient-level value in each state | percentage points or protein groups |
| mean_on_minus_off | mean of patient-level on-minus-off differences | as above |
| median_on_minus_off | median of the same differences | as above |
| n_patients_on_greater_than_off | patients with a positive difference | patients |
| CI95_low, CI95_high | 95% percentile bootstrap interval of the mean difference (50,000 patient resamples) | as above |
| p_signflip_MC | two-sided Monte Carlo sign-flip P (200,000 draws, plus-one) | - |
| q_BH_two_endpoints_within_definition | Benjamini-Hochberg q over the two endpoints of this definition | - |
| n_acquisitions_off, n_acquisitions_on | acquisitions of the paired patients in each state | acquisitions |
| units | unit of the endpoint | - |

Patient-level value: human-assigned share is the mean over the patient's acquisitions in
that state of the human-assigned fraction of total LFQ intensity, times 100; protein-group
depth is the median of quantified protein groups over those acquisitions.

## endpoints_producer_output.tsv

The same statistics as written by the analysis script, at full float precision, with the
producer's column names: `n_patients_increase` (= n_patients_on_greater_than_off),
`q_two_clinical_endpoints` (= q_BH_...), `acq_off`, `acq_on` (= n_acquisitions_...);
`metric` is `Human` or `pg_parsimony`. The recomputation in the notebook is compared with
this file bit for bit.

## composition_by_definition.tsv

The Composition sheet: mean origin composition of total LFQ intensity (fractions, 0 to 1)
of the paired patients, by definition and state (8 rows). `state` is `Off` or `On` (the
sheet calls this column `antibiotic_state`). Columns `Human`, `Microbial_species`,
`Microbial_above_species`, `Other_annotated_origins`, `Unresolved_origin` are the five
origin classes (Supplementary Methods 5).

## primary_paired_statistics.tsv

The primary (definition A) statistics of Fig. 4e-g as produced by the primary analysis
(2 rows). Definition A of this table must reproduce it (relative tolerance 1e-10). Columns
as in `endpoints_producer_output.tsv`; `units` names the endpoint unit.

## definition_counts.tsv

Aggregate counts, `count` (key) and `value`.

| key | meaning |
|---|---|
| acquisitions | acquisitions in the PREDICT roster (306) |
| HE_flag_values.<v> | acquisitions with HE field value `<v>` (`0`, `1`, empty = missing) |
| HE_by_current_state.<v>\|<state> | acquisitions by HE field value and state under definition A |
| HE_rows_patients | patients with at least one acquisition with HE field 1 |
| HE_off_rows | acquisitions Off under A with HE field 1 |
| HE_off_row_patients | patients contributing those acquisitions |
| state_changes_A_to_B.<a>-><b> | acquisitions whose state changes from A to B |
| HE_with_core_rifaximin | HE acquisitions that also have the core rifaximin field 1 |
| HE_with_treatment_rifaximin | HE acquisitions that also have the treatment rifaximin field 1 |
| paired_<scenario>.n | paired patients under the definition |
| paired_<scenario>.dropped_vs_A, .added_vs_A | paired patients lost or gained relative to A |
| paired_<scenario>.HE_patients_in_pairs | paired patients with any HE record |
| HE_patients_in_A_pairs | paired patients under A with any HE record |
| A_pairs_with_HE_off_rows | paired patients under A with an Off acquisition carrying an HE record |
| gate_A_reproduces_clinical_refit | definition A reproduced the primary statistics |
