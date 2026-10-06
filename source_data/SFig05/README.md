# Source data: Supplementary Figure 5

Acquisition names are the raw-file names of the public acquisitions. No clinical column
is present.

## fusion_pooled_current.tsv
sha256: 0f90dc0a3003ae3f5097832ec9c5cd55404647d065593adda72d04bc77ba121e
rows: 120 (one row per cohort, control seed, q threshold and length window)
produced by: healthy_stool_donor_entrapment (healthy stool, PREDICT, CapScan); campi_fusion_375m (CAMPI)
Drawn in panels a-c. Byte-identical to Source Data package sheet `Supplementary_Figure_05/fusion_pooled_current.tsv`.

| column | unit | meaning |
|---|---|---|
| seed | id | fusion control permutation seed |
| q_threshold | dimensionless | peptide q threshold |
| min_residues | count of residues | lower peptide length bound |
| max_residues | count of residues | upper peptide length bound |
| acquisitions | count of acquisitions | acquisitions summed (filled for CAMPI rows only) |
| n_canonical | count of peptides | N, accepted canonical sequences summed over acquisitions |
| paired_numerator | count of peptides | E + E₀ + 2B summed over acquisitions |
| n_unpaired_accepted | count of peptides | U, accepted sequences outside the audited one-to-one paired subset |
| diagnostic | fraction | (paired_numerator + n_unpaired_accepted) / n_canonical |
| paired_component | fraction | paired_numerator / n_canonical |
| unpaired_component | fraction | n_unpaired_accepted / n_canonical |
| median_acquisition_diagnostic | fraction | median over acquisitions of the per-acquisition diagnostic (CAMPI rows only) |
| cohort | id | CAMPI, HSTOOL (healthy stool), PREDICT, CAPSCAN |
| n_entrapment_q_tied | count of peptides | accepted entrapments with peptide q tied to their accepted target (filled for HSTOOL, PREDICT and CAPSCAN rows only) |

## fusion_acquisition_components.tsv
sha256: 84243602ce1bdac6cf5270e110a6be9d2e9d751ce95525ef89eec8ba0bd5961c
rows: 21,300 (one row per acquisition, seed, q threshold and length window)
produced by: healthy_stool_donor_entrapment; campi_fusion_375m
Byte-identical to package sheet `Supplementary_Figure_05/fusion_acquisition_components.tsv`. CAMPI rows carry the cohort label `CAMPI_MOCK`.

| column | unit | meaning |
|---|---|---|
| q_threshold, min_residues, max_residues | as above | stratum |
| n_canonical | count of peptides | N, accepted canonical sequences |
| n_target_only, n_entrapment_only, n_shared | count of peptides | accepted sequences found only in targets, only in controls, or in both |
| n_unpaired_accepted | count of peptides | U |
| n_paired_entrapment | count of peptides | E, accepted paired entrapments |
| n_entrapment_alone | count of peptides | E₀, accepted entrapments whose paired target was not accepted |
| n_entrapment_q_strict_better | count of peptides | accepted entrapments with strictly better q than their accepted target |
| n_entrapment_q_tied | count of peptides | as above, tied q; B = strict_better + tied |
| n_entrapment_score_strict_better, n_entrapment_score_tied | count of peptides | the same by score (sensitivity) |
| n_q_score_order_discordant | count of peptides | pairs whose q order and score order disagree |
| paired_numerator | count of peptides | E + E₀ + 2B |
| score_order_numerator_sensitivity | count of peptides | numerator with score ordering |
| paired_diagnostic | fraction | (paired_numerator + U) / N for the acquisition |
| score_order_diagnostic_sensitivity | fraction | the same with score ordering |
| paired_component_fraction, unpaired_component_fraction, q_tie_component_fraction | fraction | shares of the acquisition diagnostic |
| n_multiple_peptidoforms | count of peptides | sequences with several accepted peptidoforms |
| n_unknown_shared_target_control_sequence, n_unknown_partner_shared, n_unknown_nonbijective_mapping, n_unknown_outside_audited_digest | count of peptides | reasons a sequence was counted in U |
| cohort | id | CAMPI_MOCK, HSTOOL, PREDICT, CAPSCAN |
| sample | id | acquisition (raw-file name) |
| seed | id | fusion control permutation seed |

## fusion_decomposition_current375.tsv
sha256: 1b13ec931e8b0bdcd18e07a1c8341327662ec323688853682ae7964159e2374c
rows: 6 (one row per seed and length window, CAMPI, q≤0.01)
produced by: campi_fusion_375m (recount step)
Drawn on the continuation page; the exact input of its producer. The package sheet
`Supplementary_Figure_05/fusion_decomposition_current375.tsv` has the same values and differs only in the last printed digit of four floats (≤ 2 × 10^-15).

| column | unit | meaning |
|---|---|---|
| seed | id | fusion control permutation seed |
| min, max | count of residues | length window |
| n | count of peptides | N |
| E, E0, B | count of peptides | paired entrapments, entrapments without accepted target, entrapments with better or tied q |
| q_ties | count of peptides | tied part of B |
| paired_numerator | count of peptides | E + E0 + 2B |
| U | count of peptides | unresolved accepted sequences |
| paired_pct, unresolved_pct, total_pct | % | 100 × (E + E0 + 2B)/N, 100 × U/N, their sum |
| unresolved_share_pct | % | 100 × U / (E + E0 + 2B + U) |

## donor_entrapment_seed_means.tsv
sha256: f1b4658e64b3a0e4b854ad38fb5b38145b30f6ec0438d6891bdaa0703a43581a
rows: 630 (one row per target acquisition, length and q threshold; 63 targets)
produced by: healthy_stool_donor_entrapment
Drawn in panel d. Byte-identical to package sheet `Supplementary_Figure_05/donor_entrapment_seed_means.tsv`.

| column | unit | meaning |
|---|---|---|
| target | id | acquisition whose spectra were searched |
| donor | id | acquisition the database was selected from (equal to target for T0_self) |
| tier | id | T0_self (self-selected) or T1_within_group (donor-selected) |
| group | id | healthy-stool experimental design |
| min_length | count of residues | 5 or 9 (to 50) |
| q_threshold | dimensionless | peptide q threshold |
| n_target, n_entrapment | count of peptides | accepted targets T and entrapments E, mean of 3 seeds |
| paired_estimate | fraction | (E + E₀ + 2B)/(T + E), mean of 3 seeds |
| lower_fraction | fraction | E/(T + E), mean of 3 seeds |
| combined_estimate | fraction | 2E/(T + E), mean of 3 seeds |

## donor_recount.tsv
sha256: 5510b874ec3fa9559597edbf8a9ef430e62a5b9d72e612d941f4a72a3ead7e0e
rows: 1,890 (one row per target, seed, length and q threshold)
produced by: healthy_stool_donor_entrapment
Per-seed records from which `donor_entrapment_seed_means.tsv` is recomputed. The original file's leading pandas `index` column was removed (`cut -f2-`); all other bytes are unchanged.

| column | unit | meaning |
|---|---|---|
| target, donor, tier, group, min_length, q_threshold | as above | stratum |
| seed | id | FDRBench entrapment seed (20260910, 20260911, 20260912) |
| n_target, n_entrapment | count of peptides | T and E |
| partner_undiscovered | count of peptides | E₀ |
| entrapment_ranks_higher | count of peptides | B |
| exact_score_ties | count of peptides | entrapments with exactly the target's score |
| paired_numerator | count of peptides | E + E₀ + 2B |
| lower_fraction, combined_estimate, paired_estimate | fraction | as above, per seed |
| paired_ties_upper | fraction | paired estimate counting exact ties as B |

## simulation_mixture_outcomes.tsv
sha256: 783e4aa37060ded3d937be84495cfdb2c38d08fcfb07411b1ca963c9722c8387
rows: 72 (one row per mixture and ordinary arm)
produced by: restricted_reference_simulation
Dots of panels e,f. Byte-identical to package sheet `Supplementary_Figure_05/simulation_mixture_outcomes.tsv`.

| column | unit | meaning |
|---|---|---|
| arm | id | fixed, independent (selected from independently noised replicate), self (selected from the searched measurement) |
| discoveries | count of peptides | target canonical peptides at q≤0.01 |
| false_discoveries | count of peptides | discoveries not among the known generating sequences |
| realised_FDP | fraction | false_discoveries / max(discoveries, 1) |
| decoy_winners | count of peptides | decoys passing the threshold |
| empty_discovery_set | boolean | no discoveries |
| unit_id | id | synthetic mixture |
| regime | id | baseline or stress |

## simulation_all_outcomes.tsv
sha256: e6427d549e37a18d60d26447c962c274f08c2b41a58594a66edbadc7029067a1
rows: 288 (one row per mixture and arm, all 12 arms)
produced by: restricted_reference_simulation (collect step, `primary.tsv`)
Columns as `simulation_mixture_outcomes.tsv`; its 72 ordinary rows are identical to that file. Used to rerun the bootstrap in its original row order.

## simulation_uncertainty.tsv
sha256: f267eee37eaa0cf29e122b82a27914abdb309a8edc31a8d4fd9c397d20b0e0ec
rows: 8 (one row per regime and endpoint)
produced by: restricted_reference_simulation (collect step)
Squares and intervals of panels e,f. Byte-identical to package sheet `Supplementary_Figure_05/simulation_uncertainty.tsv`.

| column | unit | meaning |
|---|---|---|
| regime | id | baseline or stress |
| endpoint | id | fixed, independent, self, or self_minus_independent (paired difference) |
| n_units | count of mixtures | 12 |
| mean | fraction | unweighted mean over mixtures |
| bootstrap_low, bootstrap_high | fraction | 2.5% and 97.5% quantiles of 20,000 whole-mixture bootstrap means |
| t_low, t_high | fraction | t interval (not printed in the manuscript) |
| total_discoveries, total_false | count of peptides | sums over mixtures |
| pooled_FDP | fraction | total_false / total_discoveries |

## simulation_design.tsv
sha256: 04cb9549e146ddb018aee85b662a6e9de6deb223739acadc3384f874ec278aab
rows: 6 (one row per design constant)
produced by: written for this repository from the simulation code and reference FASTA (restricted_reference_simulation); origin files and lines are given in each row and in `figures/SFig05/PROVENANCE.md`

| column | unit | meaning |
|---|---|---|
| parameter | id | spectra_per_acquisition, present_spectra, generating_min_length, generating_max_length, reference_proteins |
| regime | id | baseline, stress or all |
| value | count | value in the producing code, or the number of FASTA records |
| origin_file | id | file of the simulation code (or the reference FASTA) |
| origin_line | count | line number in that file |
| origin_sha256 | id | sha256 of that file |

### simulation_design_by_mixture.tsv

Design values of the restricted known-truth simulation, read from each of the 24 validation mixtures (Supplementary Methods 6: 50% stress, 25% baseline absent; 7–30 residues; ~30,000-protein reference; map RES-47k). Produced by `pipelines/restricted_reference_simulation/sim_design.py`.

| column | unit | meaning |
|---|---|---|
| mixture | id | validation unit u01–u24 |
| regime | text | baseline or stress |
| seed | integer | mixture seed |
| generating_peptides | peptides | generating peptides in the mixture (2,048) |
| absent_from_search_reference | peptides | generating peptides from proteins omitted from the search reference |
| fraction_absent | fraction | absent / generating |
| min_length, max_length | residues | shortest and longest generating peptide |
| reference_proteins | proteins | CAMPI reference proteins before omission (29,673) |
| search_reference_proteins | proteins | proteins left in the search reference |
| omitted_proteins | proteins | proteins removed from the search reference |
| slurm_job | id | generation job |
