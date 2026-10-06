# Source data: Figure 2

Tables used by `figures/Fig2/`. Two further inputs live under other units and are read from there:
`source_data/SFig03/CAMPI_membership_<sample>.tsv` (the 15 A14 peptide membership tables) and
`source_data/Fig1/acquisition_sizes_726.tsv` (CAMPI reservoir size). Their sha256 values are checked in
`figures/Fig2/panels.py` (`MEMBERSHIP_SHA256`, `SIZES_SHA256`).

No table contains clinical or personal data; CAMPI SIHUMIx is a defined bacterial community.

## alphanovo_accuracy_per_acquisition.tsv
sha256: 82bc0030c9cafaa8b68ea79d00037dc2f44e62374f1801cbb361eeb95ad9a18c
rows: 30 (one per acquisition x length window)
produced by: campi_alphanovo_accuracy

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| window | residues | 9-50 (printed) or 5-50 (context) |
| n_reference_spectra | count of spectra | spectra with a rank-1 target PSM in the REF search at peptide q<=0.01 and a peptide in the window |
| n_with_prediction | count of spectra | of these, spectra with an AlphaNovo top prediction |
| n_with_nonempty_prediction | count of spectra | with a non-empty prediction |
| n_exact_all | count of spectra | top prediction equals the reference peptide after I->L |
| n_retained_top30 | count of spectra | prediction passes the production filter (top 30% by score within the acquisition) |
| n_exact_retained | count of spectra | passes the filter and is exact |
| n_exact_shuffled | count of spectra | exact after permuting predictions across spectra (negative control) |
| frac_with_prediction | fraction | n_with_prediction / n_reference_spectra |
| frac_exact_all | fraction | n_exact_all / n_reference_spectra |
| frac_exact_among_predicted | fraction | n_exact_all / n_with_prediction |
| frac_exact_retained_of_reference | fraction | n_exact_retained / n_reference_spectra |
| frac_exact_among_retained | fraction | n_exact_retained / n_retained_top30 |
| frac_exact_shuffled | fraction | n_exact_shuffled / n_reference_spectra |
| join_offset | count of scans | scan offset used to join predictions and PSMs (0 for all) |
| join_mass_agree_20ppm | fraction | joined pairs whose precursor masses agree within 20 ppm |
| prediction_rows | count of spectra | prediction rows in the acquisition |
| production_retained_rows | count of spectra | prediction rows retained by the production filter |

## alphanovo_accuracy_pooled.tsv
sha256: 6ad34ac899e76f30e1f01f5f9105ada2b85bd6996afd6dd4f14f50bcb48397b3
rows: 2 (one per length window)
produced by: campi_alphanovo_accuracy

| column | unit | meaning |
|---|---|---|
| window | residues | 9-50 or 5-50 |
| n_acquisitions | count of acquisitions | 15 |
| n_reference_spectra | count of spectra | spectra with a rank-1 target PSM in the REF search at peptide q<=0.01 and a peptide in the window |
| n_with_prediction | count of spectra | of these, spectra with an AlphaNovo top prediction |
| n_with_nonempty_prediction | count of spectra | with a non-empty prediction |
| n_exact_all | count of spectra | top prediction equals the reference peptide after I->L |
| n_retained_top30 | count of spectra | prediction passes the production filter (top 30% by score within the acquisition) |
| n_exact_retained | count of spectra | passes the filter and is exact |
| n_exact_shuffled | count of spectra | exact after permuting predictions across spectra (negative control) |
| pooled_frac_* | fraction | pooled counterparts of the per-acquisition fractions (sums of counts over 15 acquisitions) |
| median_frac_*, min_frac_*, max_frac_* | fraction | median, minimum and maximum of the per-acquisition fraction over 15 acquisitions |

## arm_CAMPI_FL_REF_MG__direct_acquisition_summary.tsv
sha256: e685e28241ad22d097525e3d4d07d50ee410479187f10ee66a287e387807dabd
rows: 15 (one per acquisition)
produced by: campi_simplified_grouping

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| quantified_study_groups | count of protein groups | quantified groups from the direct export; must equal quantified_groups of the simplified rows |
| canonical_peptides | count of peptides | canonical peptides |
| lfq_features | count of features | LFQ features |
| intensity | summed intensity | summed LFQ intensity |

## arm_CAMPI_FL_REF_MG__reaggregated_per_acquisition.tsv
sha256: 7e83bf1ab6548b4e2f5ad156ce94e3ca910270ce89739864a4371fd7081e5adb
rows: 30 (one per acquisition x grouping method)
produced by: campi_simplified_grouping

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| method | id | grouping method; panel c uses 'simplified' |
| quantified_groups | count of protein groups | quantified protein groups in the acquisition |
| lfq_features | count of features | quantified LFQ features |
| intensity | summed intensity | summed LFQ intensity |
| median_assigned_peptides_per_group | count of peptides | median assigned peptides per group |
| mean_assigned_peptides_per_group | count of peptides | mean assigned peptides per group |
| groups_with_at_least_two_assigned_peptides | count of protein groups | groups with two or more assigned peptides |
| groups_with_at_least_two_assigned_peptides_pct | % | previous column as % of quantified_groups (panel c percentages) |
| median_compatible_assigned_coverage_pct | % | median sequence coverage by compatible assigned peptides |
| incompatible_features | count of features | features incompatible with the group assignment |
| incompatible_feature_pct | % | incompatible features as % of features |
| incompatible_intensity_pct | % | intensity share of incompatible features |

## arm_CAMPI_FL_REF__direct_acquisition_summary.tsv
sha256: 9ddf09061b1f092444d8f6eebf3ae389b85b8516e1c8041beab44f154d4081f4
rows: 15 (one per acquisition)
produced by: campi_simplified_grouping

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| quantified_study_groups | count of protein groups | quantified groups from the direct export; must equal quantified_groups of the simplified rows |
| canonical_peptides | count of peptides | canonical peptides |
| lfq_features | count of features | LFQ features |
| intensity | summed intensity | summed LFQ intensity |

## arm_CAMPI_FL_REF__reaggregated_per_acquisition.tsv
sha256: 130d4b89fa15644a2fcbf75dd30cbc4b39dc8cf76e8abb0bf66a9481ab25be39
rows: 30 (one per acquisition x grouping method)
produced by: campi_simplified_grouping

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| method | id | grouping method; panel c uses 'simplified' |
| quantified_groups | count of protein groups | quantified protein groups in the acquisition |
| lfq_features | count of features | quantified LFQ features |
| intensity | summed intensity | summed LFQ intensity |
| median_assigned_peptides_per_group | count of peptides | median assigned peptides per group |
| mean_assigned_peptides_per_group | count of peptides | mean assigned peptides per group |
| groups_with_at_least_two_assigned_peptides | count of protein groups | groups with two or more assigned peptides |
| groups_with_at_least_two_assigned_peptides_pct | % | previous column as % of quantified_groups (panel c percentages) |
| median_compatible_assigned_coverage_pct | % | median sequence coverage by compatible assigned peptides |
| incompatible_features | count of features | features incompatible with the group assignment |
| incompatible_feature_pct | % | incompatible features as % of features |
| incompatible_intensity_pct | % | intensity share of incompatible features |

## arm_CAMPI_FL_UNIPROT__direct_acquisition_summary.tsv
sha256: 9180536d0a2e4989f54c6fe8f7c024b9b60cfbb90a0c0d678f921598b0f60d12
rows: 15 (one per acquisition)
produced by: campi_simplified_grouping

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| quantified_study_groups | count of protein groups | quantified groups from the direct export; must equal quantified_groups of the simplified rows |
| canonical_peptides | count of peptides | canonical peptides |
| lfq_features | count of features | LFQ features |
| intensity | summed intensity | summed LFQ intensity |

## arm_CAMPI_FL_UNIPROT__reaggregated_per_acquisition.tsv
sha256: 22b986eb35f3086641f23c70d6151effe43b29386beadd9adede3853b9205fac
rows: 30 (one per acquisition x grouping method)
produced by: campi_simplified_grouping

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| method | id | grouping method; panel c uses 'simplified' |
| quantified_groups | count of protein groups | quantified protein groups in the acquisition |
| lfq_features | count of features | quantified LFQ features |
| intensity | summed intensity | summed LFQ intensity |
| median_assigned_peptides_per_group | count of peptides | median assigned peptides per group |
| mean_assigned_peptides_per_group | count of peptides | mean assigned peptides per group |
| groups_with_at_least_two_assigned_peptides | count of protein groups | groups with two or more assigned peptides |
| groups_with_at_least_two_assigned_peptides_pct | % | previous column as % of quantified_groups (panel c percentages) |
| median_compatible_assigned_coverage_pct | % | median sequence coverage by compatible assigned peptides |
| incompatible_features | count of features | features incompatible with the group assignment |
| incompatible_feature_pct | % | incompatible features as % of features |
| incompatible_intensity_pct | % | intensity share of incompatible features |

## arm_CAMPI_MG__direct_acquisition_summary.tsv
sha256: 02b6b172fb2123384e50518df107c815c1b310a9b0f5f2a88b0feeb7bb20aba9
rows: 15 (one per acquisition)
produced by: campi_simplified_grouping

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| quantified_study_groups | count of protein groups | quantified groups from the direct export; must equal quantified_groups of the simplified rows |
| canonical_peptides | count of peptides | canonical peptides |
| lfq_features | count of features | LFQ features |
| intensity | summed intensity | summed LFQ intensity |

## arm_CAMPI_MG__reaggregated_per_acquisition.tsv
sha256: 0bf2d7e592b02eb3baeab8b0b3eb2a38f97ac27572aa84dd12d844dcf2ab2194
rows: 30 (one per acquisition x grouping method)
produced by: campi_simplified_grouping

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| method | id | grouping method; panel c uses 'simplified' |
| quantified_groups | count of protein groups | quantified protein groups in the acquisition |
| lfq_features | count of features | quantified LFQ features |
| intensity | summed intensity | summed LFQ intensity |
| median_assigned_peptides_per_group | count of peptides | median assigned peptides per group |
| mean_assigned_peptides_per_group | count of peptides | mean assigned peptides per group |
| groups_with_at_least_two_assigned_peptides | count of protein groups | groups with two or more assigned peptides |
| groups_with_at_least_two_assigned_peptides_pct | % | previous column as % of quantified_groups (panel c percentages) |
| median_compatible_assigned_coverage_pct | % | median sequence coverage by compatible assigned peptides |
| incompatible_features | count of features | features incompatible with the group assignment |
| incompatible_feature_pct | % | incompatible features as % of features |
| incompatible_intensity_pct | % | intensity share of incompatible features |

## arm_CAMPI_REF__direct_acquisition_summary.tsv
sha256: 9b797660a6b296cf1d464e5413bdf9e51ec19daa507f5bc36f9d7c55f1ad8b32
rows: 15 (one per acquisition)
produced by: campi_simplified_grouping

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| quantified_study_groups | count of protein groups | quantified groups from the direct export; must equal quantified_groups of the simplified rows |
| canonical_peptides | count of peptides | canonical peptides |
| lfq_features | count of features | LFQ features |
| intensity | summed intensity | summed LFQ intensity |

## arm_CAMPI_REF__reaggregated_per_acquisition.tsv
sha256: f0dc2fb2ae86939beba26d53fa68ce74d746203eb4198c077195ecefbe89643c
rows: 30 (one per acquisition x grouping method)
produced by: campi_simplified_grouping

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| method | id | grouping method; panel c uses 'simplified' |
| quantified_groups | count of protein groups | quantified protein groups in the acquisition |
| lfq_features | count of features | quantified LFQ features |
| intensity | summed intensity | summed LFQ intensity |
| median_assigned_peptides_per_group | count of peptides | median assigned peptides per group |
| mean_assigned_peptides_per_group | count of peptides | mean assigned peptides per group |
| groups_with_at_least_two_assigned_peptides | count of protein groups | groups with two or more assigned peptides |
| groups_with_at_least_two_assigned_peptides_pct | % | previous column as % of quantified_groups (panel c percentages) |
| median_compatible_assigned_coverage_pct | % | median sequence coverage by compatible assigned peptides |
| incompatible_features | count of features | features incompatible with the group assignment |
| incompatible_feature_pct | % | incompatible features as % of features |
| incompatible_intensity_pct | % | intensity share of incompatible features |

## campi_FL_5_frequencies_S01.tsv
sha256: c0d526797a8a37e0dfaf8ad1b558405a58ea392434148ffef9c2134d351f1704
rows: 21467 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_FL_5_frequencies_S02.tsv
sha256: 9795c176835d765651fa943353bbd045ea0a8741a9f8e45803b667aeb43d07d7
rows: 18965 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_FL_5_frequencies_S03.tsv
sha256: 9c5a707bde552d265d46898990468467e39170705484da5a1904e01fd3e7336e
rows: 27861 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_FL_5_frequencies_S04.tsv
sha256: 14dc7f62df8d261acd35f0e63d3d62809981d0a44218528a59b9eebd15f44eaa
rows: 27197 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_FL_5_frequencies_S05.tsv
sha256: 48728ef46b545ff7c0338557a2d020a7b1a75fb93931f6ccb23680766c1860fc
rows: 32259 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_FL_5_frequencies_S06.tsv
sha256: b4a1c0b12b786eb9a6f87bad60708a5b83234b56526b419cbcfe12760d6505ed
rows: 32788 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_FL_5_frequencies_S07.tsv
sha256: de20e52928289c9bff8e16b1539366465aba5c91ebf94fd8572014df9cf4a21a
rows: 17026 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_FL_5_frequencies_S08.tsv
sha256: 9d62091f53afff537b5ff705e0a805f5f2fbc1bcfeba8d07bd26e67d0ebd6a09
rows: 20896 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_FL_5_frequencies_S09.tsv
sha256: 2858e3696ab140545a585066c423f92c254a1fa1bf5dea0fdde378344178a0cb
rows: 20483 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_FL_5_frequencies_S10.tsv
sha256: a03fd3f01f52c3401baafa59b04cb5f46c8c4535293e38bc80c3e69281cd55f9
rows: 28175 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_FL_5_frequencies_S11_Fraction1.tsv
sha256: baf6369c3cb1494870aa4d962fe1b107919468f45833ebc62f7e31b50a62aca9
rows: 17660 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_FL_5_frequencies_S11_Fraction2.tsv
sha256: b4391c0ce6d49f2676faa7b44e401b71b00d4e762ddbb41efde74843debc36a3
rows: 18253 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_FL_5_frequencies_S11_Fraction3.tsv
sha256: 49a1cfed6485759ee9d69220b2a5e38ed7e901d4eb03eeeaada3cc4d0e15fb45
rows: 18088 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_FL_5_frequencies_S11_Fraction4.tsv
sha256: bec3e19b0ff908d09e2cfaff13c30fae6cd75280ba4c632ba7f01c2d8b2ae6f9
rows: 9688 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_FL_5_frequencies_S12.tsv
sha256: 3c5a6e270afa56dba8009c1cfd9902e8dd508c05309b374ced5b4ff7fb1963a1
rows: 19004 (one per distinct accepted FL peptide (5-50 residues))
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| peptide | id | peptide sequence |
| accepted_psms | count of PSMs | accepted PSMs for the peptide in this acquisition |

## campi_depth_standardization.tsv
sha256: 7e2073608b4317cb6b0b642d755deebbd5217716029b6a697a8e507a96db03f5
rows: 150 (one per acquisition x arm x length window)
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| arm | id | FL, REF, MG, FLR, FLRM |
| min_length | residues | lower peptide length bound of the population (5 = 5-50 residues, 9 = 9-50 residues) |
| ms2 | count of spectra | MS2 spectra in the acquisition |
| reported_candidates | count of PSMs | PSMs reported by the search |
| accepted_psms | count of PSMs | accepted PSMs |
| peptides | count of peptides | distinct accepted peptides |
| identified_fraction | fraction | accepted_psms / ms2 |
| distinct_per_accepted | fraction | peptides / accepted_psms |
| equal_ms2_draw | count of spectra | k, the smallest ms2 over all rows (38,254) |
| expected_peptides_equal_ms2 | count of peptides | expected distinct peptides when k of ms2 spectra are drawn without replacement (hypergeometric) |
| equal_accepted_psm_draw | count of PSMs | smallest accepted_psms over all rows |
| expected_peptides_equal_accepted_psms | count of peptides | same expectation at equal accepted PSMs (not printed) |

## campi_discrepant_prediction_pooled_summary.tsv
sha256: 42e3fd63a3562062981d16f45d749833913514f3191d0b09f9af867e1472e2d7
rows: 18 (one per population x prediction-reference category)
produced by: campi_near_miss_audit

| column | unit | meaning |
|---|---|---|
| all_candidates_annotated | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| any_lookup | count of spectra | spectra whose prediction retrieved at least one candidate protein (candidate_proteins > 0) |
| category | id | relation of the retained de novo prediction to the reference peptide: exact_IL, adjacent_swap, length_difference, other_disagreement, rearrangement, substitution |
| family_any_compatible | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| family_any_compatible_evaluable | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| family_exclusively_compatible | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| family_exclusively_compatible_evaluable | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| n | count of spectra | spectra in the population and category |
| population | id | all_reference_paired; production_retained_reference_paired (prediction retained by the production filter; used here); production_retained_reference_paired_both_q01 |
| reference_identified_fl_anywhere | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| reference_identified_fl_anywhere_evaluable | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| reference_identified_fl_same_spectrum | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| reference_identified_fl_same_spectrum_evaluable | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| reference_present_final | count of spectra | reference peptide present in the selected final FASTA (searchable) |
| reference_present_final_evaluable | count of spectra | denominator for reference_present_final |
| same_candidate_contains_reference | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| same_candidate_contains_reference_evaluable | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| same_selected_contains_reference | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| same_selected_contains_reference_evaluable | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| selected_family_any_compatible | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| selected_family_any_compatible_evaluable | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| selected_family_exclusively_compatible | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |
| selected_family_exclusively_compatible_evaluable | count of spectra | not used by Figure 2; defined in the upstream A03 report_v2.py |

## campi_gain_availability.tsv
sha256: c9bc2d53c7af91c1346f837e529f95af300a50b3b5bd67bf13a5fb62c76fd74e
rows: 90 (one per acquisition x add-back destination x availability class)
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| destination | id | FLR (FL + REF) or FLRM (FL + REF + MG) |
| availability | id | whether the gained peptide was searchable in the selected FL FASTA: sequence_absent_from_FL_fasta, sequence_present_digest_incompatible, digest_compatible_in_FL_fasta |
| gained_peptides | count of observations | peptides accepted in the destination but not in FL, in this class |
| total_gained | count of observations | all gained peptides of the acquisition |
| total_lost | count of observations | FL peptides not accepted in the destination |
| net_change | count of observations | total_gained - total_lost |

## campi_new_quantified_peptides.tsv
sha256: 609937efda054880dce03eecb1cd2a61522a7850b12fa18c0bb0cbff7e063801
rows: 112 (one per acquisition x destination x category x availability)
produced by: campi_variability

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| destination | id | FLR or FLRM |
| category | id | newly_identified_and_quantified or already_identified_FL_newly_quantified |
| availability | id | as in campi_gain_availability.tsv |
| peptides | count of observations | peptides quantified only in the destination |
| total_destination_only_quantified | count of observations | all destination-only quantified peptides of the acquisition |

## length_by_member_class_a14.tsv
sha256: 76216e6553833ed793497c320d2298ce641bb8b81820a997c54d33dc2a926c9c
rows: 123 (one per member class x peptide length)
produced by: campi_reference_gap

| column | unit | meaning |
|---|---|---|
| class | id | One member, Shared by members (several members) or No member |
| length | residues | peptide length |
| n_peptides | count of peptides | pooled FL peptides of this class and length |
| pct_within_class | % | n_peptides as % of the class (panel h) |

## member_summary.tsv
sha256: bdf6d5f0ded8bc8455f95bd1ac8c4ab52af38294ac7030959a269715873591b8
rows: 80 (one per arm x length window x member)
produced by: campi_reference_gap

| column | unit | meaning |
|---|---|---|
| arm | id | database arm: FL (FastaLake), REF (supplied reference proteomes), MG (metagenome co-assembly proteins), FL_R (FL + REF), FL_R_M (FL + REF + MG) |
| min_length | residues | lower peptide length bound of the population (5 = 5-50 residues, 9 = 9-50 residues) |
| member | id | community member (species name as in the supplied reference) |
| pooled_specific_peptides | count of peptides | distinct pooled peptides unique to this member among the eight proteomes (panel f) |
| acquisitions_ge1 | count of acquisitions | acquisitions with at least one member-specific peptide |
| acquisitions_ge10 | count of acquisitions | acquisitions with at least ten |
| median_specific_peptides | count of peptides | median per-acquisition member-specific peptides |

## metrics_per_acquisition.tsv
sha256: db3499ea77963196082fab113695260bd6d3c4642c71c5811b4ae639186462e9
rows: 75 (one per arm x acquisition)
produced by: campi_database_metrics

| column | unit | meaning |
|---|---|---|
| strategy | id | database arm (FL, REF, MG, FLR = FL + REF, FLRM = FL + REF + MG) |
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| db_records | count of proteins | protein entries in the search database (panel a medians) |
| peptides | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| peptides_per_protein_mean | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| peptides_per_protein_median | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| matched_peaks_median | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| matched_intensity_pct_median | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| psms | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| proteins | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| protein_coverage_pct_median | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| protein_coverage_first_occurrence_pct_median | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| protein_coverage_pct_mean | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| protein_ge50pct_coverage_pct | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| protein_ge2peptide_pct | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| series_coverage_pct_median | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| backbone_coverage_pct_median | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| eligible_coverage_pct_median | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| fragment_rows | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| matched_peaks_mismatches | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| shared_all_psms | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| shared_all_series_coverage_pct_median | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| peptide_length_median | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| mean_assigned_peptides_per_group | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| quantified_groups | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |
| assigned_coverage_pct | see upstream | not used by Figure 2; defined in the upstream benchmark.py (A41) |

## paired_reference_recovery.tsv
sha256: c446dabab7ad3c39c0d055424f60e19c50411eaebe82ee06fae52c54a4f74919
rows: 150 (one per acquisition x arm x length window)
produced by: campi_reference_gap

| column | unit | meaning |
|---|---|---|
| sample | id | CAMPI SIHUMIx acquisition (S01-S12; S11 measured as four fractions) |
| arm | id | database arm: FL (FastaLake), REF (supplied reference proteomes), MG (metagenome co-assembly proteins), FL_R (FL + REF), FL_R_M (FL + REF + MG) |
| min_length | residues | lower peptide length bound of the population (5 = 5-50 residues, 9 = 9-50 residues) |
| peptides | count of peptides | distinct accepted peptides in the acquisition (panel b, n = 15 per arm) |
| reference | count of peptides | distinct peptides of the REF search of the same acquisition |
| overlap | count of peptides | peptides shared with the REF search |
| gained | count of peptides | peptides not in the REF search |
| lost | count of peptides | REF peptides not found by this arm |
| count_ratio | dimensionless | peptides / reference |
| reference_recovery | fraction | overlap / reference (within-acquisition recovery) |

## recovery_summary.tsv
sha256: 74658e932aea9f40a62245b6f3153f578289343182a02581ecca83736c060b13
rows: 10 (one per arm x length window)
produced by: campi_reference_gap

| column | unit | meaning |
|---|---|---|
| arm | id | database arm: FL (FastaLake), REF (supplied reference proteomes), MG (metagenome co-assembly proteins), FL_R (FL + REF), FL_R_M (FL + REF + MG) |
| min_length | residues | lower peptide length bound of the population (5 = 5-50 residues, 9 = 9-50 residues) |
| acquisitions | count of acquisitions | acquisitions pooled |
| pooled_peptides | count of peptides | distinct peptides pooled over 15 acquisitions |
| pooled_reference | count of peptides | distinct pooled REF peptides (76,127 at 5-50) |
| pooled_overlap | count of peptides | pooled peptides also in the pooled REF set |
| pooled_reference_recovery | fraction | pooled_overlap / pooled_reference (panel e right column) |
| median_acquisition_count_ratio | dimensionless | median over acquisitions of peptides / reference |
| median_acquisition_reference_recovery | fraction | median over acquisitions of reference_recovery |
| known_member_specific | count of peptides | pooled peptides found in exactly one of the eight member proteomes |
| known_member_shared | count of peptides | pooled peptides found in two or more member proteomes |
| absent_eight_member_reference | count of peptides | pooled peptides found in none of the member proteomes |

## reservoir_protein_presence_by_proteome.tsv
sha256: 032c84a3383dbd16ecccb95d586c3fc3b2f131d95c4f1757d4ab0983005c68a8
rows: 16 (one per member proteome, plus one shuffled control per member)
produced by: campi_reservoir_containment

| column | unit | meaning |
|---|---|---|
| label | id | member reference proteome; _SHUFFLED = shuffled-sequence control |
| n_proteins | count of proteins | reference records |
| n_unique_il_sequences | count of sequences | distinct sequences after I->L |
| n_present_exact | count of proteins | records whose full I/L-equivalent sequence is in the 375,333,654-protein reservoir |
| fraction | fraction | n_present_exact / n_proteins, rounded to 4 decimals |
| n_present_exact_raw | count of proteins | records present without I/L equivalence |
| n_present_il_variant_only | count of proteins | records present only via I/L equivalence |
| n_kept_UHGG | count of proteins | present records found in UHGG |
| n_kept_GMGC | count of proteins | found in GMGC |
| n_kept_GMSC | count of proteins | found in GMSC |
| n_kept_UNIPROT_SP | count of proteins | found in UniProtKB/Swiss-Prot |
| n_kept_UNIPROT_TR | count of proteins | found in UniProtKB/TrEMBL |

## stage_summary.tsv
sha256: bfd51133b19b7c4bb6a5594ae4efe84aac9f8b6a1c632ff9c25f9874b587bd05
rows: 12 (one per length window x stage)
produced by: campi_reference_gap

| column | unit | meaning |
|---|---|---|
| min_length | residues | lower peptide length bound of the population (5 = 5-50 residues, 9 = 9-50 residues) |
| stage | id | stage at which a REF peptide-acquisition observation missed by FL was lost (recovered_FL = not lost); the first three stages mean the sequence was not in the selected FASTA |
| peptide_acquisition_observations | count of observations | peptide-acquisition observations at this stage, summed over 15 acquisitions (panel d) |
| no_direct_prediction_observations | count of observations | of these, observations whose spectrum had no exported de novo prediction |

## upset_intersections_a14.tsv
(5-50 residues; the v28f panel g. v30 panel g draws upset_intersections_FL_9_50.tsv.)
sha256: 7f0405f42961e168297824f1118e573453a6edc3fe3cb464792c063ac79c20b5
rows: 230 (one per member combination)
produced by: campi_reference_gap

| column | unit | meaning |
|---|---|---|
| rank | dimensionless | rank of the combination by size |
| members | id | member combination, '|' separated; '(none of the eight)' for degree 0 |
| degree | count of members | members in the combination |
| n_peptides | count of peptides | pooled FL peptides whose member set is exactly this combination |
| drawn | boolean | True for the 14 combinations drawn in panel g |

## upset_intersections_FL_9_50.tsv
sha256: 1f7e4d0e4b01a59eb11a9059641a26274134e48adf83f85db032092c3004035b
rows: 126 (one per member combination)
produced by: recount of the 15 membership tables (`source_data/SFig03/CAMPI_membership_<sample>.tsv`), FL arm, peptides of 9-50 residues; v30 len950 rebuild (`panels.upset_at`). The notebook recounts it and asserts equality (panel g, v30).

| column | unit | meaning |
|---|---|---|
| rank | dimensionless | rank of the combination by size |
| members | id | member combination, '|' separated; '(none of the eight)' for degree 0 |
| degree | count of members | members in the combination |
| n_peptides | count of peptides | pooled 9-50-residue FL peptides whose member set is exactly this combination |
| drawn | boolean | True for the 14 combinations drawn in panel g (8 single-member sets and the 6 largest shared sets) |

### reservoir_member_presence.tsv

Exact presence of each CAMPI member reference proteome in the 375,333,654-entry reservoir (Results: 91.3–100%; Source Data map RES-24). Byte-identical copy of `protein_presence_by_proteome.tsv` from `pipelines/campi_reservoir_containment` (sha256 032c84a3...). Rows ending `_SHUFFLED` are the shuffled-proteome negative controls (all 0).

| column | unit | meaning |
|---|---|---|
| label | text | community member (or its shuffled control) |
| n_proteins | proteins | proteins in the supplied member reference proteome |
| n_unique_il_sequences | sequences | distinct sequences after I/L collapse |
| n_present_exact | proteins | proteins whose I/L-equivalent sequence occurs exactly in the reservoir |
| fraction | fraction | n_present_exact / n_proteins |
| n_present_exact_raw | proteins | exact presence without I/L collapse |
| n_present_il_variant_only | proteins | present only as an I/L variant |
| n_kept_UHGG … n_kept_UNIPROT_TR | proteins | present proteins by the reservoir source that holds them |
