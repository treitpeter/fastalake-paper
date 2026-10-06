# Source data: Supplementary Figure 12

Tables read by `figures/SFig12/SFig12.ipynb`. Origins and hashes are in
`figures/SFig12/PROVENANCE.md`. Every table is on the v28f roster (347 CapScan acquisitions,
15 participants); the 331-acquisition versions in the earlier Source Data package are replaced.

Two tables the figure also reads live under `source_data/Fig5/` and are not duplicated here:
`contributor_cells.tsv` (panels c, d and the panel f recomputation; sha256
`6a704174763ce561aa2189fba7b137f3dc480edeebf7f52bdd16cabbb492dd84`) and
`participant_leaders_resolved.tsv` / `participant_leaders_strict_both_methods.tsv` (the numbers
on v28d line 157). Their data dictionaries are in `source_data/Fig5/README.md`.

Sample names are acquisition (raw file) names. No clinical covariate is in any table here.

## HSTOOL_KO_all.tsv
sha256: e443ff28bb1ee485de65c55489f5f98d9b4c86037d3f392e8385edf7dc8248af
rows: 3254 (one row per KO, plus one `FunctionallyUnassigned` row); 58 acquisition columns
produced by: peptide_centric_functional_allocation

| column | unit | meaning |
|---|---|---|
| feature | id | KEGG ortholog (`K#####`), or `FunctionallyUnassigned` for accepted feature intensity not allocated to any KO |
| one column per acquisition | LFQ intensity | intensity allocated to that KO in that acquisition (peptide-centric equal split over distinct eligible KOs); the column sum is the acquisition's accepted feature intensity |

## PREDICT_KO_all.tsv
sha256: 663e1213de6c7fcdec2304bef076d420a4f5bf7320908c4d369136e16f454892
rows: 7450 (one row per KO, plus `FunctionallyUnassigned`); 306 acquisition columns
produced by: peptide_centric_functional_allocation (repaired PREDICT lineage)

Columns as in `HSTOOL_KO_all.tsv`.

## CAPSCAN_KO_all.tsv
sha256: 1fcc2eb3fa9a634e9ef2f4f6f7630b8d835a3743c7afcf027312373af507eb23
rows: 7818 (one row per KO, plus `FunctionallyUnassigned`); 346 acquisition columns (the 347-acquisition roster less one acquisition with no accepted features)
produced by: peptide_centric_functional_allocation, rerun on the 347 roster

Columns as in `HSTOOL_KO_all.tsv`.

## KO_allocation_coverage.tsv
sha256: 05a9a9d09028d3835706be3fad3d4f71903c59288cba70056fd92744b1648b01
rows: 710 (one row per acquisition: 58 HSTOOL, 306 PREDICT, 346 CAPSCAN)
produced by: `figures/SFig12/panels.py:coverage_from_matrices` from the three matrices above (the notebook re-derives it byte for byte)

| column | unit | meaning |
|---|---|---|
| cohort | id | HSTOOL (healthy stool), PREDICT or CAPSCAN |
| sample | id | acquisition |
| total_intensity | LFQ intensity | accepted feature intensity (column sum of the KO matrix, unassigned remainder included) |
| unassigned_intensity | LFQ intensity | the `FunctionallyUnassigned` row |
| allocated_fraction | fraction | 1 - unassigned_intensity / total_intensity (panel a, shown as %) |

## healthy_stool_leader_diversity.tsv
sha256: 93ed4ead5ea85f8d4fd8bf09750bb07dc8a66c0fe8565129f20c85988f496147
rows: 3862 (one row per function x acquisition design with a resolved leader in at least 80% of that design's acquisitions)
produced by: peptide_centric_functional_allocation (derive_hstool_pc.py, function_tier_summary.tsv)

| column | unit | meaning |
|---|---|---|
| tier | id | acquisition design (TechnicalReplicates, WorkflowReplicates, OnePiece, Intra, Interindividual) |
| ko | id | KEGG ortholog |
| label | text | KO name |
| eligible_in_all_five_designs | 0/1 | 1 if the function is eligible in every design |
| n_observed_resolved | count of acquisitions | acquisitions of the design with a resolved leading species for this function |
| n_distinct_leaders | count of species | distinct leading species among them |
| expected_distinct_leaders_in_8_resolved_acquisitions | count of species | analytical rarefaction: sum over leading species of 1 - C(n - v, 8) / C(n, 8), n = n_observed_resolved, v = acquisitions led by the species (panel b) |
| mean_intensity_pct_equal_three_designs | % | mean normalised function intensity over the three biological designs (ordering only) |
| species_names | text | the distinct leading species, `; `-separated |
| main_top12 | 0/1 | 1 if among the twelve functions shown in the main-text healthy-stool figure |

## healthy_stool_grouping_leader_sensitivity.tsv
sha256: fb4e1e3c34a233f627aad13fd78e33af77c74d9b785959ce7f934cbd8d589716
rows: 16 (grouping rule x design x subset x filter)
produced by: peptide_centric_functional_allocation (S10 summary)

| column | unit | meaning |
|---|---|---|
| method | id | `original` protein grouping or `simplified` (direct peptide assignment); the v30 figure draws only `simplified` |
| tier | id | TechnicalReplicates or WorkflowReplicates |
| subset | id | `all_eligible` (plotted in panel e) or `author_top12` |
| filter | id | `primary` or `well_resolved` |
| n_kos | count of functions | functions compared |
| n_pairs | count of acquisition pairs x functions | pairwise comparisons |
| equal_KO_agreement | fraction | pairwise same-leader agreement averaged equally over functions (panel e, shown as %) |
| pooled_pair_agreement | fraction | same-leader fraction over all pairs pooled |

## CapScan_grouping_leader_sensitivity.tsv
sha256: 593d0251c4db0947172fd0437b342908546ea1910cac6cdfa18dd801fe5af9b2
rows: 8 (grouping rule x eligibility rule)
produced by: peptide_centric_functional_allocation (check_s10c_eligibility.py four_rules on the 347-roster selections); recomputed in the notebook

| column | unit | meaning |
|---|---|---|
| source | id | `A33_CAPSCAN__simplified` (direct peptide assignment) or `A33_CAPSCAN__original` (the v30 figure draws only the `simplified` rows); the `A33_` prefix is the producer's label and is kept because the plotting code selects rows by it; the rows are from the 347-roster rerun |
| rule | id | primary, strict_pair_filter, strict_donor_cell_eligibility, strict_paired_donor_eligibility |
| donor_balanced_agreement | fraction | same-leader fraction averaged within participant-region, then participant, then over participants (panel f, shown as %) |
| pooled_cell_agreement | fraction | same-leader fraction over all pairs pooled |
| n_pairs | count of function x participant x region pairs | pairs compared |
| n_kos | count of functions | functions compared |
| n_donors | count of participants | participants |
| n_donor_regions | count of participant x region | participant-region combinations |

## capscan_contributor_region_cells_direct.tsv
sha256: e7cb87a78397858a6487ab4c81475f3dfdc7b465cf9fe8643963af4f6e6ee78c
rows: 104348 (one row per function x participant x sampling region; direct peptide assignment)
produced by: peptide_centric_functional_allocation (C2 participant-15 rerun, fig5ef_run.py)

| column | unit | meaning |
|---|---|---|
| KO | id | KEGG ortholog |
| participant | id | CapScan participant (1-15) |
| region | id | duodenum, jejunum, ileum, colon (intestinal sampling regions) or stool |
| leading_species | text | largest carrier-supported species contribution (lexical tie-break) |
| leader_share | fraction | leading species' share of species-resolved signal |
| gap | fraction | leader share minus runner-up share, of species-resolved signal |
| species_resolved_fraction | fraction | species-resolved share of non-host functional signal |

