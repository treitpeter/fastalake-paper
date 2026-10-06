# Source data: Figure 5

CapScan, 347 LC-MS/MS acquisitions from 308 specimens of 15 participants. Region ids:
`duodenum` = Region 1, `jejunum` = Region 2, `ileum` = Region 3, `colon` = Region 4
(device-defined sampling regions), `stool` = stool. Participants are study pseudonyms 1-15.
All tables replace the 14-participant / 331-acquisition sheets of the 4 October Source Data
package (see `figures/Fig5/PROVENANCE.md`). Withheld from the design table: device ids,
swallow and recovery dates and times, transit time, pH, O2, DNA and protein yields and prior
antibiotic exposure. No Figure 5 computation uses them.

## design_347.tsv
sha256: d825accc884bdba5ea1ed7cd5d0606ee93ffa30cc13b3751e7c56c2a4444829e
rows: 347 (one row per acquisition)
produced by: capscan_participant_roster

| column | unit | meaning |
|---|---|---|
| sample_id | id | raw-file name of the acquisition |
| participant | id | participant pseudonym (1-15) |
| region | id | sampling group (duodenum, jejunum, ileum, colon = regions 1-4; stool) |
| specimen_id | id | physical specimen; acquisitions of one specimen share it (308 specimens) |
| og_profile_positive | boolean | acquisition has a positive non-host OG profile (346 of 347); roster of panels b-d |

## taxonomy_signal.tsv
sha256: 6e4160e3bbe601802c3ee64b7c573fb1f0f9408c37e74b3f344c182ca3e8c568
rows: 347 (one row per acquisition)
produced by: capscan_functional_profiles

| column | unit | meaning |
|---|---|---|
| sample | id | acquisition (= design sample_id) |
| total_LFQ | LFQ intensity | total accepted LFQ intensity |
| Human-assigned | LFQ intensity | intensity of features assigned to human proteins |
| Microbial above species | LFQ intensity | microbial intensity resolved above species rank |
| Microbial species | LFQ intensity | microbial intensity resolved at species rank |
| Unresolved origin | LFQ intensity | intensity of unresolved origin |

## origin_fraction_cells.tsv
sha256: 3d73defc53ccab9c2f717638c4dc96365256644d1dcd34dd107290bcffc4a7c6
rows: 75 (one row per participant and sampling group)
produced by: capscan_functional_profiles (recomputed by `panels.stat_b_region_tests`)

| column | unit | meaning |
|---|---|---|
| participant | id | participant |
| region | id | sampling group |
| Human-assigned, Microbial above species, Microbial species, Unresolved origin | fraction of total LFQ | origin share, averaged within specimen then participant-region |

## four_region_tests.tsv
sha256: 885df8ef9cef703756cab97240e7662963c3bc6d7dc8f8e555ab842fa90e9736
rows: 4 (one row per capsule region versus stool)
produced by: capscan_functional_profiles (recomputed by `panels.stat_b_region_tests`)

| column | unit | meaning |
|---|---|---|
| region | id | capsule region compared with stool |
| n | count of participants | paired participants |
| mean_region, mean_stool | % of total LFQ | mean human-assigned share |
| mean_difference | percentage points | mean paired difference, region minus stool |
| ci_low, ci_high | percentage points | 95% percentile-bootstrap interval |
| p | dimensionless | exact two-sided sign-flip P |
| q | dimensionless | Benjamini-Hochberg q over the four regions |

## similarity_summary.tsv
sha256: ac61627e0eb03c7893603819a7f7b590172cfec50b8d59eb7c30699699f69a39
rows: 3 (one row per OG prevalence filter)
produced by: capscan_functional_profiles (recomputed by `panels.stat_cd_profiles`)

| column | unit | meaning |
|---|---|---|
| prevalence | fraction of profiles | minimum fraction of the 75 profiles in which an OG is present |
| participants | count of participants | participants |
| features, positive_features | count of OGs | OGs retained; OGs with positive signal |
| own_closer | count of participants | participants closer to own stool |
| mean_own, mean_other | Euclidean distance (Hellinger profiles) | mean distance to own stool, to other stools |
| mean_advantage | Euclidean distance | mean of (other minus own) |
| ci_low, ci_high | Euclidean distance | 95% percentile-bootstrap interval |
| p | dimensionless | two-sided permutation P (99,999 permutations) |

## own_stool_participants.tsv
sha256: 67f22f573d418bdee6e7b2a7ed0bc1ccc2fc140ee5793ab10db2455f6002ff0e
rows: 45 (one row per participant and prevalence filter)
produced by: capscan_functional_profiles

| column | unit | meaning |
|---|---|---|
| prevalence | fraction of profiles | prevalence filter |
| participant | id | participant |
| own_stool_distance | Euclidean distance | mean distance of the four capsule profiles to own stool |
| other_stool_mean_distance | Euclidean distance | mean distance to the 14 other stools |

## orthogonal_partition.tsv
sha256: 035057c305e7196665b92f91525f4038176d7d41f25a4f4b09cfab51a1921319
rows: 24 (one row per population, component and prevalence filter)
produced by: capscan_functional_profiles

| column | unit | meaning |
|---|---|---|
| population | id | All five regions (75 profiles) or Capsules only (60 profiles) |
| component | id | Participant, Capsules versus stool, Between capsule regions, Residual |
| fraction | fraction of total sum of squares | share of profile variation |
| prevalence | fraction of profiles | prevalence filter |

## mosaic_donor_regions.tsv
sha256: c3b50d0be5b0b4d6256e96c56c675d411db75ae0f9cc9aed62d4dcbbb3404d9e
rows: 1800 (one row per function, participant and sampling group; 24 x 15 x 5)
produced by: peptide_centric_downstream

| column | unit | meaning |
|---|---|---|
| row | id | display row of the function (1-24) |
| KO | id | KEGG ortholog |
| participant | id | participant |
| region | id | sampling group |
| species | id | leading annotated species (GTDB name), empty if not resolved |
| phylum | id | phylum of the leading species |
| resolved_fraction | fraction | species-resolved share of the function's non-host signal |
| top_share | fraction | leading species' share of species-resolved signal |
| status | id | Resolved, or the reason no leader is shown |

## mosaic_abundance.tsv
sha256: 8bf92d581260a971c3aff4824aee008ff7d64d2804789ebbb829dc467a71879f
rows: 24 (one row per displayed function)
produced by: peptide_centric_downstream

| column | unit | meaning |
|---|---|---|
| row | id | display row |
| KO | id | KEGG ortholog |
| mean_share_percent | % of total LFQ | mean share of the function in total LFQ intensity |

## heatmap_function_key.tsv
sha256: d0ee5641797fd588ba159362019515b340efed901aa03490bd241ea74ff0493c
rows: 24 (one row per displayed function)
produced by: peptide_centric_downstream (fixed display key)

| column | unit | meaning |
|---|---|---|
| row | id | row in the key file |
| KO | id | KEGG ortholog |
| display_name | id | label printed in panel e |
| official_name | id | KEGG name |
| url | id | KEGG entry |

## participant_leaders_resolved.tsv, participant_leaders_strict_both_methods.tsv
sha256: 0db1bb548130d1d2443ce353dd0657a9c67b5f8664172b7b75f3f0126178d478 and de2ce4d3935611877671b7414876e31cd9152a31fe3bfc106de1cf5a189e3746
rows: 15 each (one row per participant)
produced by: peptide_centric_downstream (v30 recount, direct peptide assignment only)

v30 tables. The original protein grouping was removed from the manuscript (v30), so these tables
hold direct peptide assignment only. `participant_leaders_resolved.tsv`: all functions eligible
under direct assignment, comparisons with a resolved leader in both the intestinal region and
stool (58,208 comparisons, 1,442 functions). `participant_leaders_strict_both_methods.tsv`: the
strict subset (species-resolved fraction >= 0.20 and leader margin >= 0.05 in both the region and
the stool cell) under direct assignment only (24,537 comparisons, 1,264 functions); the file name
is kept from v28f so that the plotting code is unchanged. Comparisons are averaged within
participant and region, then over regions within participant.

| column | unit | meaning |
|---|---|---|
| participant | id | participant |
| agreement_simplified | fraction | share of function-region comparisons in which the capsule leader equals the own-stool leader, direct peptide assignment |

## contributor_cells.tsv
sha256: 6a704174763ce561aa2189fba7b137f3dc480edeebf7f52bdd16cabbb492dd84
rows: 21168 (one row per eligible function and participant; 1,467 functions)
produced by: peptide_centric_downstream

| column | unit | meaning |
|---|---|---|
| KO | id | KEGG ortholog |
| participant | id | participant |
| leading_species | id | leading annotated species |
| leader_share | fraction | leader's share of species-resolved signal |
| gap | fraction | leader share minus second share |
| species_resolved_fraction | fraction | species-annotated share of the function's non-host signal |
| function_share_total | fraction | function's share of total LFQ |

## evenness_per_acquisition.tsv
sha256: 3840cef3b6c9cd22c0f6650c86f297137872a16d8c8ced1504c38f38bc515c31
rows: 347 (one row per acquisition)
produced by: capscan_dominance_evenness (recomputed by `panels.stat_evenness_per_acquisition`)

| column | unit | meaning |
|---|---|---|
| sample_id, participant, region, specimen_id | id | as in design_347.tsv |
| source | id | species annotation source (a33) |
| compartment | id | fluid (regions 1-4) or stool |
| top1 | fraction | leading named species' share of named-species signal |
| pielou | dimensionless | Pielou evenness of named-species signal |
| richness | count of species | named species with positive signal |
| resolved_mass | LFQ intensity | named-species signal |
| species_resolved_fraction | fraction | named-species signal over non-host functional signal |

## evenness_summary_tests.tsv
sha256: c78776e1691342ca4c497e7039d9eb7559f10527c6f0ac549a0250faab3f8323
rows: 2 (one row per metric)
produced by: capscan_dominance_evenness (recomputed by `panels.stat_evenness`)

| column | unit | meaning |
|---|---|---|
| source, metric | id | annotation source; top1 or pielou |
| acq_n_fluid, acq_n_stool, acq_n_undefined | count of acquisitions | acquisitions used and undefined |
| acq_median_fluid, acq_median_stool | as metric | acquisition-level medians |
| paired_participants | count of participants | participants with fluid and stool |
| participant_median_fluid, participant_median_stool | as metric | median over participants of participant medians |
| median_paired_difference, mean_paired_difference | as metric | fluid minus stool |
| participants_fluid_gt_stool | count of participants | participants with fluid above stool |
| wilcoxon_W, wilcoxon_p_exact | dimensionless | exact Wilcoxon signed-rank statistic and two-sided P |
| signflip_p_exact, signflip_patterns | dimensionless; count | exact sign-flip P and number of patterns |

## native_correlations.tsv
sha256: 0004f35d522422580e1d7cce3bbcbc2d2967a8fbc8af28668b73408d52305375
rows: 12 (one row per feature filter and layer pair)
produced by: capscan_crossomics (rows of `sensitivity_all_features` recomputed by `panels.stat_g_native`)

| column | unit | meaning |
|---|---|---|
| cohort, analysis | id | CAPSCAN; feature filter (sensitivity_all_features = all features, shown in Fig. 5g) |
| layer_a, layer_b | id | DNA, Non-host protein, Host protein, Metabolites |
| rho | dimensionless | Spearman correlation of distance upper triangles |
| p, q | dimensionless | permutation P and Benjamini-Hochberg q (6 pairs) |
| permutations, n_profiles, n_participants | count | 1,999; 65; 13 |

## functional_correlations.tsv
sha256: 4c3468fc11199d32c9dd99fbda823bb536f624d01102df87562b6e29906de3f8
rows: 64 (one row per scope, representation, layer pair and prevalence)
produced by: capscan_crossomics (rows of scope `matched12`, prevalence 0, recomputed by `panels.stat_g_functional`)

| column | unit | meaning |
|---|---|---|
| cohort, scope | id | CAPSCAN; matched12 = 13 matched participants (label from the 12-participant version), displayed_mixed12_14 = 15/13 participants |
| level | id | OG, KO or Pathway |
| layer_a, layer_b | id | metaG, microbial_metaP, host_metaP, metabolites |
| prevalence | fraction of profiles | feature prevalence filter |
| n_units, n_donors | count | profiles; participants |
| rho, p, q | dimensionless | Spearman rho, permutation P, BH q (16 per scope and prevalence) |

## native_distances_long.tsv, functional_distances_long.tsv
sha256: 07bda50dc99acf5fc5e6a77ec91d9fd71bc90451a924c2413cac193de6de4530 and 006d3a0fa3201603f15f329846c11cc0b4f0f65da3584483eb51c363523353ad
rows: 8320 and 24960 (one row per profile pair, i < j, per layer; 2,080 pairs of 65 profiles)
produced by: capscan_crossomics (upper triangles of the producer's distance matrices; read with `float_precision='round_trip'`)

| column | unit | meaning |
|---|---|---|
| analysis / scope, prevalence, level | id | which matrix (native all features; functional matched12, prevalence 0, OG/KO/Pathway) |
| layer | id | data layer |
| participant_i, region_i, participant_j, region_j | id | the two participant-region profiles |
| distance | Euclidean distance | distance between region-adjusted transformed profiles |

## Larger files (Zenodo)

Read from `build/zenodo/Fig5/`; `panels.zenodo_ready` checks each sha256.

| file | bytes | sha256 | content | used by |
|---|---|---|---|---|
| OG_other_unresolved.npz | 42421948 | e277cb600c7ee2b6f6e94b0bbd0f86a32a90878cd1b082045c77ed55b7b946a1 | sparse acquisition x OG non-host functional intensities (347 acquisitions) | panels.stat_cd_profiles |
| OG_other_unresolved_axes.json | 1809017 | 3c3544a93f32ad1eacd55e3f1a9758265cf724fa3004aed0f13c3999e81af955 | row (acquisition) and column (OG) labels | panels.stat_cd_profiles |
| function_species.npz | 21910543 | 4182b86aeb96de12e95c0572c039445a00e3fd70abe91b3addcf7c585e76817d | KO x species intensities per acquisition | panels.stat_evenness_per_acquisition |
