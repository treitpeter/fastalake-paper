# peptide_centric_downstream

Internal ID: A33 (previous folder name A33_peptide_centric_downstream)

Does the biological conclusion (leader species, function, region contrasts, own-stool similarity) hold when proteins are grouped peptide-centrically, for HSTOOL, PREDICT and CapScan?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `build_design_347.py` | call chain of the used outputs |
| `build_gated_joint.py` | call chain of the used outputs |
| `build_supp_9_10_19.py` | call chain of the used outputs |
| `collector.py` | call chain of the used outputs |
| `common_pepcentric.py` | call chain of the used outputs |
| `compare_control.py` | call chain of the used outputs |
| `derive_hstool_pc.py` | call chain of the used outputs |
| `fig5ef_run.py` | call chain of the used outputs |
| `hstool_pc.py` | call chain of the used outputs |
| `make_leader_tables.py` | call chain of the used outputs |
| `scripts__build_design_347.py` | call chain of the used outputs |
| `search_pc.py` | call chain of the used outputs |

Launchers: `slurm/hstool.sbatch`.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: Merged A33_hstool_leading_species and C2_capscan_participant15_contributors (same A33 peptide-centric code | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/Fig3, figures/Fig5, figures/SFig12, text_statistics/supplementary_methods_6_7, text_statistics/supplementary_methods_8_11. What is used: Merged A33_hstool_leading_species and C2_capscan_participant15_contributors (same A33 peptide-centric code; C2 is its 347-roster rerun). HSTOOL simplified pairwise protein-group table and design (Fig3); 347-acquisition run347 mosaic_donor_regions.tsv, mosaic_abundance.tsv, participant_leaders_*.tsv, contributor_cells.tsv, function_species.npz (Fig5 e,f); HSTOOL function x design leader diversity with rarefaction and CapScan function x participant(x region) leading-species cells under both grouping rules (SFig12, Slurm 6058098); representative_presence.tsv 43.2% (line 751)

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

| jobs | cores | requested memory | peak RSS | longest task | summed task time | states |
|---|---|---|---|---|---|---|
| 11 tasks | 2 | 16G, 2450M, 32G, 4G | 2.3 GB | 3 min 43 s | 15 min 24 s | COMPLETED 9, FAILED 2 |

Partitions: p.hpcl92. Job ids: 5973723, 5973693, 5973716, 5973758, 6058096, 6058097, 6058098.

## Choices

Approach: peptide-centric collector, gated joint, leader tables and the CapScan 347-acquisition rerun (fig5e_f).

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| grouping | peptide-centric (original and simplified variants compared) | tests sensitivity of downstream biology to the grouping rule | both variants are written to *__original and *__simplified |
| CapScan roster | fixed 347 acquisitions, 15 participants | roster fixed by A80 | 331-acquisition control run (job 6058096) |
| gating | gated joint of species and KO evidence | code in build_gated_joint.py | ungated variant not run here |

## Known limitations

See `OPEN_ITEMS.md`.

## Direct-assignment recount of Fig. 5f (v30)

The v30 manuscript removed the original protein grouping; production is direct peptide
assignment. The published Fig. 5f leading-species agreement had been restricted to
functions eligible under both grouping rules and comparisons resolved under both.
`strict_subset_direct_recount.py` recomputes it from the `fig5ef_run.py` outputs
(`results/contributor_method_cells.tsv`, `results/contributor_eligibility.tsv`) with the
direct-assignment criteria only. Same estimator: comparisons averaged within participant
and region, then over regions, then over the 15 participants (unit: participant).

| definition | mean | participants | comparisons | functions |
|---|---|---|---|---|
| control, all eligible (both rules) | 45.58% | 15 | 51,484 | 1,335 |
| control, strict (both rules) | 59.99% | 15 | 18,137 | 1,135 |
| direct, all eligible | 46.43% | 15 | 58,208 | 1,442 |
| direct, strict | 57.57% | 15 | 24,537 | 1,264 |

The controls reproduce the published 45.6% and 60.0%; the v30 values are 46.4% and 57.6%.
Launcher: `slurm/strict_subset_direct_recount.sh` (environment variables
`FASTALAKE_FIG5EF_RUN_DIR`, `FASTALAKE_OUT_DIR`, `FASTALAKE_PYTHON`). Resources: 2 cores,
8 GB requested, peak RSS 3.4 MB, 2 s. Job id and hashes in RECEIPTS.md.
