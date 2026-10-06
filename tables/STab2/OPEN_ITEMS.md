# Open items: STab2

- MISMATCH: v28d line 743 cites Supplementary Table 2 for "2 of 376,412" and "12 of 411,354", but neither denominator occurs in any cell of the published workbook; the values come from the healthy-stool pooling experiment (`source_data/STab2/hstool_pooling_per_acquisition.tsv`, package file `Figure3/per_acquisition.tsv`). For v28f either the sentence cites that Source Data file, or the table gains a sheet with these 360 rows (NUMBERS.tsv rows 1-4; notebook cell Methods line 743).
- AUTHOR: the full-route control experiment covers 710 acquisitions (CapScan 331); the v28f roster is 347 CapScan acquisitions and 726 FASTAs. It was not rerun for the 16 acquisitions of CapScan participant 15, and no text says that it was not rerun by design. Decide whether the legend (v28d line 477) or Methods states "710 acquisitions (before participant 15 was added)" or the experiment is extended (Read_me sheet, Cohort_overview n column).
- AUTHOR: the `fusion_*` sheets are the 2026-09-08 fusion run (710 acquisitions, two seeds). v28d Supplementary Figure 5 and Methods line 715 use a later reconciliation (unit SFig05, healthy_stool_donor_entrapment and campi_fusion_375m). Keep these sheets in Supplementary Table 2, relabel them as superseded, or remove them for v28f (sheets fusion_Read_me .. fusion_All_acquisitions; fusion_Read_me still names "SupplementaryFigure6c–f").
- AUTHOR: the registry note for this workbook (TABLE_REGISTRY.tsv, key CONTROL_EXPOSURE) says the earlier CAMPI-reservoir stratum is superseded for current CAMPI inference; the workbook itself carries no such note (Read_me sheet).
- TEXT-GAP: the manuscript does not state how acquisition quartiles in Quantile_summary were computed; the code uses Python `statistics.quantiles(method="exclusive")` across acquisitions with undefined fractions excluded (STATS.tsv row 2).
- TEXT-GAP: the manuscript does not state the seed of the shuffled-control construction (20260906, experiment.json `control_source_stats.seed`) or the fusion control seeds (20260907, 20260908; fusion sheets column `seed`) (STATS.tsv row 1).
- UNTRACED: the bytes of the published file (`4768d39f...`) differ from the `package_v15.py` output (`5eb86111...`); the step between them was not traced to a script. Cell values of the published file equal the rebuild (ASSEMBLY.md step 4).
- UNTRACED: interpreter and package versions of the original producer runs (build_package.py, package_edit3.py, package_v15.py, reduce.py) were not recorded (PROVENANCE.md, Panels table).
- UNTRACED: Slurm job ids of the 2026-09-08 fusion run (A92) and of the CAMPI independent-selection run (A93) were not found in the archived inputs (UPSTREAM.json).
- MISMATCH: `fusion_cohort_threshold_summary` and `fusion_matched_seed_summary` are byte-identical source tables (sha256 bb656c79...), so two sheets of the workbook repeat the same 80 rows under different names (source_data/STab2/README.md).
- shared requests: PIPELINE: the four upstream analyses (cohort_full_route_controls, cohort_fusion_controls, campi_independent_selection, healthy_stool_pooling_controls) have no `pipelines/` directory; their original code contains absolute paths (experiment.json, sbatch log paths) and was not moved in this unit. Pipeline ids were chosen by this unit and should be checked against ids other units create.
- shared requests: the notebook does not call `fl.init()` because the unit draws nothing; the setup cell otherwise follows CONVENTIONS.md 4.2.

## v30 reference (added 2026-10-06)

- none: every row is located in the v30 text and agrees.
