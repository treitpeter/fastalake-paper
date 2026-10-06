# capscan_participant_roster

Internal ID: A80 (previous folder name A80_capscan_participant15_roster)

Which 347 CapScan acquisitions (15 participants, 16 added acquisitions) form the fixed roster used by every CapScan statistic?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `build_fixed_roster.py` | call chain of the used outputs |
| `make_design347.py` | call chain of the used outputs |

Launchers: none, the step ran interactively or inside another launcher.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: fixed 347-acquisition CapScan roster (15 participants): design_347.tsv (participant, region, specimen) and participant15_added.tsv (16 added acquisitions and jo | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/Fig5, figures/Fig6, figures/SFig15, tables/SData1, text_statistics/methods, text_statistics/supplementary_methods_8_11. What is used: fixed 347-acquisition CapScan roster (15 participants): design_347.tsv (participant, region, specimen) and participant15_added.tsv (16 added acquisitions and join evidence); Fig5 a counts, roster for every CapScan statistic, SFig15 CapScan task selection, SData1, Methods line 269

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

No Slurm job ids are recorded for this step (see OPEN_ITEMS.md). It ran without a recorded accounting entry.

## Choices

Approach: fixed roster, design table (participant, region, specimen) and join evidence for the 16 added acquisitions.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| roster size | 347 acquisitions, 15 participants | participant 15 added to the earlier 331 | 331 kept as control |
| join | by recorded specimen and region identifiers, evidence written per row | no neighbour substitution | none |

## Known limitations

See `OPEN_ITEMS.md`.

## CapScan DNA count reconciliation

The manuscript quotes five CapScan DNA counts. They count different things and are
reconciled here. Unit throughout: one specimen = one BioSample = one ENA sample alias whose
first token is the CapScan `Main` specimen id. Proteomic acquisitions link to specimens
through `aliquot_main_id` (331-acquisition roster, `source_data/SData1/v28d_capscan_acquisitions.tsv`)
and `workbook_Main` (16 participant-15 additions, `participant15_added.tsv`).

| count | what it is | source |
|---|---|---|
| 335 | assembly FASTAs in the CapScan reservoir, one per specimen, independent of the proteomic roster | `LAB_ROOT/canonical/publication_data/cohorts/CAPSCAN/0_lake/CAPSCAN_lake_stats.txt`; `source_data/methods/reservoir_sources.tsv` |
| 302 of 335 | specimens with a public WGS run in PRJNA822660; the other 33 (32 capsule, 1 stool) are amplicon-only in ENA and were dropped in the original study | `LAB_ROOT/analyses/capscan_metag_availability_2026-09-09/ena_run_table_2026-09-10.tsv` |
| 301 | completed shotgun runs, one per specimen, 15 participants; SRR18585253 (participant 15, jejunum) excluded because mate 2 is unavailable | `LAB_ROOT/analyses/capscan_metag_analysis_2026-09-10/aggregate/selection.json`; `source_data/supplementary_methods_8_11/capscan_dna_runs.tsv` |
| 245 runs, 13 participants | native cross-omics profile: gut runs (272) of participants with DNA in all five groups (participants 5 and 7 drop; 257 runs), restricted to the protein and metabolite cross-omics specimen set (12 dropped) | `text_statistics/supplementary_methods_8_11/stats.py:stat_dna_runs` |
| 282 | 331-roster acquisitions whose specimen has one of the 301 runs (250 BioSamples, 15 participants, 4 stool acquisitions of participant 15) | `source_data/SData1/v28d_capscan_mg_pairs.tsv` |

347 − 282 = 65 = 36 (331 roster, specimen with assembly only) + 13 (331 roster, no
assembly) + 16 (participant-15 additions, not searched against specimen metagenomes; of
these 13 have a specimen with a completed run, 1 the excluded run, 2 none). Under the 347
roster, 295 acquisitions have a specimen with a completed run; the matched benchmark ran on
282 of the 331 roster. The 295-acquisition benchmark is `pipelines/capscan_mg_benchmark_295`.
Record: `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/DNA_COUNTS.md`.
