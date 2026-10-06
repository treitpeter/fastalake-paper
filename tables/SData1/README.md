# Supplementary Data 1

Supplementary Data 1 (`Sample_rosters.xlsx`) is the sample index of the paper. For every
acquisition of the four datasets it lists the analysis rosters the acquisition belongs to,
its recorded specimen and participant links, and the documented exclusions, so that every
n in the paper can be traced to named acquisitions. The notebook rebuilds the workbook as
published with v28d (CapScan on 331 acquisitions, 710 acquisitions in all), checks it
against the published file, builds the v28f workbook with CapScan on the
347-acquisition roster (726 acquisitions in all; five CapScan analysis rosters still on 331), and
builds the v30 workbook, in which every CapScan analysis is on the 347-acquisition roster. The v30
workbook is the Supplementary Data 1 of the v30 manuscript.

Run: `python paper.py build SData1`. Outputs go to `build/SData1/`:
`Sample_rosters_v28d.xlsx`, `Sample_rosters_v28f.xlsx`, `Sample_rosters_v30.xlsx`,
`ANALYSIS_INDEX_v28d.tsv`, `ANALYSIS_INDEX_v28f.tsv`, `ANALYSIS_INDEX_v30.tsv`,
`CapScan_specimen_number_discrepancies_v28f.tsv`. The 347-roster inputs of the v30 build were copied
into `source_data/SData1/v30_*.tsv` once by `prepare_v30_inputs.py` (Slurm job 6075320).

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| workbook v28d | 34 sheets as published | "rosters ... for all 710 acquisitions" (489) | 710 acquisitions | none (counts) | Table v28d | `source_data/SData1/v28d_*.tsv` | sample_rosters_workbook, source_data_package |
| workbook v28f | 34 sheets, CapScan on 347 | "347 acquisitions from 308 specimens and fifteen participants" (501) | 726 acquisitions; CapScan 347 acquisitions, 308 specimens, 15 participants | none (counts) | Table v28f | `v28d_*.tsv`, `CAPSCAN_participant15_join.tsv`, `CAPSCAN_eligibility_audit_347.tsv`, `source_data/Fig5/design_347.tsv` | capscan_participant_roster, capscan_crossomics |
| workbook v30 | 34 sheets, every CapScan analysis on 347 | "analysis-specific acquisition rosters ... for all 726 acquisitions" (v30 text paragraph 259, after edit V30S-1) | 726 acquisitions (CAMPI 15, healthy stool 58, PREDICT 306, CapScan 347 from 308 specimens and 15 participants) | none (counts) | Table v30 | `v28d_*.tsv`, `CAPSCAN_*_347.tsv`, `v30_*.tsv`, `source_data/Fig5/design_347.tsv`, `source_data/supplementary_methods_8_11/capscan_dna_runs.tsv` | capscan_participant_roster, capscan_crossomics, terminal_histidine_inference |
| specimen-number census | acquisitions whose file name and instrument record differ | "Sixteen acquisitions, fourteen of them in this roster" (503) | 368 RAW files audited; 16 discordant; 14 in the 347 roster | none (counts) | Table v28f | `CAPSCAN_raw_header_audit.tsv`, `CapScan_specimen_number_discrepancies.tsv` | capscan_raw_header_audit |
| roster numbers | cohort and analysis sizes | lines 489, 501, 503 | see NUMBERS.tsv | none (counts) | Numbers | as above | as above |

No statistical test is computed in this unit; `STATS.tsv` has the header only.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `ROSTER_V28D` | CAMPI 15, healthy stool 58, PREDICT 306, CapScan 331 | rosters of the published workbook (v28d lines 489, 501) | none |
| `CAPSCAN_ROSTER_V28F` | 347 | participant-15 join fix adds 16 intestinal-fluid acquisitions (v28d line 501, 503) | 331 (v28d workbook, rebuilt and checked here) |
| `CAPSULE_REGIONS` | duodenum, jejunum, ileum, colon | regions 1-4 (v28d line 201) | none |
| `P15_ALIQUOT` | pellet | aliquot type set for participant-15 fluid acquisitions by the producer (make_design347.py) | none |
| `WITHHELD_COLUMNS` | PREDICT visit codes (august_visit_id, later_visit_number, recorded_visit_number) and recorded group code (august_group) | clinical covariates may not enter the repository (CONVENTIONS.md section 8) | complete workbook with `FASTALAKE_CONTROLLED_DIR` |
| `WITHHELD_ROWS` | PREDICT_ANTIBIOTIC42 rows of PREDICT membership (306) and PREDICT patients (93) | membership reveals antibiotic exposure per patient | aggregate counts stay in the Analyses sheet (155 acquisitions, 153 specimens, 42 patients) |
| `V28F_REBUILT` | CAPSCAN_SEARCH_ALL, CAPSCAN_QUANT_NONEMPTY, CAPSCAN_REGIONAL_DESCRIPTIVE15, CAPSCAN_REGIONAL_PAIRED14 to PAIRED15 | their inclusion rules use only the roster and the nonempty flag, and regenerating them on 331 reproduces the published rows exactly | the other CapScan rosters stay on 331 (`V28F_NOT_REINDEXED_IDS`, OPEN_ITEMS.md) |
| `V28F_NOT_REINDEXED_IDS` | functional non-DNA, native, functional DNA, sensitivity and histidine rosters | their membership needs specimen-level and DNA eligibility of the 347 reruns, not indexed yet | OPEN_ITEMS.md |
| `V30_REINDEX` | NONDNA14 -> CAPSCAN_FUNCTIONAL_NONDNA15; NATIVE12 -> CAPSCAN_NATIVE13; FUNCTIONAL_DNA12 -> CAPSCAN_FUNCTIONAL_DNA13; FUNCTIONAL_NONDNA12_SENS -> CAPSCAN_FUNCTIONAL_NONDNA13_SENS | the 347 cross-omics run (capscan_crossomics, Slurm 6063849) has protein and metabolite profiles in all five groups for 15 participants and DNA profiles for 13 (participants 5 and 7 lack colon and duodenum DNA) | the 331 rosters (v28d and v28f workbooks) |
| valid cross-omics specimen | positive microbial-OG acquisition (`CAPSCAN_eligibility_audit_347.tsv`) and a metabolite row (`v30_capscan_metabolite_specimens.tsv`) | rule of build_rosters.py (metspec); recomputed set must equal the 347 cross-omics specimen roster (303 specimens) | none |
| DNA runs of the DNA analyses | gut/stool runs of the 301 whose specimen is valid and whose participant is DNA-complete (245) | rule of build_rosters.py; must equal the native_profile runs of `capscan_dna_runs.tsv` | none |
| `V30_H_QUANT_COLUMN` | in_quantitative_nonempty_roster (was in_quantitative330_roster) | on 347 the quantitative roster has 346 acquisitions | none |
| histidine sheets | finite n and finite participant sets of the fixed-347 A28 run (Slurm 6064031), summarised with the build_histidine_scope.py code | same summary rule as v28d | 331 run (v28d, v28f) |

## What changes from v28f to v30

Receipt `checks/SData1_v30.json` (`changes_from_v28f`, `counts`, controls):

- CAPSCAN membership: the four cross-omics analyses re-indexed. NONDNA15: 342 acquisitions, 303
  specimens, 15 participants, 75 cells. NATIVE13 and FUNCTIONAL_DNA13: 310 acquisitions, 277
  specimens, 13 participants, 65 cells, 245 DNA runs and 245 DNA specimens. NONDNA13_SENS: as
  NATIVE13 without DNA runs.
- CAPSCAN H acquisitions: 331 to 347 rows (16 participant-15 rows added); quantitative-roster column
  renamed and evaluated on 347. CAPSCAN H finite n and CAPSCAN H participants: from the fixed-347
  A28 run (fluid endpoints now up to 15 finite participants; 44,010 saved tests, as before).
- CAPSCAN exclusions and MS nonmembership: regenerated for the four analyses with the published
  reason texts; the participant-15 exclusions are gone (participant 15 now has all five groups).
- CAPSCAN DNA eligibility: 210 to 225 rows (participant 15 added); column in_complete12 renamed
  in_complete13. The 210 rows of participants 1 to 14 are unchanged.
- Analyses: the five rows updated (ids, counts, rules, provenance); the v28f "not yet indexed" note removed.
- Sheet directory: counts and the CAPSCAN H acquisitions description updated.
- Positive control: the membership generator reproduces the published 331 rows of the four
  cross-omics analyses from their published members (row construction and order).

## What changes from v28d to v28f

All in CapScan sheets and the computed sheets (receipt `checks/SData1_v28f.json`,
`changes_from_v28d`):

- CAPSCAN acquisitions and CAPSCAN ID crosswalk: 331 to 347 rows (16 participant-15
  acquisitions appended; the 331 published rows unchanged).
- CAPSCAN membership: the four re-indexed analyses regenerated on 347 (CAPSCAN_SEARCH_ALL
  347 acquisitions, 308 specimens, 15 participants, 75 cells; CAPSCAN_QUANT_NONEMPTY,
  CAPSCAN_REGIONAL_DESCRIPTIVE15 and CAPSCAN_REGIONAL_PAIRED15 346, 307, 15, 75).
- CAPSCAN MS nonmembership: regenerated for the four; 16 participant-15 rows added to
  CAPSCAN_PUBLIC_MG_282 (no specimen-matched metagenome database for them).
- CAPSCAN exclusions: the participant-15 row of CAPSCAN_REGIONAL_PAIRED14 removed.
- CAPSCAN regional pairs: 56 to 60 rows (fifteen participants).
- Analyses: the four rows updated; five rows still on 331 carry a note in `limitations`.
- Cohorts: CapScan 347 acquisitions, 308 specimens, 15 participants; design text "All
  fifteen participants contribute capsule and stool".
- Sheet directory: counts and three CapScan descriptions updated; the two rows naming
  sheets that are not in the published workbook (Sources, Author confirmations) dropped.

## Known limitations

- The public build (no `FASTALAKE_CONTROLLED_DIR`) lacks the withheld PREDICT cells; its
  Sheet directory therefore counts 2,142 PREDICT membership rows and 372 PREDICT patient
  rows instead of 2,448 and 465, and the Analyses row PREDICT_ANTIBIOTIC42 points to rows
  that are absent.
- The workbook xlsx bytes differ between runs (openpyxl stores the save time); the
  comparison is on cell values (sheet digests).
- v28f only: five CapScan analysis rosters stay on 331; v30 re-indexes them.
- `diff_builds` labels its counts rows_v28d and rows_v28f for any pair of builds (from and to).
- The Read me sheet text is reproduced as published, including its review-checkpoint
  wording (OPEN_ITEMS.md).
