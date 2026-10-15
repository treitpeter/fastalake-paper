# Discussion counts (6 October 2026)

All from existing outputs; no new searches. Every computation ran on Slurm (p.hpcl92/8/93, nice 10000): jobs 6073023 (A), 6073031 (C), dcB/dcD/dcE1/dcE2/dcBx (logs/ holds stdout/stderr). Scripts in this folder; `common.py` holds the peptide rule.

**Peptide rule (all tasks):** Sage rank 1, label 1 (target), peptide_q <= 0.01, modifications stripped, I->L, distinct; 9-50 primary, 5-50 alongside. This is `review/v15_integration/analyses/confidence/peptide.py::canonical`.

## A. Sequence accumulation
Source: `review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/{PREDICT,CAPSCAN}/accumulation.tsv` (592 and 500 rows read), `figures/SUMMARY.json`. Sequence = exact amino-acid string (SHA-256, no I/L collapse). Output: `A_accumulation_readout.tsv`, `A_RESULT.json`.

| cohort | arm | specimens | after 50: fixed / order-averaged | after 100: fixed / order-averaged | final union |
|---|---|---|---|---|---|
| PREDICT | FL_public | 296 | 180,560 / 193,638 | 297,101 / 297,761 | 535,781 |
| PREDICT | MG | 296 | 8,420,630 / 9,199,502 | 17,254,979 / 17,142,534 | 44,771,168 |
| CapScan | FL_public | 250 (282 acq.) | 167,676 / 185,842 | 261,212 / 270,726 | 416,703 |
| CapScan | MG | 250 (282 acq.) | 5,924,945 / 6,583,022 | 12,910,772 / 12,365,381 | 28,547,993 |

Combined PREDICT+CapScan (SUMMARY.json): FL 758,358 (194,126 shared); MG 71,726,951 (1,592,210 shared).
Growth (order-averaged new sequences per added specimen, first 50 vs last 50): PREDICT FL 3,873 -> 952, MG 183,990 -> 134,236; CapScan FL 3,717 -> 783, MG 131,660 -> 105,596.

Per-specimen peptide IDs (union of acquisitions within a specimen), from the `expected_sets` files worker.py verified against (2,312 files read; 5-50 counts equal `per_specimen.tsv` in 1,092/1,092 rows). Output `A_per_specimen_medians.tsv`.

| cohort | FL 5-50 | FL 9-50 | MG 5-50 | MG 9-50 |
|---|---|---|---|---|
| PREDICT (296) | 27,254.5 | 21,483.5 | 22,284.5 | 18,847 |
| CapScan (250) | 30,151 | 24,447.5 | 17,339 | 14,572 |

CapScan 347: the matched curve must stay at 282 acquisitions (MG exists only for 282). FL alone, from the A08_347 razor FASTAs (347 read): the 282 reproduce 416,703 exactly; all 347 give 449,065 (+32,362 from the 16 participant-15 acquisitions, +7.8%). Median 9-50 per acquisition (length9_50.txt): 22,237 (282), 21,744 (347), 18,972 (new 16). The specimen-level curve for 347 is not built: the 16 have no specimen IDs in the accumulation manifest. A paired FL/MG plot is unchanged; an FL-only statement over 347 would say 449,065.

## B. smORFs (GMSC)
**Source rule** (first header token; reservoir kept header is the highest-priority source, UHGG > GMGC > GMSC > HUMAN_SP, `canonical/publication_data/cohorts/HSTOOL/0_lake/HSTOOL_lake_stats.txt`): `MGYG*` = UHGG, `GMGC10.*` = GMGC, `GMSC10.*` = GMSC, `sp|*` = Swiss-Prot human. 0 headers fell outside these in all 711 searched FASTAs. Because of the priority, a GMSC-tagged sequence is never identical to a UHGG or GMGC sequence. At most 89 HUMAN_SP sequences duplicate a higher-priority sequence across the whole reservoir.

**Reservoir** (233,497,702 unique sequences; PREDICT_lake.fasta has the same stats): unique sequences by kept source are UHGG 95,303,634 (40.8%), GMSC 91,315,439 (39.1%), GMGC 46,858,287 (20.1%), Swiss-Prot human 20,342 (0.009%). Input header records were UHGG 95,303,635, GMSC 94,490,666, GMGC 52,544,908, HUMAN_SP 20,431 (242,359,640 total). 7,767,598 sequences carry more than one source header. The same per-source record counts appear in `final/public347/A08_347/evidence/COMPLETE.json`.

**Searched FASTAs**: HSTOOL 58 = `arms/noclust_2026-08-07` parsimony.fasta (roster `simplified_grouping/.../hstool_samples.txt`). PREDICT 306 = `arms/abl_nometag_noclust`, with results from the accumulation manifest for the 296. CapScan 347 = `final/public347/A08_347`. The searched FASTA was read from each search's results.json. GMSC-only peptide = every protein in the searched FASTA containing it (I/L-equivalent exact substring) is GMSC-tagged. "Absent from MG" = not an exact I/L-equivalent substring of any sequence in the matched `metag_hashed_plus_human` FASTA actually searched. Output: `B_per_acquisition.tsv`, `B_summary.tsv`, `B_extra.txt`, peptide lists `B_gmsc_only/`. 0 accepted peptides were unmapped. The 9-50 set equals the reference set for every PREDICT/HSTOOL/A08_347 acquisition. For 4 of 282 CapScan acquisitions the A08_347 set differs from the production A08 set (…1254, …1251, …1368, …1815). A08_347 is used throughout.

| | HSTOOL 58 | PREDICT 306 (296 matched) | CapScan 347 (282 matched) |
|---|---|---|---|
| GMSC proteins, pooled / median per acq. | 24,600 of 709,676 (3.47%) / 426 (3.40%) | 89,860 of 2,370,126 (3.79%) / 316 (3.78%) | 96,427 of 2,759,980 (3.49%) / 287 (3.60%) |
| GMSC-only 9-50 peptides, median per acq. (share) | 420 (1.52%) | 312 (1.42%) | 279 (1.28%) |
| GMSC-only 5-50, median (share) | 464.5 (1.40%) | 363.5 (1.30%) | 322 (1.17%) |
| distinct GMSC-only 9-50 peptides | 4,784 | 19,568 | 13,372 |
| GMSC-only 9-50 absent from matched MG: summed / share | n/a (no MG) | 26,808 of 89,504 (30.0%); per-acq. median 87 (27.5%) | 37,898 of 76,952 (49.2%); per-acq. median 127 (46.8%) |
| same, 5-50 | n/a | 28,855 of 103,387 (27.9%); median 94.5 (25.6%) | 41,035 of 88,058 (46.6%); median 136 (43.4%) |
| distinct 9-50 GMSC-only peptides absent from every matched MG in which they were identified | n/a | 3,685 | 4,080 |

As a share of all accepted 9-50 peptides per acquisition, GMSC-only and absent from MG is a median of 0.39% for PREDICT and 0.60% for CapScan.

## C. Parsimony traceability
Code (canonical @ 66f1e6ae):
- `publication/chain/03_persample.sh:101-108`: the extractor writes the per-acquisition candidate set `3_persample/<S>/sample_lake.fasta` and `sample_hits.tsv`. The hits file holds peptide, protein_id and position for every **de novo** peptide x every candidate protein containing it, written at `rust/fasta_extractor_v2/src/main.rs:1221-1241`. No chain stage deletes either file; the only `rm` is the temp peptide dir at line 111. Both are present on disk for this acquisition.
- `publication/chain/04_parsimony.sh:110,126-127`: parsimony reads sample_lake.fasta and the de novo CSV. `rust/parsimony_engine/src/main.rs:1511-1512` builds pep2prot in memory. Only the selected FASTA and stats JSON are written (`:1594-1632`). The peptide-to-all-proteins graph is not saved.
- After search, Sage's `proteins` column lists only proteins in the searched (parsimony) FASTA.
- So for an accepted peptide, alternative homologs are recoverable after the search by re-matching it against the retained candidate FASTA (all 24,430 accepted 9-50 peptides occur in it) or the reservoir. They are not stored as a ready-made accepted-peptide map. `sample_hits.tsv` covers only accepted peptides that were also de novo evidence peptides: 12,466 of 24,430 (51.0%).

Worked acquisition 20181231_QX7_PhGe_LC95_SA_WorkflowReplicates_2: 126,152 PSM rows read, 24,430 accepted 9-50 (29,663 at 5-50). Records scanned: 12,606 searched, 6,671,447 candidate, 26,881,038 evidence, 233,497,702 reservoir. Unit = distinct sequences containing the peptide (I/L-equivalent). Output `C_per_peptide.tsv`, `C_RESULT.json`.

| 9-50 | median | IQR | 1 only | >=10 | >=100 | >=1000 | max |
|---|---|---|---|---|---|---|---|
| searched FASTA | 1 | 1-3 | 13,118 | 1,592 | 28 | 0 | 217 |
| candidate set | 174 | 51-493 | 1,069 | 22,030 | 15,297 | 3,069 | 51,467 |
| reservoir | 201 | 66-556 | 771 | 22,717 | 16,296 | 3,484 | 51,467 |

At 5-50 the medians are 1 / 195 / 253 (reservoir max 345,364). Parsimony assigns one razor protein, while the searched FASTA itself still holds a median of 1 (mean 3.4) containing proteins.

## D. Clustering cost (CAMPI, 15 acquisitions)
The printed 0.9/2.2/2.5/4.1/7.2/8.8% are **5-50**: `final/clustering_check/CHECK.md` recount rule "5-50 residues", and `summary.tsv` = 0.871/2.157/2.457/4.149/7.226/8.774. Recomputed here from the same 105 Sage outputs (`review/v15_integration/analyses/clustering_90_95_2026-09-29` + prior_search; 7,725,921 PSM rows read). Every 5-50 count equals `per_acquisition_recount.tsv`. Output `D_summary.tsv`.

| identity | 99% | 97% | 95% | 90% | 75% | 50% |
|---|---|---|---|---|---|---|
| median reduction 5-50 | 0.87 | 2.16 | 2.46 | 4.15 | 7.23 | 8.77 |
| median reduction 9-50 | 0.86 | 2.07 | 2.35 | 3.70 | 6.90 | 8.34 |

At 9-50, all 15 acquisitions are below baseline at every threshold.

## E. Candidate competition (Sage 0.14.6 `scored_candidates`)
**E1** (`analyses/sequence_accumulation_2026-09-29/candidates/`): recomputed from all 578 manifest searches (150,322,532 PSM rows read). 0 of 578 acquisition medians differ from `full/per_acquisition.tsv`. **Window:** no peptide-length filter on spectra. The population is every rank-one spectrum ID reported by both arms, target or decoy, before q filtering. The listed peptide gain uses 5-50. **Design:** same-spectrum, frozen databases, not cross-fitted. All paired spectra already have a rank-one length of 5-50 (Sage min_len 5).

| cohort | FL median / MG median (ratio) all paired | 9-50 (both arms' rank-one 9-50) | paired spectra (9-50) |
|---|---|---|---|
| PREDICT 296 | 228 / 5,128 (21.2x) | 277 / 6,082.5 (21.2x) | 41,669,631 (31,355,421) |
| CapScan 282 | 241.5 / 2,924.5 (14.3x) | 277.5 / 3,322.75 (14.2x) | 27,105,787 (20,327,993) |

These are cohort medians of per-acquisition medians. Output `E1_per_acquisition.tsv`, `E1_summary.tsv`.

**E2** A97 (`sample_specific_merge_all_raw/2026-09-11`): the PREDICT (306) and CapScan (331) joint union searches have no Sage output. Both the original and the 700 GiB retries ended OUT_OF_MEMORY in LFQ tracing (`STATUS_CHECK_2026-09-14.md`, jobs 5907491/5911555, 5907492/5911554). There is therefore no scored_candidates for the FL study union in PREDICT or CapScan. Union sizes are FL 555,987 (PREDICT 306), 452,255 (CapScan 331) and 460,977 (CapScan 347, A97 README). The MG unions, as database size only, are 44,771,168 (PREDICT 296) and 28,547,993 (CapScan 282).

Only HSTOOL 58 has an FL study-union search: the A06 one-acquisition searches, each acquisition against its own FASTA and against the 153,267-sequence union, plus the 58-file joint search (9,001,927 rows). On 7,751,917 paired rank-one spectra the median is 237 for per-sample FL vs 1,745.5 for the study union (7.05x). At 9-50 it is 278 vs 2,063.5 (7.16x). The joint search gives identical per-spectrum values to the A06 union searches, which share the same database. Output `E2_hstool_per_acquisition.tsv`, `E2_RESULT.json`.

## CapScan matched roster 282 -> 295 (7 October 2026; v33)

The v33 matched CapScan roster adds the 13 participant-15 acquisitions of `pipelines/capscan_mg_benchmark_295` (295 acquisitions, 263 specimens). A and E1 were re-run for that roster with byte copies of `common.py`, `a_accumulation.py` and `e_candidates.py`; the only code difference is `common.py` line 8 (`ACC`, the tree being read; `run295/common_r295.diff`, `run295/common_c282.diff`). The upstream producers (sequence_accumulation `worker.py`, sha 331ad8e4..., and `candidates/worker.py`, sha e82ada22...) ran unchanged in mirror trees on the production manifest (578 rows unchanged) plus 13 rows (`run295/build_manifest295.py`; FL_public = A08_347, MG = mg295 databases and same-spectrum searches; specimen = `15|main:<workbook_Main>`). Scripts and launchers: `run295/`.

| claim | 282 (archived) | 295 |
|---|---|---|
| CapScan MG union (printed "28.5 million" -> "29.2 million") | 28,547,993 | 29,184,216 |
| CapScan FL union ("0.42 million") | 416,703 | 424,792 (unchanged as printed) |
| CapScan specimens / acquisitions | 250 / 282 | 263 / 295 |
| fold MG/FL ("69- to 84-fold") | CapScan 68.5, PREDICT 83.6 | CapScan 68.7; range unchanged |
| E1 CapScan MG candidates, 9-50 median | 3,322.75 (printed 3,323) | 3,317.5 (printed 3,318) |
| E1 CapScan FL candidates, 9-50 median | 277.5 (printed 278) | 276.0 (printed 276) |
| PREDICT (44.8 M / 0.54 M / 6,083 / 277) | unchanged | unchanged |

Controls (COLLECT.json, 23/23 true): the 282 re-run reproduces the archived accumulation, per_specimen, per_acquisition, header_audit, candidate per_acquisition and medians and the A_* and E1_* outputs byte for byte; in the 295 run all 578 pre-existing E1 rows and 564 accumulation acquisition-arm rows equal the 282 control. Negative controls: the roster gate rejects a dropped acquisition, merged specimens and swapped MG arms (job 6087486); swapping the MG results of two new rows makes worker.py fail its set assertion (job 6087505). A first producer negative control that swapped whole MG arms consistently passed the worker (job 6087489): the worker cannot detect a consistent arm-to-acquisition mis-assignment, which only the roster gate catches.

Sections B (GMSC) and C, D were not re-run: B uses the 282 matched CapScan acquisitions for its "absent from matched MG" rows, and C, D do not involve CapScan.
