# Discussion counts, CapScan 282 -> 295 (7 Oct 2026): COMPLETE, all checks PASS

Roster change only. `common.py`, `a_accumulation.py`, `e_candidates.py` are byte copies of `discussion_counts/`; the upstream producers they read
(`sequence_accumulation_2026-09-29/worker.py`, sha 331ad8e4…, and `candidates/worker.py`, sha e82ada22…) are byte copies run inside mirror trees
`trees/{c282,r295}/` (symlinks only for `confidence/` and the archived human Swiss-Prot FASTA). The roster is the manifest:
`trees/r295/.../manifest.json` = production manifest (578 rows unchanged) + 13 ACQ13 rows (`build_manifest295.py`; FL_public = public347/A08_347,
MG = mg295/db/_db + mg295/mg_full/out; specimen = `15|main:<workbook_Main>`). The only code difference: `common.py` line 8 `ACC` points to the
tree being run (`common_c282.py`, `common_r295.py`; `common.py` restored to the original afterwards). PREDICT FL/MG for the accumulation readout come
from the archived `full/PREDICT` (symlink; PREDICT rows of the manifest are identical).

| claim | 282 (archived) | 295 | control (282 rerun) | job IDs |
|---|---|---|---|---|
| CapScan MG union (printed "28.5 million") | 28,547,993 -> 28.5 million | 29,184,216 -> **29.2 million** | 28,547,993, byte-identical accumulation.tsv | 6087487 (c282), 6087490 (295), producers 6087492/6087493 |
| CapScan FL union (printed "0.42 million") | 416,703 -> 0.42 million | 424,792 -> 0.42 million (unchanged as printed) | 416,703 | same |
| CapScan specimens / acquisitions | 250 / 282 | 263 / 295 | 250 / 282 | same |
| Fold MG/FL (printed "69- to 84-fold") | CapScan 68.5 -> 69, PREDICT 83.6 -> 84 | CapScan 68.7 -> 69; range unchanged "69- to 84-fold" | 68.5 | same |
| E1 CapScan MG candidates, 9-50 (printed "3,323") | 3,322.75 | 3,317.5 -> **3,318** | 3,322.75 | 6087488 (c282), 6087491 (295), producers 6087492/6087493 |
| E1 CapScan FL candidates, 9-50 (printed "278") | 277.5 | 276.0 -> **276** | 277.5 | same |
| PREDICT (44.8 M / 0.54 M / 6,083 / 277) | 44,771,168 / 535,781 / 6,082.5 / 277 | identical | identical | — |

Controls (COLLECT.json, job 6087499; 23/23 checks true):
- Positive: 282 rerun reproduces the archived accumulation.tsv, per_specimen, per_acquisition, header_audit, candidates per_acquisition and cohort medians, and the discussion_counts A_* and E1_* outputs byte for byte.
- In the 295 run, all 578 pre-existing E1 rows and 564 accumulation acquisition-arm rows equal the 282 control; 13 / 26 new rows; worker COMPLETE 295 acq / 263 specimens / 15 participants.
- Negative (roster gate, job 6087486, MANIFEST295.json): dropping one of the 13, merging two specimens, or swapping two MG arms is rejected.
- Negative (producer, job 6087505): MG results of two new rows swapped with their expected peptide sets kept -> worker.py fails its set assertion (NEGATIVE_CONTROL_PASS).
- Note: a first producer negative control (job 6087489, `trees/neg_armswap`) swapped whole MG arms and the worker PASSED; the worker cannot detect a consistent arm-to-acquisition mis-assignment. That is caught by the roster gate only.

Memory: archived MaxRSS x1.05 (accum 18.6 GB -> 20G; 295 scaled -> 21G, used 19.2 GB; producers 5.6 -> 6G/7G).

Edits: `edits_discussion_295.jsonl` (2 pairs, anchored to v301/vis.txt wording; the text says "candidate peptides per spectrum").
Figures: the accumulation curves (`sequence_accumulation_2026-09-29/figures/sequence_accumulation_*`) and candidate figure (`candidates/figures/candidate_competition_and_yield.*`) are built from the 282 outputs; neither is a release figure unit (no accumulation/candidate table in release_prep/fastalake-paper_v28 source_data). Not rebuilt.
