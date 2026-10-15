# mg295 phase 3: items 1–4 of PHASE2 "Not recomputed", at 295 CapScan acquisitions (7 Oct 2026)

No manuscript text, DOCX or release repository was edited. Each item has its own folder with ITEMn.md, code, sbatch, logs and JOBS.txt.
All jobs ran on p.hpcl92 with --nice=10000 (no hpcl8/93).

Status: items 1, 3 and 4 are COMPLETE. Item 2: the grouping at 295 is complete and its 282 control is byte-identical. The export, the Supp. Fig. 8b–d ablation, the full 282 control and the negative control are still queued or running (see below).

| claim (para) | 282 value | 295 value | producer | control result | job IDs |
|---|---|---|---|---|---|
| 1. Min-7 CapScan gain (398, SFig4a,b) | 45.0% (44.972) | **43.8%** (43.80) | M09 `run_task.py`/`task.sbatch` (byte-identical), M09 collector rows; `stats.py:stat_same_spectrum`, `panels.py:stat_ab_matched` | PASS: 1311 min-7 accepted sets identical (8,499/6,745); 2,890 matched rows = M09 table, 0 mismatches; producers return 44.9721893631811 / 84.01297879503241 exactly; negative controls caught | 6085669, 6085717, 6085718, 6085721, 6085785, 6085786 |
| 1. Min-7 CapScan identity recovery (398) | 84.0% (84.013) | 84.0% (84.03), unchanged | same | same | same |
| 1. PREDICT 14.5% / 85.4% | 14.5 / 85.4 | unchanged | same | same | same |
| 1. Min-7 search counts (126, 225) | 2,022 (2,017 standard) | 2,035 (2,030) | search inventory + 13 MG min-7 searches | receipts 13/13 PASS | same |
| 1. SFig4 legend CapScan n (225) | 282 | 295 | sd295 tables | — | — |
| 2. Matched-CapScan study dictionary / SFig8b–d (129, 235, 331) | no printed values; roster "282" | roster "295" once collected; FL/metaG median quantified groups 8,519.5/6,028.0 → 8,477.0/6,179.0 (not printed) | `simplified_grouping_results/CAPSCAN_two_arm_282_2026-09-11_immediate` `dispatch.py`/`study.py` (byte copy); `export_assignments.py` (path-only diff); `h09_evidence_ablation.py` (byte copy); `plot_corrected_ablation.py` transform | grouping: 282 SUMMARY.json and matched_protein_group_counts.tsv byte-identical. Full 282 control (export+ablation) and negative control PENDING | 6085760–6085762, 6085764–6085771 (tail submits 295 export/ablation/negctrl, IDs appended to item2_dictionary/JOBS.txt) |
| 3. Name collision, sequences (288) | 2,327,521 of 28,527,640 (8.2%) | **2,393,928 of 29,163,863 (8.2%)** | `analyses/sequence_accumulation_2026-09-29/worker.py` (sha 331ad8e4, byte copy) | PASS: 282 header_audit.tsv byte-identical (d685580b) and all other outputs identical; negative control (swapped assembly) stopped at worker gate | 6085728, 6085729, 6085730, 6086073 |
| 3. Name collision, shared names (288) | 2,414,558 of 2,414,668 | **2,494,126 of 2,494,236** (>99.9%) | same | same | same |
| 3. Multi-accession sequences (288) | 53 | 53, unchanged | same | same | same |
| 4. Matched pairs (326) | 578 | 591 | `Review_v16/Confidence/lineage/collect_v2.py` summary (verbatim copy) | PASS: 282 summary byte-identical to archived ordinary_summary.tsv (34bf8961) | 6085666, 6085691, 6085704, 6085731 |
| 4. PEP≤0.01 retention CapScan 5–50 (326) | 84.6% | **84.7%** (0.846912) | same | same | same |
| 4. PEP retention PREDICT | 85.7% | unchanged | same | same | same |
| 4. Cross-fit PEP 9–50 intestinal (326) | 28.3% → 36.3% (n=226) | **26.2% → 34.7%** (24.6–44.3), n=239 | `C1_pep_primary/crossfit/cf_pep.py` option (a) + `cf_analyze.py` (byte copies; analyze.py parameterises paths only) | PASS: cf_pep rerun on rows 0/296/335 = C1 JSONs; analyze control reproduced all three C1 tables byte for byte; negative control (13 altered rows) caught | same |
| 4. Cross-fit PEP PREDICT / CapScan stool | −3.5→1.7 / −5.5→−1.1 | unchanged | same | same | same |

## Item 6 (fold-cache pool, not rerun)
The 282 cross-fit rows (and the 282 part of the item-4 cross-fitted PEP pooling) use the 331-pooled iteration4 fold caches. The 13 new acquisitions use the 347-pooled public347 caches. The 295 numbers therefore mix the two pools, as in phase 2. Item 4 checked …1311 from both pools: baseline pooled sets and MG under PEP identical; FL/FULL under PEP differ by ordinary Sage per-fold variation (7,395 vs 7,401; 8,306 vs 8,292 at 5–50).

## Edits
`edits_v301_mg295_phase3.jsonl`: 8 edits (items 1, 3, 4), each anchored once in v30/visible_v30f.txt, no overlap with phase-2 anchors. Item 2's three roster edits (paras 129, 235, 331: 282 → 295) are written by `item2_dictionary/collect.sh` only if both its controls pass; append them afterwards:

```
bash item2_dictionary/collect.sh   # prints control results, 282 vs 295 SFig8b–d medians, writes item2_dictionary/item2_edits.jsonl
cat item2_dictionary/item2_edits.jsonl >> edits_v301_mg295_phase3.jsonl
```

## Not done / notes
- No figures rebuilt (SFig04 needs the 13 MG rows from item1 `per_acquisition_MG_new13.tsv`; SFig08 needs item 2 output).
- Release README `pipelines/metagenome_header_audit` and Source Data `supplementary_methods_1_5/header_audit.tsv` still state 250/282; the 295 table is `item3_collision/lab/review/v15_integration/analyses/run295/full/CAPSCAN/header_audit.tsv`.
- Item 2 FL arm follows the original producer: canonical `noclust_2026-08-07` CAPSCAN `5_sage/razor_hash-acc` searches (which contain the 13), not A08_347.
- Not printed in para 326: same-spectrum CapScan PEP gain at 5–50 goes from 61.7% to 59.9%.
- Item 3 specimen IDs for the 13 are `15|<workbook_Main>` (aliquot_main_id basis, as production uses for participant 15), giving 263 specimens.
- Rule breach (item 3 worker): one recursive grep under review/v15_integration on the login node, stopped before it finished with no output.
