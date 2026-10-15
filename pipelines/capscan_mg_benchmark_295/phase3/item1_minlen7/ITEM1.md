# Item 1: min-length-7 matched benchmark (Supp. Fig. 4a,b; Supp. Methods 7 para 398) at 295 CapScan acquisitions

Status: COMPLETE (7 Oct 2026). All jobs on p.hpcl92, --nice=10000.

## Method (same producer as the 282 values)
- MG min_len-7 searches: `code/run_task.py` and `code/task.sbatch` are byte copies of M09 `run_2026-09-24_minlen7/code/` (sha 500793ea.../bba43f16...; same Sage 0.14.6 binary sha d3820543..., receipt checks). `code/prepare13.py` applies the M09 transformation (copy source results.json, drop version/output_paths, min_len 5->7, output_directory) to `mg295/mg_full/out/<s>/results.json`.
- FL_public min_len-7 for the 13: existing `public347/minlen7_public347/searches/` (unchanged).
- Collector `code/collect13.py`: read() and the matched-row formula are verbatim M09 `collect_peptides.py`/`summarize_peptides.py`; statistics are taken by the release producers themselves, `text_statistics/supplementary_methods_6_7/stats.py:stat_same_spectrum` and `figures/SFig04/panels.py:stat_ab_matched` / `stat_gain_5_8`, on 295-row copies of the two Source Data tables (`sd295/`, 282 committed rows unchanged + 13 appended).

## Results
| claim | 282 value | 295 value | producer | control result | job IDs |
|---|---|---|---|---|---|
| min-len-7 same-spectrum gain, CapScan (para 398 "45.0% (CapScan)") | 44.97% | **43.80% -> "43.8%"** | stats.py:stat_same_spectrum (= panels.py:stat_ab_matched fresh7) | 282 reproduced exactly (44.9721893631811) from raw searches and from committed table | 6085717, 6085785, 6085786 |
| min-len-7 identity recovery, CapScan ("84.0%") | 84.01% | 84.03% -> "84.0%" (unchanged) | same | 84.01297879503241 reproduced exactly | same |
| min-len-7 PREDICT gain / recovery ("14.5%", "85.4%") | 14.46 / 85.39 | unchanged (PREDICT rows untouched) | same | reproduced exactly | 6085721 |
| SFig04 5-8 residue gain CapScan (panels.stat_gain_5_8; text "96%") | 95.56% | 93.75% -> "94%" | panels.py:stat_gain_5_8 | 95.56216698309652 reproduced | consistent with phase2 edit MG295 "96%->94%" (para 398) |
| SFig04 a, CapScan n (para 225 legend) | 282 | 295 | panels.py:stat_ab_matched | | |
| min-len-7 searches repeated (paras 126, 225) | 2,022 (2,017 standard) | 2,035 (2,030 standard) | panels.py:stat_search_count (inventory + 13 MG receipts) | 13/13 receipts PASS | 6085785 |
| not in v30 text, for Source Data only: CapScan nested9 n_more_fl | 281/282 | 294/295 | panels.stat_ab_matched | | |
| not in text: panel c CapScan MG median change, re-search at >=7 | +0.90% | +0.86% (filtered: -0.26 -> -0.26) | M09 wide formula | | |

## Controls
- Config control (job 6085669, `logs/control282_configs.json`): M09 transformation reproduces all 282 archived CapScan MG configs (0 differ); the 1311 control config built from mg295/mg_full differs from M09's only in the fasta path (sha-identical file) and output_directory. PASS.
- Search control (task 0, job 6085717_0): 1311 re-searched at min_len 7 with the copied code; accepted canonical sets at >=5/>=7/>=9 identical to the archived M09 search (8,499 / 8,499 / 6,745). PASS.
- Collector control (job 6085721, `CONTROL_282.json`): all 2,890 FL_public/MG matched rows (296 PREDICT + 282 CapScan x 5 populations) recomputed from the archived results.sage.tsv equal M09 matched_comparisons.tsv (0 mismatches); producers give 44.9721893631811 / 84.01297879503241 / 14.455684375777889 / 85.39185872471847 from both committed and rebuilt tables. PASS.
- Gates for the 13: FL_public rows equal public347 per_acquisition_new16.tsv; MG min-5 counts equal mg_full COUNTS.json (13/13). PASS.
- Negative controls: task with corrupted source sha refused by run_task.py (job 6085718, `negative/run.log`, AssertionError, exit 1); a single perturbed reference count is detected by the row comparison (CONTROL_282.json). PASS.

## Memory
Pilot: largest DB (1803, 60.8 MB) MaxRSS 19.03 GB -> 20G for the remaining 12 (max observed 18.8 GB).

## Outputs
`RESULTS_295.json`, `CONTROL_282.json`, `matched_comparisons_new13.tsv`, `per_acquisition_MG_new13.tsv`, `sd295/{supplementary_methods_6_7/minlen7_benchmark_per_acquisition.tsv, SFig04/benchmark_matched_per_acquisition.tsv}`, `item1_edits.jsonl` (4 edits, all found once), `searches/`, `receipts/`, `logs/`, `JOBS.txt`.
Source Data note: SFig04 `search_inventory.tsv` and `per_acquisition_length_counts*.tsv` would need the 13 MG rows (`per_acquisition_MG_new13.tsv`, receipts) for the figure rebuild; the figure itself was not rendered.
