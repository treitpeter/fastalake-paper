# Item 4: PEP-gate summaries (para 326) at 295 CapScan acquisitions (7 Oct 2026)

Status: COMPLETE. Edits are in `item4_edits.jsonl` (3 edits, each anchored once in para 326, with no overlap with phase-2 MG295-23).

## Producers (copied byte for byte into `c1_copy/`, with sha256 in `COPY_SHA256.txt`)
- **Same-spectrum PEP retention ("84.6%").** This is the `ordinary_summary.tsv` step of `Review_v16/Confidence/lineage/collect_v2.py`, run by `conf_summary.py` with the statements copied verbatim.
  - Input at 295: the archived 578×36 `ordinary_all_pairs.tsv` plus `phase2/per_acq13/confidence_rows.tsv`. Those 13×36 rows come from the frozen `worker_v2` reader; their phase-2 control on 1311 matched 36 of 36.
- **Cross-fitted PEP gains ("−1.1%", "36.3%").** The source is `coauthor_comment_check_v25_2026-10-05/compute/C1_pep_primary/crossfit/cf_pep.py`, option (a): pooled q≤0.01 and Sage per-fold PSM PEP≤0.01.
  - Summaries are produced by `cf_analyze.py`.
  - `run_pep.py` imports `cf_pep` unchanged and only redirects `cf_pep.IT4` to `mg295/crossfit` for the 13.
  - `analyze.py` is `cf_analyze.py` with inputs and outputs parameterised; the diff is in `analyze.diff`, and the statistics are unchanged.
  - The 282 rows are the archived C1 `per_acq/*.json` (578 records). The 13 rows are in `per_acq13/578..590.json`.

## Results

| claim | 282 value | 295 value | producer file | control result | job IDs |
|---|---|---|---|---|---|
| matched pairs | 578 | 591 | conf_summary.py | n/a | 6085691 |
| PEP≤0.01 FL retention 5–50, CapScan | 84.6% (0.846194) | **84.7%** (0.846912) | collect_v2.py → conf_summary.py | Control: the summary is byte-identical to the archived `ordinary_summary.tsv` (sha 34bf8961…). Negative control: 13 planted `removed_fl` changes were detected. A first version that planted only a single row did not move the median, so it was replaced (kept as `negctrl_v1_single_row_not_detected/`). | 6085691, 6085731 (collect) |
| same, PREDICT | 85.7% | 85.7% (unchanged) | same | same | same |
| cross-fit 9–50, PREDICT, base → PEP | −3.5% → 1.7% | unchanged (−3.49 → 1.71) | cf_pep.py / cf_analyze.py | Control: rerunning `cf_pep` on rows 0, 296 and 335 gave records equal to the C1 `per_acq` JSONs. `analyze.py control` gave `CF_PEP_PER_ACQUISITION`, `SUMMARY` and `SUMMARY_BY_REGION` byte-identical to C1. | 6085666 (pilot and control), 6085691 |
| cross-fit 9–50, CapScan stool | −5.5% → −1.1% | unchanged (−5.51 → −1.12; CI −4.41 to 2.28) | same | same | same |
| cross-fit 9–50, intestinal fluid | 28.3% → 36.3% (n=226) | **26.2% (17.0–34.3) → 34.7% (24.6–44.3)**, n=239 | same | The base 26.15 (16.97–34.30) equals the phase-2 Fig. 4 value of 26.2% (17.0–34.3). For all 13, the baseline pooled sets equal the mg295 frozen POOLED_ACCEPTED sets. | 6085666, 6085704 (array 1–12), 6085731 (collect) |

Other 295 values (not printed in para 326), from `out295/CF_PEP_SUMMARY.tsv`:
- CapScan all at 9–50: 17.61 → 25.11.
- Intestinal fluid at 5–50: 38.34 → 47.20.
- Intestinal acquisitions with more FL at 9–50 under PEP: 223 of 239 (was 214 of 226).
- The same-spectrum CapScan PEP gain at 5–50 is 59.9% (was 61.7%). It is not printed in para 326.

Memory: the pilot (6085666) peaked at 696 MB, so the array requested 800M (696 MB × 1.05). Partition p.hpcl92 with `--nice=10000`.

## Fold-cache pool (item 6)
- **282 rows:** the C1 PEP recount of iteration4 cross-fit outputs, which are the **331-pooled** fold caches.
- **13 rows:** the mg295 cross-fit outputs, which are the **347-pooled** public347 caches.
- **Diagnostic (`control/pool1311.json`):** I compared 1311 run from both pools.
  - Baseline pooled sets are identical in all arms.
  - MG under PEP is identical.
  - Under PEP, FL changes from 7,395 to 7,401 and FULL from 8,306 to 8,292 at 5–50. Per-fold Sage PEP varies at this level between reruns.
  - The 295 numbers therefore mix the two pools, as phase 2 does.

## Outputs
- `n295/ordinary_summary.tsv` (sha 5108ef1a…) and `n295/ordinary_all_pairs_recomputed.tsv`.
- `out295/CF_PEP_{PER_ACQUISITION,SUMMARY,SUMMARY_BY_REGION}.tsv`.
- `control/` and `logs/`.
