# Item 2: matched-CapScan study dictionary and Supp. Fig. 8b–d ablation at 295

## What the 282 values came from (producer chain, all re-run here unchanged)

1. **Study dictionaries (FL and metaG):** `simplified_grouping_results/CAPSCAN_two_arm_282_2026-09-11_immediate`.
   - `code/dispatch.py` (sha256 e068992d…) runs `code/study.py` (5f69c457…), `group_searches(q=0.01, config results.json)`.
   - FL = `canonical/publication_data/arms/noclust_2026-08-07/cohorts/CAPSCAN/5_sage/razor_hash-acc/<s>`.
   - metaG = `canonical/publication_data/arms/metag_hashed_human/CAPSCAN/<s>`.
   - Production jobs were 5910681/5910682/5910683.
2. **Assignment export:** `reviews/supplementary_completion_2026-09-13/export_assignments.py` (fd88b9e4…, unchanged since the 282 run), tasks CAPSCAN_282_FL and CAPSCAN_282_metaG.
3. **Evidence ablation:** `Hartout_edits/H09_evidence_ablation/frozen/h09_evidence_ablation.py` (161aefe8…), job 5928243, giving `full_CAPSCAN/acquisition_summary.tsv`.
4. **SFig08 b–d table:** `MM_review_analysis/end_of_day_2026-09-15/code/plot_corrected_ablation.py` line 13 transform, giving `clinical_fixed_dictionary_ablation.tsv`. That table is a byte copy in the release repository.

The text prints no ablation or grouping values. The only printed quantities that depend on this chain are the roster sizes in paras 129, 235 and 331. The figure's dots and median bars change; the medians are reported below.

## What was done here

- **Code copies:**
  - `code/study.py`, `dispatch.py`, `h09_evidence_ablation.py` and `h01_eligibility.py` are byte-identical (sha256 checked).
  - `code/export_assignments.py` changes only its input and output locations; see `EXPORT_CHANGES.diff`.
  - `ablation_values.py` holds the plot_corrected_ablation transform for the CapScan rows. The PREDICT rows are unchanged.
- **Run folders:** `code/build_runs.py` builds `runs/ctl282` (the original manifest relocated, keeping the original search_inputs records) and `runs/r295`.
  - r295 is the same roster plus the 13 from `mg295/inputs/ACQ13.tsv`.
  - FL for the 13 comes from the same noclust razor_hash-acc root, which already holds them (2026-08-07).
  - metaG for the 13 comes from `mg295/mg_full/out/<s>`, the same-spectrum MG search.
- **Preflight:** the unchanged `dispatch.py preflight` passes for both manifests.
  - Negative preflight (`controls/NEGATIVE_PREFLIGHT.txt`): a changed input size and a changed roster hash both fail, as required.
- **282 control** (`code/compare_control.py` → `controls/CONTROL_282.json`): every grouping, export, H09 and SFig08-table output must equal the archived production files.
  - Comparison uses decompressed sha256.
  - Path- or timing-bearing JSON is compared without those fields.
- **Negative control** (`code/negative_control.py` → `controls/NEGATIVE_CONTROL.json`):
  - A planted one-digit change in the archived acquisition summary is detected.
  - The 295 output differs from 282.
- **Memory:**
  - Grouping: measured production MaxRSS × 1.05 (× 295/282 for r295). FL 3G; metaG 19G for the control and 20G for r295.
  - Export: the ctl282 export is the pilot at 48G. `jobs/tail295.sbatch` reads its MaxRSS and submits the r295 exports at MaxRSS × 1.05 × 295/282.
  - H09: 4G (production 2.9 GB).
- **Partition:** all jobs run on p.hpcl92 with `--nice=10000` and no `--begin`.

## Status

See JOBS.txt for job IDs. Collect with:

```
bash LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/phase3/item2_dictionary/collect.sh
```

This writes `item2_edits.jsonl`, but only if the control and the negative control PASS and r295 is complete.

## Results

Interim results, 7 Oct. The grouping stage is complete. The export, H09 and r295-tail jobs are queued or running (see JOBS.txt).

**Grouping control at 282:** the relocated run reproduces the production comparison exactly. `comparison/SUMMARY.json` and `matched_protein_group_counts.tsv` are byte-identical. The full per-file control will be in `controls/CONTROL_282.json` (job 6085767).

| claim | 282 value | 295 value | producer | control | job IDs |
|---|---|---|---|---|---|
| Matched-CapScan study dictionary: median quantified groups, FL / metaG (not printed in text) | 8,519.5 / 6,028.0 | 8,477.0 / 6,179.0 | dispatch.py collect | SUMMARY.json and counts byte-identical at 282 | 6085760-2 (282), 6085768-70 (295) |
| Ratio of medians FL/metaG; median paired ratio; FL higher (not printed) | 1.413; 1.245; 239/282 | 1.372; 1.233; 252/295 | same | same | same |
| Supp. Fig. 8b–d CapScan medians (not printed) | pending | pending | h09 + ablation_values.py | CONTROL_282.json | 6085766 (282); r295 IDs appended to JOBS.txt by tail 6085771 |
| Paras 129, 235, 331: roster "282" becomes 295 | 282 | 295 once H09 r295 completes | make_edits.py | requires CONTROL PASS and NEGATIVE PASS | - |
