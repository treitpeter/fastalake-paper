# mg295 phase 2: recount at 295 CapScan acquisitions, v30.1 edit list, figure rebuild (6 Oct 2026, 19:40)

No manuscript text, DOCX or release repository was edited. The edits are in `edits_v301_mg295.jsonl`: 32 edits, each anchored once in `v30/visible_v30f.txt`.

## How the numbers were recomputed

1. **Source tables.**
   - `repo/` is a copy of `release_prep/fastalake-paper_v28`. Every CapScan-matched Source Data table gets 13 rows (`build_sd295.py`). The 282 existing rows are kept byte for byte.
   - Per-acquisition counts for the 13 come from the frozen Review_v16 Confidence reader (`per_acq13.py`).
   - **Control:** the reader reproduced all 36 confidence rows and all 4 SFig06 rows of 1311 exactly (`per_acq13/CONTROL.json`).
   - Cross-fit rows use the public347 `analyze_final.py` and `recovery.py` row definitions.
   - SFig09 derived tables are regenerated with the unit's own `derive_tables()`.
2. **Code.** Only the n-asserts were changed: Fig4 578/282 → 591/295; SFig09 282/250/226 → 295/263/239 and EXPECTED_BIN_N [80,79,80]; SFig06 BENCH_N. Diff: `diff -r repo_base repo`.
3. **Running the units.**
   - The unit notebooks (Fig4, SFig06, SFig09, Supp Methods 6–7, Methods, Abstract/Discussion) were executed by `extract.py`, which harvests the values in NUMBERS.tsv order.
   - **Control:** the harvest of the unmodified copy (`repo_base`) reproduces every committed recomputed value: Fig4 88/88, Methods 31/31, Supp Methods 6–7 228/228, Abstract 18/18.
   - The SFig09 baseline cannot be harvested because its pixel check needs the Söhne fonts. Its committed values serve as the 282 reference.
4. **Outputs.**
   - `NUMBER_CHANGES.tsv` (58 changed values, 46 printed).
   - `MIXED_MODEL_282_vs_295.json`, `SECOND_ENGINE_295.json`.
   - `crossfit_summary_591.tsv` and `crossfit_summary_by_region_591.tsv`, which are the Source Data summaries at 591.

Participant 15 has no DNA quantity recorded for its capsule specimens, so the "99 acquisitions with recorded quantities" figure is unchanged. Database size is the sanitizer n_kept; all 13 are present.

## What changes (9–50 residues unless stated)

| claim | 282 | 295 |
|---|---|---|
| Cross-fit intestinal gain | 28.3% (20.6–35.3) | 26.2% (17.0–34.3) |
| Cross-fit regions | 21–36% | 20–35% |
| Cross-fit stool (unchanged, no new stool) | −5.5% (−8.4 to −2.3) | −5.5% |
| Cross-fit recovery, intestinal | 75.4% | 75.5% |
| Same-spectrum CapScan gain | 42.4% | 40.1% |
| Same-spectrum CapScan, acquisitions with higher FL counts | 281 | 294 of 295 |
| Same-spectrum CapScan recovery | 86.4% | 86.5% |
| Same-spectrum regions | 47–64% | 46–62% |
| Same-spectrum intestinal, acquisitions with higher FL counts | 225/226 | 238/239 |
| Mixed model, fluid vs stool (cross-fit) | 1.40 (1.23–1.60), P=4×10⁻⁷ | 1.37 (1.20–1.55), P=2×10⁻⁶ |
| Mixed model, region P | 0.63 | 0.59 |
| Mixed model, same-spectrum | 1.39 | 1.35 (P=3×10⁻⁶; same pattern) |
| Supp. Fig. 9 tertiles | 70.4 / 58.3 / 42.0% (76 / 75 / 75) | 69.6 / 52.9 / 42.4% (80 / 79 / 80) |
| Specimens | 250 | 263 |
| Fluid database median | 107,000 | 105,000 |
| Participants with both compartments | 14/14 | 15/15 |
| Pooled yield, MG / FL | 34% / 53% | 36% / 54% |
| Second engine CapScan, MSFragger | 1.36 | 1.35 |
| Second engine CapScan, Sage on the same acquisitions | 1.43 | 1.41 |
| Second engine CapScan, n | 281/282 | 294/295 |

PREDICT values and the stool values are unchanged; this is checked in the harvest, where the PREDICT bootstrap stream is identical.

**For Peter:** the Abstract (para 17), "21–36% more in intestinal fluid", becomes "20–35%". The Abstract is frozen.

## Not recomputed

These still say 282. Each needs an extra analysis on the 295 roster, or the text must state that it stays on the 282 roster.

1. **Supp. Fig. 4a,b and Supp. Methods 7 (para 398), min-length-7 results.**
   - Affected values: "45.0% (CapScan)" and "84.0%".
   - Needs 13 MG Sage searches at min_len 7. FL_public min7 for the 16 already exists in public347/minlen7.
   - Cost: about 30 s each.
2. **Study-dictionary / grouping comparisons and Supp. Fig. 8b–d KO ablation (paras 129, 235, 331).** These need the matched-CapScan study dictionary re-inferred on 295.
3. **Assembly name-collision audit (para 288).**
   - Current values: 2,327,521 of 28,527,640, and 2,414,558 of 2,414,668.
   - Needs a recount over 263 specimens. The hashed FASTAs exist, so this is cheap.
4. **PEP-gate summaries (para 326).**
   - Affected values: "84.6% in CapScan"; cross-fitted PEP gains "−1.1%" and "36.3%".
   - They come from precomputed summary tables and per-fold PEP pooling. The ordinary rows are in `confidence_ordinary_all_pairs.tsv` at 295; the summaries and the cross-fitted PEP pooling were not rerun.
5. **Supp. Methods 7 singular-fit paragraph.** No new fallback occurred (`CROSSFIT13_NOTES.json`), so "three acquisitions" stays.
6. **The 13 cross-fits used the 347-pooled fold caches; the 282 Source Data rows come from the 331-pooled run (iteration4).** public347 showed the pool size does not change selection: pooled sets were identical apart from Sage per-fold nondeterminism in 2 rows. State this in Methods, or switch all 282 rows to public347. The latter changes a few union values only.

## Figure rebuild instructions

Copy these source tables from `mg295/phase2/repo/source_data/` into the release repository (each table carries 13 appended rows):
- `Fig4/{crossfit_per_acquisition,same_spectrum_benchmark,capscan_dna_yield_by_acquisition}.tsv`
- `SFig06/same_spectrum_benchmark_1_capscan.tsv`
- `SFig09/{capscan_assembly_dna_yield_gain_by_acquisition,capscan_acquisitions_9_50,a08_benchmark282_old,a08_benchmark282_new,capsule_assembly_gain_*}.tsv`
- `supplementary_methods_6_7/{crossfit_per_acquisition,confidence_ordinary_all_pairs}.tsv`
- `supplementary_methods_6_7/crossfit_summary*.tsv`, taken from `phase2/crossfit_summary*_591.tsv`
- `methods/crossfit_lda_fallback_347.tsv`
- `SData1/v28d_capscan_mg_pairs.tsv`

Then:
1. Apply the n-assert patch.
2. Update every README sha256/rows line and every NUMBERS.tsv v30 column.
3. Rebuild Fig4 (b,c point clouds and per-arm medians), SFig06 (d,e), SFig09 (a,b), Supp Methods 6–7, Methods and Abstract/Discussion with `paper.py build`, with `FASTALAKE_FONT_DIR` set. Fonts are absent here, so no render was made.
4. Run `paper.py check` and the figure-consistency check on the built DOCX.
