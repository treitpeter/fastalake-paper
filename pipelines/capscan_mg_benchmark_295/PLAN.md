# mg295: CapScan paired metagenome benchmark, 282 → 295 acquisitions (6 October 2026)

Everything is in `final/mg295/`. Production files were not changed. No manuscript text or DOCX was edited.
Job IDs: `JOBS_PILOT.txt` and `JOBS.txt`. Controls: `controls/CONTROL.json` and `db/control/CONTROL.json`.

## Roster

`inputs/build_inputs.py` writes `inputs/ACQ13.tsv`, `inputs/EXCLUDED3.tsv` and `inputs/CONTROL.tsv`.

An acquisition is MG-matched when its specimen (`workbook_Main`) has one of the 301 active WGS runs (`analyses/capscan_metag_analysis_2026-09-10/aggregate/active_runs.tsv`). This is the same rule as for the 282.

- **13 in:** participant 15 on runs 1800–1810 and 1812–1813 (all intestinal: duodenum 3, jejunum 3, ileum 4, colon 3). Each has 1 BioSample, SAMN27256486–98.
- **3 out:**
  - **1799:** run SRR18585253 is in `excluded_runs.tsv`: mate 2 unavailable, skip authorised 2026-09-11. An assembly does exist for SAMN27256485, but it is not used, by the same rule as the 282.
  - **1811 and 1814:** SAMN31647763 and SAMN31647764 are amplicon-only in ENA. There is no WGS run and no assembly.

The new total is 295 = 282 + 13. The 13 are all participant 15 and all intestinal, so stool stays at n = 56 and intestinal becomes 239.

## Inputs (sha256)

- **Human Swiss-Prot:** `databases/isoforms/human/uniprot/swissprot_human_reviewed_2026_03.fasta`, b205ed99…3e63.
- **Sage binary:** sha256 d3820543…5064625.
  - The binary reports `sage 0.14.6`, and all 282 production MG `results.json` say version 0.14.6.
  - This is the same binary as A08_347, split-half and cross-fit. No 0.14.7 binary is involved anywhere.
- **MSFragger 4.1 jar:** 9a16aa5a….
- **Percolator 3.06.5:** 1f1034b4….
- **mzML:** `capscan_unified/mzml_combined/<s>.mzML`, which is a symlink to `Dari_RAWs_Main_human/raw/mzML/<s>.mzML`.
  - This is the same file that A08_347, spectrum_split347 and the second-engine FL_public tasks use. The sha256 equals the second-engine record for all 13.

| acquisition | assembly (prefixed .faa) | faa sha256 | hashed sha256 (= production assemblies_hashed) | MG DB sha256 (`db/_db/<s>.metag_hashed_plus_human.fasta`) | MG entries | mzML sha256 |
|---|---|---|---|---|---|---|
| …1800 | SAMN27256486 | 81a6f4f6aed4c714 | 34325332667b84e4 | b4a1006d49dc94cc2f930f707b324971b9f222623e80c4e3ba866ce3304ddc90 | 56,179 | 9eca0ff07b51dae5 |
| …1810 | SAMN27256496 | 3d11419e79aa414a | 25bc115f6d0f3192 | 93b798ca9eeeb12d40c90245fa905e817d0f12eff8bf1b9c4d736bc75df10104 | 123,299 | 3a49c830410eaa05 |
| …1812 | SAMN27256497 | 9328e52bb3675bf8 | 843cf1faf0d9e909 | be785de4b70b9e260d069510fcdc8d7e80b4e2e26969a598ba9edbf2fbe6f427 | 14,431 | 6bf63da1faaae68a |
| …1801 | SAMN27256487 | 52d3f8c8033d0a49 | fe76c9399f399669 | 952ab56c4c3954c13b9f1f3a2cd3e97136b471a67bba598c943db7584af01813 | 104,126 | 134418496a93b229 |
| …1802 | SAMN27256488 | 2cf2cce4fa75f855 | 54e23e3bfcccb653 | 18299ddb2e3321802e06ec5296334f2e12f92594f25fc233540bc436d02ab22f | 130,266 | bb61c55b0c1ab965 |
| …1803 | SAMN27256489 | aea61ebe2dff035e | df5e65bf2f505cea | 855161087e3df4cac82dcdd6489346e7f3a5baf5453a8090fe3fa9077d8a3503 | 157,159 | 930c6e009673c97e |
| …1813 | SAMN27256498 | 59d37adcf7da1b98 | 62befc82a778cd4c | 9ef3216cf74854e32407c146b0cecce99ae614ddf95a8128fa6d780c1002aa4a | 93,922 | 5912434d9928659f |
| …1807 | SAMN27256493 | a6f8b04a1036489f | df9950a55205483a | 34d94dcb9fa2ada2bdc571db11504302295ebacf5faed17d455dfd8df8114932 | 100,715 | ac1459576d0f1775 |
| …1808 | SAMN27256494 | b3d34fc8bebcecc1 | 90c31aae7442c067 | 7cb8f9a2b7d3de13604275bedab9329425cf3e838e0c1bdb9d9aae2d9e280333 | 138,323 | f8d49457c0b680dd |
| …1809 | SAMN27256495 | 9eaad94356994983 | 4313f0d4291d6ae9 | 96924eacd4462b0e30cfd81ac2f90fed87b01cc0717125f16e9c7572cb32e5ac | 82,465 | 134a9f951dde4611 |
| …1804 | SAMN27256490 | 3295215006062e49 | 9527a7807f424dd9 | 02505cbd2c3a6b00ea4456d33ac4305cfb28aa1be6ef465d752d60d286cb6ebe | 84,934 | bde9cb73d172ea27 |
| …1805 | SAMN27256491 | 32fb9efe149dbae8 | ac221b3418d014ea | a7cd04bc9d3ea61da448b58c85fd650aa4e0ed3d4efae2fa93db7bab366785f8 | 83,131 | 1f5c12698069ddd3 |
| …1806 | SAMN27256492 | 899b27c56110f4fb | 6d005aae6fc96215 | d98d18abb9fa597f919bc8696054e3461ef8c31f11fb1d0d77bdb2a6853a0da6 | 85,194 | 14c3f7e80c05fa7a |

Full sha256 values and byte counts are in `db/BUILD13.json`.

## What was built, and how it relates to production

**(a) MG FASTA: `db/build_mg_fasta.py`**
- The script follows the production chain:
  - `rehash_metag_headers.py` pass 2 for one sample, using production `read_fasta`, imported.
  - Then `build_metag_human_fasta.py`, which concatenates the human Swiss-Prot file and checks the header count with the imported `count_headers`.
- **Control:** I rebuilt the production specimen SAMN27256279 (control acquisition …1311).
  - The rehashed FASTA sha a7ab1199… equals the production file.
  - The final DB sha 41793854… equals the production `metag_hashed_human/CAPSCAN/_db` file.
  - Status: PASS.
- **Negative control:** a rebuild with one residue changed fails the sha gate, as it should.
- **The 13:** for each, the rehashed FASTA equals the production `assemblies_hashed` file, which already existed for all 335 specimens.

**(b) Same-spectrum MG Sage search: `mg_full/search_mg.sh`**
- This is a copy of `canonical/publication/arms/search_metag_hashed.sbatch`.
- The config-copy block is byte-identical to production. It changes only fasta, mzML and output.
- The template config is the control's `metag_solo` config. All 283 `metag_solo` CAPSCAN configs are identical apart from those three fields, which I checked.
- `count.py` applies the frozen accepted rule.
- Outputs: `mg_full/out/<s>/{results.sage.tsv, results.json, ACCEPTED_{5,9}_50.txt, COUNTS.json}`.

**(c) 10-fold cross-fit: `crossfit/design/`**
- These are copies of `public347/crossfit_capscan347/design`. `worker.py`, `sage_q.py` and `sage_lda.py` are byte-identical.
- `common.py` changes only two things (see `crossfit/CHANGES.diff`):
  - ROSTER points to `disjoint/roster_mg295.tsv`.
  - MAN is read in place from the public347 manifest.
- The fold rule is unchanged: SHA-256 of `FastaLake-D11-crossfit-20261004:<scan>`, mod 10.
- The caches are the public347 fold caches, pooled over 347 and already containing the participant-15 evidence, used through a symlink.
- Arms are FL, FULL and MG. Pooled rescoring (LDA plus picked-peptide q) runs per acquisition inside the worker.
- Outputs: `crossfit/outputs/CAPSCAN/<s>/`.

**Roster: `disjoint/roster_mg295.tsv`, sha 58079533…**
- The production 578-row roster, byte-identical lines, plus 13 rows with i = 578–590.
- Built by `inputs/build_roster.py` with the same derivation as `iteration3/disjoint/build_roster.py`.

**(d) Spectrum-split MG_ODD (Supp. Fig. 6d,e)**
- This arm is per-acquisition MG, so it is run here.
- It uses `iteration3/disjoint/worker.py` copied byte-identically, through `run_disjoint.py` with `FL_ROSTER`.
- For the 13 it runs the `odd` arm only. Their MG_FULL search is (b), because no archived frozen MG set exists for them.
- The odd-scan mzML must be byte-identical to `spectrum_split347` PARITY: 13/13 pass.
- Outputs: `disjoint/outputs/CAPSCAN/<s>/MG_ODD/`.

**(e) Second engine (MSFragger plus Percolator)**
- FL_public for these 13 already exists: `second_engine_all` tasks 1610–1625.
- MG tasks were built here by `second_engine/build_se_manifest.py`, using the production rule. FASTA and source_config come from the (b) `results.json`, and raw is the FL_public task's mzML.
- The tasks run through the unchanged `second_engine_all/code/task.sh` and `frozen/` code.
- Outputs: `second_engine/results/CAPSCAN__<s>__MG/`.

## Control: existing pair …1311 (roster i = 335) re-run with the mg295 scripts (`controls/CONTROL.json`: PASS)

| step | archived | rerun | identical set |
|---|---|---|---|
| (a) MG DB sha | 41793854… | 41793854… | yes |
| (b) MG same-spectrum 5–50 / 9–50 | 8,433 / 6,595 | 8,433 / 6,595 | yes / yes (psms q0.01 15,650 = 15,650) |
| (d) MG_ODD 5–50 (iteration3) | 4,947 | 4,947 | yes; odd mzML identical; MG_FULL equals frozen 5–50 and 9–50 |
| (c) cross-fit pooled FL / FULL / MG (public347, same caches) | 9,461 / 10,620 / 8,369 | 9,461 / 10,620 / 8,369 | yes |
| (c) union, diagnostic only | FL 9,486 / FULL 10,605 | 9,483 / 10,610 | no; known Sage per-fold variation, as in the public347 control |
| (e) second-engine MG 5–50 / 9–50 | 8,237 / 6,722 | 8,237 / 6,722 | yes |

## Pilot: …1800 (duodenum), all steps run end to end, outputs complete

| step | result | MaxRSS | elapsed |
|---|---|---|---|
| (b) | MG 23,097 / 18,911 (FL_public A08 29,208 / 23,199) | 7.0 GB | 1:49 |
| (d) | MG_ODD 13,229 | 9.7 GB | 1:06 |
| (c) | pooled FL 25,514, FULL 29,233, MG 23,112 (5–50) | 15.5 GB (control 14.5 GB) | 26:53 |
| (e) | MG 22,229 / 19,056 | 13.2 GB (array max 21.2 GB, below the production 26G request) | 1:53 |

Memory requests:
- **(b) and (d):** 24G, from the iteration3 MEM_LAW at the largest new DB (≈22.4 GiB).
- **(c):** 18G = ceil2(15.5 × 1.05 + 0.62). There is an automatic retry at 36G that runs only for a task missing COMPLETE.json.
- **(e):** 26G, the production request.

## Jobs

Pilot jobs: 6074569 (b), 6074570 (d), 6074571 (c), 6074572/6074573 (e), 6074574 (compare). DB build: 6074553.

Remaining 12: b 6075257, d 6075258, e manifest 6075259, e 6075260. All of these are COMPLETED, 13/13 outputs present.

Cross-fit:
- 6075396, array 0–11 %8, 18G.
- 6075405, the retry-if-missing at 36G, afterany.
- Expected finish about 17:40–18:00 CEST, given about 27 min per task with 8 concurrent.

Partitions are p.hpcl92, p.hpcl8 and p.hpcl93, all with `--nice=10000` and no `--begin`, as Peter instructed via the coordinator on 6 Oct.

## Phase 2: what must be recomputed at n = 295 (282 + 13)

1. **Fig. 4b,c, cross-fitted gains and recovery.**
   - Producer: `analyze_final.py` and `recovery.py` (public347 copies). Merge the 282 public347 cross-fit outputs with `mg295/crossfit/outputs` for the 13.
   - Figures: SUMMARY CapScan all / intestinal / stool, by region, and recovery 72.0/75.4%.
   - Also the abstract cross-fit gains, Methods line 267, and Suppl. Methods 6–7.
2. **CapScan same-spectrum gains.**
   - Figures: 51.7% / 42.4%, intestinal 64.4 / 54.9, regions 47–64%, "282/282", "281 of 282", "225 of 226", identity recovery 84.3%.
   - Inputs: FL_public = A08_347 for the 13; MG = `mg_full/out`. Use the `public347/collect` benchmark code path.
3. **Region mixed model.** `capscan_region_mixed_model`, ln(FL/MG) by site with a participant random intercept, n = 295.
4. **Supp. Fig. 6d,e.** MG_ODD-based gains (FULL / JOINT / DISJOINT vs MG_ODD), using spectrum_split347 FL halves for the 13 and `iteration3/disjoint/analyze.py` at n = 295.
5. **Supp. Fig. 9.**
   - DNA tertiles and yield ratios: `capscan_metagenome_yield_join` plus SFig09, adding 13 rows.
   - Join inputs: DNA amount per specimen, and MG DB size from `assemblies_prefixed/CAPSCAN_sanitize_report.tsv` n_kept. All 13 SAMN are present there.
6. **Fig. 4d audit.** This concerns PREDICT only (curated MG), so it does not change.
7. **Second-engine CapScan ratios.** `second_engine_all/code/collect_all.py`, paired FL_public vs MG, 282 → 295, using the `mg295/second_engine` results for the MG of the 13.
8. **Spectra / database / ID table.** The CapScan MG rows: database sizes, IDs and n.
9. **Every "282" in the text, legends and Source Data.**
   - Also: 250 BioSamples → 263, intestinal 226 → 239, "fourteen participants (intestinal)" → check.
   - Also the DNA_COUNTS / V30N-2 sentence stating that the 16 participant-15 acquisitions are not in the 282: it becomes 13 in and 3 out.

## Status 19:45: phase 1 COMPLETE, phase 2 done (see phase2/PHASE2.md)
- All 13 outputs exist for (a)–(e).
- Cross-fit: two tasks ran out of memory at 18G (…1803, …1807). The 36G retry job 6075405 completed them at MaxRSS 17.5 GB. For future runs, request 20G.
- The Sage binary reports 0.14.6. v30 Methods already describes it as the "v0.14.7 release binary, reporting version 0.14.6".
- **Phase 2 outputs:**
  - `phase2/edits_v301_mg295.jsonl`: 32 anchored v30.1 edits.
  - `phase2/NUMBER_CHANGES.tsv`
  - `phase2/PHASE2.md`: what changed, what was not recomputed, and the figure rebuild steps.
