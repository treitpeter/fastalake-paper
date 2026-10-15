# capscan_mg_benchmark_295

Internal ID: mg295 (working-directory name; no lab code assigned)

Question. Does the public-sequence FastaLake versus matched-metagenome comparison in CapScan hold when the 13 participant-15 acquisitions whose specimen has a completed public shotgun run are added, giving 295 matched acquisitions (263 specimens; 56 stool, 239 intestinal fluid) instead of 282 (250 specimens)? The 295 roster is the CapScan matched benchmark of the v33 manuscript (Fig. 4b,c, Supplementary Figs 4, 6, 8 and 9, Methods, Supplementary Methods 6-7, Discussion counts). See the DNA count reconciliation in `pipelines/capscan_participant_roster/README.md`.

Status: complete (6-7 October 2026). All Slurm jobs are in RECEIPTS.md.

## Roster

`inputs/build_inputs.py`: an acquisition is MG-matched when its specimen has one of the 301 active WGS runs, the same rule as for the 282. 13 acquisitions enter (participant 15, runs 1800-1810 and 1812-1813, all intestinal: duodenum 3, jejunum 3, ileum 4, colon 3; BioSamples SAMN27256486-98). Three participant-15 acquisitions stay out: 1799 (run excluded, mate 2 unavailable) and 1811, 1814 (amplicon-only BioSamples, no WGS run). `inputs/build_roster.py` appends the 13 rows (i = 578-590) to the production 578-row roster, whose lines are kept byte for byte.

## Code (as it ran)

| step | files | what |
|---|---|---|
| (a) MG FASTA | `db/build_mg_fasta.py`, `jobs/db_build13.sbatch` | production rehash plus human Swiss-Prot concatenation for the 13; control rebuilds an existing pair byte for byte, a one-residue change fails the sha256 gate |
| (b) same-spectrum MG search | `mg_full/search_mg.sh`, `mg_full/count.py` | copy of the production `search_metag_hashed.sbatch` (fasta, mzML, output changed only); frozen accepted-peptide rule |
| (c) 10-fold cross-fit | `crossfit/design/*`, `crossfit/retry_if_missing.sh`, `crossfit/CHANGES.diff`, `crossfit/COPY_SHA256.txt` | public347 cross-fit design; worker, sage_q and sage_lda byte-identical; roster and manifest paths changed only |
| (d) spectrum split MG_ODD | `disjoint/*` | iteration3 disjoint worker byte-identical, odd arm only for the 13 |
| (e) second engine | `second_engine/build_se_manifest.py`, `second_engine/task.sh` | MSFragger 4.1 and Percolator 3.06.5 MG tasks for the 13 through the unchanged `clinical_second_engine` task and frozen code |
| control | `controls/compare.py` | existing pair T760_1311 re-run through (a)-(e); see RECEIPTS.md |
| phase 2 | `phase2/*` | `per_acq13.py` (frozen confidence reader for the 13), `build_sd295.py` (13 rows appended to every CapScan-matched Source Data table; 282 rows byte for byte), `extract.py`/`harvest.sh`/`run_units.sh` (harvest of unit notebooks), `extra295.py` (mixed model, summaries), `secondengine295.py`, `regen_sfig09.py`; `PHASE2.md` |
| phase 3, item 1 | `phase3/item1_minlen7/` | min-length-7 MG searches for the 13 (byte copies of the minlen7_search_sensitivity task code) and collector |
| phase 3, item 2 | `phase3/item2_dictionary/` | matched-CapScan study dictionary, assignment export and evidence ablation at 295 (Supplementary Fig. 8b-d) |
| phase 3, item 3 | `phase3/item3_collision/` | manifest for the metagenome_header_audit worker at 295 (worker unchanged) |
| phase 3, item 4 | `phase3/item4_pep/` | PEP-gate summaries and cross-fitted PEP pooling at 295 |
| Discussion counts | `discussion_counts_295/` | sequence-accumulation and candidate-competition counts at 295 with the `discussion_counts` code (see `discussion_counts_295/RESULTS.md`) |

Each step folder keeps its own record of results and controls (`PLAN.md`, `phase2/PHASE2.md`, `phase3/PHASE3.md`, `phase3/item*/ITEM*.md`, `discussion_counts_295/RESULTS.md`). These records were written in the lab before the v33 text and refer to the v30 paragraph index; values are restated against v33 in the unit NUMBERS.tsv files.

Paths. The code is the code as it ran. Absolute lab paths are written `LAB_ROOT/<relative path>` (LAB_ROOT is the lab project directory, CONVENTIONS.md), the Python interpreter is `${FASTALAKE_PYTHON:-python3}` and the Sage executable `${FASTALAKE_SAGE:-sage}`. Nothing else was changed. `final/` below means `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final`. Manuscript-edit scripts of the working directory (`make_edits*.py`, `number_changes.py`) are not included.

## Outputs used by the repository

Source Data tables with 13 appended rows: `source_data/Fig4/{crossfit_per_acquisition,same_spectrum_benchmark,capscan_dna_yield_by_acquisition}.tsv`, `source_data/SFig06/same_spectrum_benchmark_1_capscan.tsv`, `source_data/SFig09/*` (CapScan tables), `source_data/supplementary_methods_6_7/{crossfit_per_acquisition,confidence_ordinary_all_pairs,minlen7_benchmark_per_acquisition,crossfit_pep_per_acquisition}.tsv`, `source_data/methods/crossfit_lda_fallback_347.tsv`, `source_data/SData1/v28d_capscan_mg_pairs.tsv`, `source_data/SFig04/{benchmark_matched_per_acquisition,per_acquisition_length_counts,search_inventory}.tsv`. Replaced summaries: `source_data/supplementary_methods_6_7/{crossfit_summary,crossfit_summary_by_region,crossfit_pep_summary,crossfit_pep_summary_by_region,confidence_ordinary_summary}.tsv`, `source_data/SFig08/clinical_fixed_dictionary_ablation.tsv` (CapScan rows), `source_data/supplementary_methods_1_5/header_audit.tsv` (CapScan rows). Output hashes are in RECEIPTS.md.

## Results (9-50 residues unless stated; 282 -> 295)

| claim | 282 | 295 |
|---|---|---|
| cross-fit intestinal gain | 28.3% (20.6-35.3) | 26.2% (17.0-34.3) |
| cross-fit regions | 21-36% | 20-35% |
| cross-fit stool | -5.5% | -5.5% (no new stool) |
| cross-fit recovery, intestinal | 75.4% | 75.5% |
| same-spectrum CapScan gain | 42.4% | 40.1% |
| acquisitions with higher FL counts | 281 of 282 | 294 of 295 |
| same-spectrum recovery | 86.4% | 86.5% |
| mixed model fluid vs stool | 1.40 (1.23-1.60), P=4 x 10^-7 | 1.37 (1.20-1.55), P=2 x 10^-6 |
| Supp. Fig. 9 tertiles | 70.4 / 58.3 / 42.0% | 69.6 / 52.9 / 42.4% |
| second engine CapScan (MSFragger) | 1.36 (281/282) | 1.35 (294/295) |
| min-length-7 CapScan gain | 45.0% | 43.8% |
| min-length-7 searches | 2,022 | 2,035 |
| name collisions (sequences) | 2,327,521 of 28,527,640 | 2,393,928 of 29,163,863 |
| PEP<=0.01 retention CapScan 5-50 | 84.6% | 84.7% |
| cross-fit PEP intestinal | 28.3% -> 36.3% | 26.2% -> 34.7% |
| Discussion CapScan MG union | 28.5 million | 29.2 million |
| Discussion candidates MG / FL | 3,323 / 278 | 3,318 / 276 |

PREDICT and CapScan stool values are unchanged.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| roster | 295 = 282 + 13 participant-15 acquisitions with an active WGS run | same matching rule as the 282 | 1799, 1811, 1814 excluded (`inputs/EXCLUDED3.tsv`) |
| 282 rows | kept byte for byte in every Source Data table | the 282 results are unchanged | none |
| MG database | production rehash + human Swiss-Prot, Sage binary reporting 0.14.6 | identical to the 282 MG arm | none |
| cross-fit fold caches | public347 caches (347-pooled) for the 13; the 282 rows keep the 331-pooled iteration4 caches | rerunning all 282 is unnecessary: pooled selections were identical for the control acquisition apart from ordinary Sage per-fold variation | item 6 below |
| memory | measured MaxRSS x 1.05; cross-fit retried at 36G for two tasks out of memory at 18G | lab Slurm rule | RECEIPTS.md |

## Limitations and notes

- Item 6, mixed fold-cache pools. The 282 cross-fit rows (and the 282 part of the cross-fitted PEP pooling) come from the 331-pooled iteration4 fold caches, the 13 new rows from the 347-pooled public347 caches. For T760_1311 run from both pools, baseline pooled sets and MG under PEP are identical; FL and FULL under PEP differ by ordinary per-fold Sage variation (7,395 vs 7,401; 8,306 vs 8,292 at 5-50 residues).
- Item 2, 282 control. `phase3/item2_dictionary/controls/CONTROL_282.json` reports status FAIL: the re-run H09 acquisition summary and the derived SFig08 table differ from the archived files only by floating-point noise (largest relative difference 1.8 x 10^-15, in the Bray-Curtis-derived columns); grouping outputs are byte-identical. The negative control passes. In the release table the PREDICT rows were kept as committed and the 295 CapScan rows of the re-run replace the 282 CapScan rows; CapScan values of the 282 acquisitions change because the study dictionary is built over the whole roster.
- The singular-fit fallback count ("three acquisitions") is unchanged: no fallback occurred for the 13 (`phase2` CROSSFIT13_NOTES.json).
- Participant 15 has no recorded DNA quantity for its capsule specimens, so the 99 acquisitions with recorded quantities are unchanged.

## How the outputs entered the release Source Data (v33 sync, 2026-10-07)

- `source_data/supplementary_methods_6_7/crossfit_summary.tsv`: CapScan rows from `phase2/crossfit_summary_591.tsv`; the 20 PREDICT rows are kept from the 282 release, because the producer run at 591 drew different PREDICT bootstrap resamples (point estimates identical, interval bounds differ by up to 0.15 points) while the repository code (`text_statistics/supplementary_methods_6_7/stats.py`) reproduces the release PREDICT intervals. `crossfit_summary_by_region.tsv` has no PREDICT rows and is the 591 file.
- `source_data/SFig08/clinical_fixed_dictionary_ablation.tsv`: the 296 PREDICT rows are the committed rows (the re-run differs only in the last floating-point digit); the 295 CapScan rows are the phase 3 item 2 re-run.
- `source_data/supplementary_methods_1_5/header_audit.tsv`: CapScan rows replaced by the phase 3 item 3 run295 output.
- `source_data/SFig04/search_inventory.tsv`: 13 rows (run `mg295`) from the phase 3 item 1 receipts 0001-0013.
- All other tables: the 282 committed rows unchanged and 13 rows appended (phase 2 `build_sd295.py`, phase 3 items 1 and 4).
