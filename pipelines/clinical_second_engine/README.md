# clinical_second_engine

Internal ID: SE-ALL (lab folder final/second_engine_all; extends A07, campi_second_engine)

Question. Does the FastaLake peptide gain over the matched metagenome databases hold with a second search engine? MSFragger 4.1 + Percolator 3.06.5 on the main-benchmark FASTAs and mzML of PREDICT (FL_public and MG, 296 acquisitions each), CapScan (FL_public 347, FL_cohortMG 347, MG 282) and healthy stool (FL, 58). 1,626 tasks (one acquisition x arm each). Supports Results paragraph 43 and Methods paragraph 133 / Supplementary Methods paragraph 359 of the v30 text (paragraph index of visible_v30f.txt).

Pointer note (mismatch). The v30 text names `pipelines/campi_second_engine` in its analysis pointers (Methods paragraph 137 list; paragraphs 133 and 359 describe the clinical and healthy-stool searches without a separate pointer). The clinical and healthy-stool code is in this folder, `clinical_second_engine`. Either the v30 pointer adds `clinical_second_engine` or this folder is cited from `campi_second_engine`; author decision.

Code in this folder: build_manifest.py (stage inputs, hash every FASTA and mzML), bench_pick.py (benchmark task choice), task.sh (one task: frozen search, then receipt), receipt.py, append_p15.py (16 CapScan participant-15 FL_public tasks), collect_all.py (paired FL/MG and healthy-stool results), collect_fl347.py (CapScan FL_public at 347), preflight/make_preflight.py.
The search itself is the frozen CAMPI code, unchanged: search.py and closed_fragger.params from `campi_second_engine`, lib.py and peptide.py from `campi_near_miss_audit`; `slurm/stage_frozen.sh` copies them to `${FASTALAKE_OUT_DIR}/frozen/` as the run had them.
Launchers (slurm/): build_manifest.sbatch, task_array.sbatch (index = SLURM_ARRAY_TASK_ID + SE_OFFSET), collect.sbatch, append_p15.sbatch, preflight_p15.sbatch, collect_fl347.sbatch, stage_frozen.sh, truncation_check.sh; preflight/preflight.sbatch and make.sbatch.
Environment variables: FASTALAKE_DATA_DIR (data root laid out as LAB_ROOT), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON; truncation_check.sh also FASTALAKE_PREDICT_MZML_DIR and FASTALAKE_RAW_SEARCH_ROOTS. java and Percolator are taken from PATH and the run manifest; the MSFragger jar is not redistributed.

Run order: build_manifest -> preflight (CAMPI S02 FL and S03 MG reproduce the archived A07 accepted-peptide lists byte for byte; missing mzML and wrong FASTA sha256 fail as required) -> task_array per cohort x arm -> collect -> append_p15 (after its own preflight) -> task_array 1610-1625 -> collect_fl347.

## Inputs

inputs/manifest.json: per task the FASTA and mzML recorded in the Sage results.json of the main-benchmark search for that acquisition and arm, each with bytes and sha256 (manifest.tsv lists all 1,626). PREDICT MG FASTA sha256 equal the Fig. 4 benchmark record for all 296. Sage comparison counts come from the Fig. 4 benchmark Source Data tables and, for healthy stool and CapScan FL347, from each task's source results.sage.tsv (rank 1, target, peptide_q <= 0.01, canonical). Manifest hashes are in RECEIPTS.md.

## Outputs

RESULTS.md, RESULTS_CAPSCAN_FL347.md, results_per_task.tsv, results_paired.tsv, results_hstool.tsv, results_capscan_fl347.tsv and results/<task>/RECEIPT.json (hashes in RECEIPTS.md).

Headline (RESULTS.md): median FL/MG peptide-count ratio at 9-50 residues 1.051 in PREDICT (291 of 296 paired; Sage on the same acquisitions 1.095) and 1.364 in CapScan FL_public (281 of 282; Sage 1.428); healthy stool MSFragger/Sage 0.998 (58 of 58).

Failed tasks (13 of 1,626, listed in results_per_task.tsv): 5 PREDICT acquisitions in both arms (MicrobPred_13, 44, 90, 115, 181; MSFragger exits on these inputs, which are the five reconverted mzML files; see slurm/truncation_check.sh, job 6068722), and CapScan T760_1341 (FL_public and FL_cohortMG; both tasks failed the frozen Percolator minimum of 100 decoy PSMs in search.py, not memory: logs SE_FL_CAPSCAN_6064501_71.err "AssertionError: Counter({'1': 797, '-1': 77})" and _402.err "AssertionError: Counter({'1': 785, '-1': 69})"; a low-yield acquisition, not rerun because the frozen minimum would have to change) and T760_2200 (FL_public, task 906; out of memory at 14 GB).

T760_2200 FL_public was rerun with 32 GB (Slurm job 6098609, COMPLETED) and completes the CapScan matched set at 295 acquisitions: median MSFragger FL_public/MG ratio at 9-50 residues 1.3466 (Sage on the same acquisitions 1.4011); FL_cohortMG/MG 1.3475 (Sage 1.4035); FastaLake more peptides in 294 of 295 (lab file final/mg295/phase2_v2/SECOND_ENGINE_295_v2.json).

## Resources and run time (from sacct)

| stage | tasks | cores | memory request | peak RSS | longest task |
|---|---|---|---|---|---|
| FL PREDICT | 296 | 8 | 14 GB | 12.75 GB | 5.0 min |
| FL CapScan (public + cohortMG) | 678 | 8 | 14 GB | 11.32 GB | 4.0 min |
| FL CapScan participant 15 | 16 | 8 (16 threads) | 14 GB | 8.13 GB | 3.3 min |
| FL healthy stool | 58 | 8 | 14 GB | 11.69 GB | 4.1 min |
| MG PREDICT | 296 | 8 | 26 GB | 23.44 GB | 9.5 min |
| MG CapScan | 282 | 8 | 26 GB | 22.64 GB | 6.5 min |

Search arrays used 82.1 task-hours of wall time in total. Manifest build 4.2 min (32 cores, 16 GB, hashes about 1 TB). Collectors under 4 min (2 cores, 4 GB).

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| engine | MSFragger 4.1 closed search, frozen CAMPI code and params | same workflow as the CAMPI second-engine analysis | campi_second_engine |
| validation | Percolator 3.06.5 target-decoy, --seed 913, peptide q <= 0.01 | per acquisition and engine, as Sage | none |
| decoys | one fully reversed sequence per target, rev_ prefix | as campi_second_engine | Sage reversal (terminal-preserving) differs; Supplementary Methods 6 |
| digestion and ranges | 5-50 residues, 500-5,000 Da, fixed C, variable M, from the Sage search | matched chemistry | closed_fragger.params |
| inputs | exactly the FASTA and mzML of the main-benchmark Sage search per acquisition and arm | same databases and spectra as Fig. 4 | none |
| headline length | canonical 9-50-residue peptides; 5-50 alongside | Main-text headline rule | RESULTS.md |
| summary | per-acquisition FL/MG ratio, then cohort median | same estimator as the Sage benchmark | none |
| CapScan arms | FL_public/MG paired at 282 (headline); FL_cohortMG reported separately; FL_public at 347 vs Sage separately | Fig. 4b uses FL_public | RESULTS.md, RESULTS_CAPSCAN_FL347.md |

## Notes

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables, `import os` was added where a path became an environment lookup, and missing docstrings were added. Nothing was rerun for this folder.
