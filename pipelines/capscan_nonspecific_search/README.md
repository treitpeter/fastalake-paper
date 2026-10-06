# Fully non-specific Sage searches of CapScan (A05)

Internal ID: A05 (previous folder name A05_unspecific_search)

Pipeline id `capscan_nonspecific_search`. Used by: Fig6 (non-specific route, panels a to f); text_statistics supplementary_methods_8_11.

## Question

Does the compact production FASTA permit fully non-specific searching of every CapScan acquisition, and which peptide termini and cleavage contexts do the accepted peptides show? Answer: one Sage 0.14.6 search per acquisition (331 original acquisitions, plus 16 acquisitions of participant 15 added later with the same code), each against exactly its production FASTA, with the digestion changed to fully non-specific and nothing else changed.

## Inputs

- Per acquisition: the production-selected FASTA, the mzML and the production tryptic Sage configuration (`source_config`), listed with size, mtime and sha256 in `inputs/manifest.json` (manifest sha256 `0a8ef92367beb27f50b94bafd4ae25c5a32dd9a255bc85787d405020c7ce3d6c`, 331 tasks, original run; the participant-15 run has its own manifest, see `RECEIPTS.md`).
- Baseline tryptic and semi-enzymatic search results of the same acquisitions (comparators, `baseline_search`, `semi_search` in the manifest).
- Sage 0.14.6 binary, sha256 `d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625`.

## Outputs

- `results/<task>/accepted_peptide_context.tsv`: accepted canonical peptides (target, rank 1, peptide q <= 0.01) with length, C-terminal residue, parent accessions, compatible enzyme contexts and an ambiguity flag.
- `results/<task>/{peptide_q01,both_q01,length9_50,target_all,decoy_q01}.txt` and `*_gained.txt`, `*_lost.txt` against the tryptic baseline.
- `results/<task>/COMPLETE.json`: validation receipt (counts, hashes of all outputs, Slurm job id) written only after every check passed.
- `results/<task>/search/results.sage.tsv`, `results.json`, `*.pin` and the `/usr/bin/time -v` record (large; not in the repository).

## Call chain and code

`lib.py`, `peptide.py` and `search.py` are shared by A13, A20 to A28 (byte-identical in every analysis directory of the chain; one copy is kept here). `prep_A05_p15.py` and `preflight_A05.py` built and checked the participant-15 run directory.

| file | origin (LAB_ROOT relative) | sha256 at origin | sha256 here |
|---|---|---|---|
| `lib.py` | `LAB_ROOT/MM_review_analysis/A05_4_7_unspecific/frozen/lib.py` | `a46c404fc225664624942ef1e380407054bab675ee63bad96fc944f3d74f5736` | `dd09acb87797d3e577f1dbc453dc851c55fb998b5bb4da6b6d9ceb801be56308` |
| `peptide.py` | `LAB_ROOT/MM_review_analysis/A05_4_7_unspecific/frozen/peptide.py` | `40e344f8bb3bc900b1528333496a35ec203bd85b79d81d15c339c61446f03eb1` | `40e344f8bb3bc900b1528333496a35ec203bd85b79d81d15c339c61446f03eb1` |
| `search.py` | `LAB_ROOT/MM_review_analysis/A05_4_7_unspecific/frozen/search.py` | `d3011bd2701bff862dbd609fae943d056c6bf599cdb6ed0aca025eca2cec4b52` | `d3011bd2701bff862dbd609fae943d056c6bf599cdb6ed0aca025eca2cec4b52` |
| `prep_A05_p15.py` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/scripts/prep_A05_p15.py` | `56d2357df3be3b4571b7aa39d38628606e24000649d7d03745f47017acdee01e` | `56d2357df3be3b4571b7aa39d38628606e24000649d7d03745f47017acdee01e` |
| `preflight_A05.py` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/scripts/preflight_A05.py` | `b27ca44a536db641f6845eb36ae6fe1fe53b474d9bf0967ef56f082a6e111183` | `b27ca44a536db641f6845eb36ae6fe1fe53b474d9bf0967ef56f082a6e111183` |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| Sage engine | Sage 0.14.6 binary, sha256 `d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625`, version verified from the program output before each run | one engine for all searches of the review analyses | other engines: campi_second_engine |
| digestion | `cleave_at` empty (fully non-specific), `restrict` null, `c_terminal` true, `semi_enzymatic` false; asserted in `search.py` before Sage starts | the experiment asks what is observed without a protease rule | tryptic and semi-enzymatic searches of the same acquisitions are the comparators (A15/A16 register, not rebuilt here) |
| peptide length | 5 to 50 residues (Sage), nested 9 to 50 summary (`length9_50`) | production setting; 9 to 50 matches the de novo comparison population | 9 to 50 is a sensitivity population and is not recalibrated |
| peptide mass | 500 to 5000 Da | production setting | none |
| fixed / variable modification | C +57.0215; M +15.9949, at most 2 variable modifications | production setting | none |
| tolerances | precursor and fragment +/- 20 ppm, isotope errors -1 to 3, charge 2 to 4, b and y ions | production setting copied from each acquisition tryptic `results.json` | none |
| decoys | `rev_` decoys generated by Sage | target-decoy confidence | none |
| accepted peptide | rank 1 target, `peptide_q` <= 0.01; `both_q01` additionally requires `spectrum_q` <= 0.01 | the identification gate used for every route in the chain | nominal Sage q-values do not establish whole-workflow FDR; stated in the study protocol |
| quantification | `lfq` false | identification and terminal-residue experiment only; the change is after peptide search and explicit | none |
| pilot gate | first capsule and first stool acquisition of the fixed roster must complete, recover every accepted peptide inside its reported parent and show at least one fully non-tryptic compatible peptide before the other tasks may start (`require_nontryptic`) | guards against a silent search-space error | none |
| enzyme context rule | N side compatible if peptide starts the protein, follows an initiator M, or follows K/R with next residue not P; C side compatible if the peptide ends the protein or ends in K/R with next residue not P | tryptic compatibility test used to label every occurrence fully, semi or non-tryptic | A20 re-annotates contexts with missed cleavages and no initiator-M clipping; the A05 clipping assumption is superseded there |
| threads / memory | 8 Sage threads (`RAYON_NUM_THREADS=8`), 640 GB per task for the full array | measured peak 471 GB (see Resources) | 160 GB failed out of memory on the stool pilot (job 5920921) |

## Resources and run time (measured)

Measured from Slurm accounting and the per-task `/usr/bin/time -v` records. Cohort run (jobs 5920952 pilot stool, 5920953 array of 329 tasks): 8 Sage threads (16 allocated hardware threads), 640 GB requested; peak resident memory 471.4 GB (stool pilot 412.4 GB), per-task wall time median 25.5 min, maximum 137.8 min, stool pilot 127.9 min. The first submission at 160 GB (job 5920921, 5920923) ran out of memory for the stool pilot (131.5 GB peak over the two pilot tasks) and was cancelled. Participant-15 run (job 6058023, 16 tasks): 405 GB requested, peak 269.7 GB, maximum 68.0 min.

Per-job accounting is in `RECEIPTS.md`.

## How to run

Every script takes a run directory as first argument and, for array stages, a task index as second argument. The run directory holds `inputs/manifest*.json` (task list with the path, size, modification time and sha256 of every input), `frozen/` (this pipeline's code, assembled by `capscan_nonspecific_search/slurm/assemble_run_dir.sh`) and `results/`.

    export FASTALAKE_PIPELINES_DIR=<repo>/pipelines
    export FASTALAKE_OUT_DIR=<run directory>      # contains inputs/manifest*.json
    export FASTALAKE_PYTHON=<python with the packages named below>
    sbatch --array=<first>-<last>%<concurrency> <flags from slurm/submit_example.sh> slurm/run_task.sbatch

The manifests list original instrument-derived inputs by absolute path and are not part of the repository (they name lab locations; their sha256 is in `RECEIPTS.md`). A re-run therefore needs a manifest rebuilt for the new location; the frozen code verifies every input against the manifest (path, size, modification time, sha256) and stops on any difference. Per the repository conventions this pipeline is not re-run for the paper; the code is the code that produced the outputs.

## Receipts and open items

`RECEIPTS.md` lists the original job ids and the sha256 of every output. `OPEN_ITEMS.md` lists what could not be traced.
