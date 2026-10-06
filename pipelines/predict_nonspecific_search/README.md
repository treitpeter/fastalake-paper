# Fully non-specific Sage searches of the 306 PREDICT stool acquisitions (A13)

Internal ID: A13 (previous folder name A13_predict_unspecific)

Pipeline id `predict_nonspecific_search`. Used by: SFig14 (cases); text_statistics supplementary_methods_8_11.

## Question

Does the non-specific route seen in CapScan (A05) also appear in an independent stool cohort? Answer: the same search as A05, run on all 306 primary PREDICT acquisitions with each acquisition own production chemistry. Code is the A05 code; nothing is duplicated here.

## Inputs

- Per acquisition: production-selected, unclustered FASTA, mzML and production tryptic configuration (`inputs/manifest.json`, manifest sha256 `f459bff715df1961b9a0eb5378c1a7fc3984502fe1b5e8280b96797f1be9f869`, 306 tasks).
- Sage 0.14.6 binary, sha256 `d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625`.

## Outputs

- Same files as A05 per acquisition (`accepted_peptide_context.tsv`, peptide sets, gained and lost against the tryptic search, `COMPLETE.json`, raw `results.sage.tsv`, `results.json`).
- SFig14 uses `results.sage.tsv` and `results.json` of the PREDICT acquisitions for the case selection.

## Call chain and code

No code of its own. The A13 `search.py`, `lib.py` and `peptide.py` are byte-identical to the A05 copies (md5 `b47487c5d19a830284e62680b3ac9cd7` for `search.py`), so the A05 copy is used. The three `release_capacity*.py` scripts in the original directory only released held scheduler capacity and produced no used output; they are not included.


## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| Sage engine | Sage 0.14.6 binary, sha256 `d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625`, version verified from the program output before each run | one engine for all searches of the review analyses | other engines: campi_second_engine |
| digestion | empty `cleave_at`, `restrict` null, `c_terminal` true, `semi_enzymatic` false; all other configuration fields equal to the source configuration | only the digestion changes | see A05 |
| chemistry, tolerances, length, mass | per acquisition, copied from its own production configuration (identical values for all 306: 5 to 50 residues, 500 to 5000 Da, +/- 20 ppm, charge 2 to 4) | keeps the comparison to each acquisition own tryptic search valid | see A05 |
| quantification | `lfq` false, decoys generated | identification experiment | see A05 |
| pilots | first acquisition of the frozen roster and the acquisition with the largest selected FASTA (bytes, sample id breaks ties), chosen before inspecting outcomes; both must pass confidence, containment and non-tryptic checks before the other 304 tasks start | second pilot is a resource stress test | none |
| resources | 8 threads, 640 GB, 8 h, at most five concurrent cohort tasks, two concurrent pilots | capacity shared with the CapScan array (at most fifteen search tasks together) | a 480 GB stress run (job 5932999) completed at 433.8 GB peak, so 480 GB was enough; 640 GB is what the cohort ran with |
| accepted peptide | rank 1 target, `peptide_q` <= 0.01, nested `spectrum_q` <= 0.01, 9 to 50 sensitivity | as A05 | nominal q-values only |

## Resources and run time (measured)

Pilot job 5926621 (2 tasks) and array job 5926622 (304 tasks), 8 threads (16 allocated hardware threads), 640 GB requested, peak resident memory 432.3 GB, per-task wall time median 32.3 min, maximum 133.1 min. Stress run 5932999 at 480 GB: 433.8 GB peak, 127.1 min.

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
