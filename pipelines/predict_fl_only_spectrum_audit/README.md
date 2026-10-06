# Spectrum audit of Fig. 4d on the public-only and specimen-metagenome arms (A32)

Internal ID: A32 (previous folder name A32_fig4d_spectrum_audit)

Pipeline id `predict_fl_only_spectrum_audit`. Used by: Fig4 (panel d; text line 131).

## Question

Why does FastaLake identify spectra that the specimen-metagenome comparator does not, on the same two databases as Fig. 4b and 4c? Answer: for 296 matched PREDICT acquisitions both databases are re-searched with rank-5 PSM reporting and every FL-only spectrum is assigned to one of five outcomes of the metagenome search; `scoring_counts_public_vs_specimen.tsv` and `spectrum_audit_per_acquisition.tsv` carry the partition.

## Inputs

- Per acquisition: the public-only FastaLake parsimony FASTA (`FL_fasta`), the specimen metagenome plus human Swiss-Prot FASTA (`MG_fasta`) and the mzML (`mzml`), listed in `inputs/manifest.json` (sha256 in `RECEIPTS.md`) and the 296-sample roster `inputs/samples_296.txt`.
- Sage 0.14.6 binary, sha256 `d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625` (checked in every task).
- `R3_digest.py` (archived digest module, copied here).

## Outputs

- `out/<acquisition>.json` (296 audit files), `scoring_counts_public_vs_specimen.tsv` and `SUMMARY.json` from `collect.py`; raw rank-5 Sage results per acquisition (`results/<acquisition>/{FL,MG}/results.sage.tsv.gz`, not in the repository).

## Call chain and code

`task.sh` runs the two Sage searches and `audit.py`; `collect.py` writes the summary. `R3_digest.py` is the archived digest module the audit imports (location of origin in the Fig4 R8A chain).

| file | origin (LAB_ROOT relative) | sha256 at origin | sha256 here |
|---|---|---|---|
| `audit.py` | `LAB_ROOT/MM_review_analysis/A32_4_5_fig4d_public_vs_specimen/code/audit.py` | `88809c03eca2ddd94fa05b244974d01ee96bb77d9726819cf5533deead4abd58` | `88809c03eca2ddd94fa05b244974d01ee96bb77d9726819cf5533deead4abd58` |
| `collect.py` | `LAB_ROOT/MM_review_analysis/A32_4_5_fig4d_public_vs_specimen/code/collect.py` | `2d87215ce4141056d71dfa355f0ae61de43a2d48cd3c24305d2b969cbe2faa96` | `2d87215ce4141056d71dfa355f0ae61de43a2d48cd3c24305d2b969cbe2faa96` |
| `task.sh` | `LAB_ROOT/MM_review_analysis/A32_4_5_fig4d_public_vs_specimen/code/task.sh` | `810bcc2ddcfbad2ce1cc223ecc8b807bf028d42c88c504aa7e440aef4d50ea94` | `810bcc2ddcfbad2ce1cc223ecc8b807bf028d42c88c504aa7e440aef4d50ea94` |
| `R3_digest.py` | `LAB_ROOT/canonical/final_revision/thread_selection_math/R3_digest.py` | `fc2a996fe0f7adddae1a19e087c603f1172ad2b97c2398ff6e85642c71d6b105` | `fc2a996fe0f7adddae1a19e087c603f1172ad2b97c2398ff6e85642c71d6b105` |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| Sage engine | Sage 0.14.6 binary, sha256 `d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625`, version verified from the program output before each run | one engine for all searches of the review analyses | other engines: campi_second_engine |
| search settings | identical to the Fig. 4 production arms except `report_psms` 5 and `lfq` false: KR/P, 2 missed cleavages, 5 to 50 residues, 500 to 5000 Da, C +57.0215, M +15.9949 (max 2), +/- 20 ppm precursor and fragment, isotope errors -1 to 3, b and y ions, 15 to 150 peaks, at least 4 matched peaks, `predict_rt` true, `chimera` false | the audit needs the best target row per spectrum among five candidates because the rank-1 row may be a decoy | copied from the `R8A_loss_psm5` configuration; production searches used `report_psms` 1 |
| acquisitions | 296 matched acquisitions (roster of A16 PREDICT tasks that have a metagenome configuration); 5 of the mzML files are the RAW-repaired files of predict_truncated_mzml_repair | arms of panels b and c use the same set | panel d as first published used 306 acquisitions and the full-reservoir arm |
| arm definitions | FL = public-only FastaLake (`abl_nometag_noclust`), MG = specimen metagenome plus human Swiss-Prot (`metag_hashed_human`) | the same two databases as panels b and c | none |
| categories | five MG-search outcomes for each FL-only spectrum (categories A to E; C, D and E use the comparator best target for the spectrum), computed in `audit.py` with the digest logic of `R3_digest.py` | explains the identification gap by spectrum | statistic is a count partition; no test |
| resources | 16 CPUs, 52,500 MB, 4 h per acquisition, array of 296 with at most 30 concurrent | measured peak of the specimen-MG production searches 51.1 GB times 1.05 | 2 tasks ran out of memory at 52.5 GB (peak 54.6 GB) and were rerun at 64 GB (job 5972160) |

## Resources and run time (measured)

Preflight job 5966850 (1 acquisition end to end with self and swap controls): 21.3 min, 20.5 GB. Array job 5966898: 296 tasks, 16 CPUs, 52,500 MB requested, median 8.8 min, maximum 48.1 min, peak 54.6 GB; 294 completed and 2 ran out of memory. Rerun of the two tasks, job 5972160 (64 GB): 14.4 and 17.4 min. First collect job 5966899 failed; collect job 5972161 (8 GB, 1 s) wrote the summary.

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
