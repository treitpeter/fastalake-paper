# cohort_candidate_competition

Internal ID: sequence_accumulation candidates (review analysis of 29 September 2026)

Question. On identical spectra, how many candidates does Sage score with the public-sequence FastaLake database and with the matched metagenome-plus-human database, in the 296 PREDICT and 282 CapScan matched acquisitions? The metric is Sage 0.14.6 `scored_candidates`: preliminary candidates with at least one fragment match, summed over tested charge and isotope hypotheses. The population is every rank-one spectrum ID reported by both arms, target or decoy, before confidence filtering. Each acquisition's canonical peptide set (5-50 residues, peptide q <= 0.01) is checked against the frozen set.

Code in this folder: candidates/worker.py (per acquisition: paired scored_candidates medians, ratio, peptide gain; cohort medians over acquisitions and over participants), candidates/full_gate.py (requires a passing pilot from the same worker.py, then runs the full set), candidates/plot.py (figure), peptide.py (canonical peptide definition imported by worker.py).
Launchers (slurm/): stage.sh, pilot.sbatch, full.sbatch, figures.sbatch.
Environment variables: FASTALAKE_DATA_DIR (LAB_ROOT, the data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), FASTALAKE_PYTHON (interpreter). The run-time manifests list absolute input paths and are not redistributed; their sha256 is in RECEIPTS.md. candidates/worker.py reads manifest.json from the parent of its own directory (the accumulation manifest, LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/manifest.json).

The Discussion recount of these values with a 9-50 stratum and the healthy-stool study-union comparison (tasks E1 and E2) is pipelines/discussion_counts; E1 reproduced all 578 acquisition medians of this pipeline.

## Inputs

manifest.json (578 paired searches; sha256 in RECEIPTS.md); per-file receipts of every results.sage.tsv read are in candidates/full/INPUT_RECEIPTS.json. Metric definition and the Sage source it was checked against: candidates/DEFINITION.json (source scoring.rs sha256 10793f6de36a597b8ae11456de5325ed2d4c2a7e2068dc5bc006385d4f5e4b07).

## Outputs

candidates/full/per_acquisition.tsv, cohort_median_acquisition.tsv, cohort_median_participant.tsv, COMPLETE.json; candidates/figures/candidate_competition_and_yield.pdf. Hashes in RECEIPTS.md.

## Resources and run time (sacct)

pilot 6025473: 1 CPU, 4 GB request, 8 s. full 6025474: 1 CPU, 6 GB request, peak RSS 0.38 GB, 11.8 min. figures 6025475: 8 s.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| metric | Sage 0.14.6 scored_candidates | the engine's own count of preliminary candidates with a fragment match | not all theoretical precursor-compatible candidates (DEFINITION.json) |
| spectra | rank-one spectrum IDs reported by both arms, target and decoy, before q filtering | same spectra in both arms | spectra reported by one arm only are kept as coverage counts in per_acquisition.tsv |
| summary | median per acquisition, then cohort median over acquisitions; also over participants | acquisitions repeat within participants | cohort_median_participant.tsv |
| peptide gain | 5-50 residues | as in the paired comparison | 9-50 in discussion_counts task E1 |

## Known limitations

- An observational association between search competition and yield, not an isolated causal effect and not a multiple-testing factor.

## Notes

The code is the original code as it ran. Only absolute paths and interpreter or log paths were replaced by environment variables, and module docstrings were added where missing (not to worker.py files, whose sha256 is checked by the downstream scripts and recorded in COMPLETE.json). Nothing was rerun for this folder.
