# discussion_counts

Internal ID: discussion_counts (final review, 6 October 2026)

Question. The counts quoted in the Discussion, all from existing search outputs (no new search): (A) sequence accumulation readouts and per-specimen peptide medians; (B) the share of GMSC small proteins among selected proteins and accepted peptides, and GMSC-only peptides absent from the matched metagenome database; (C) for one worked healthy-stool acquisition, how many searched, candidate and reservoir proteins contain each accepted peptide; (D) the CAMPI clustering cost at 9-50 residues; (E) Sage scored_candidates per rank-one spectrum (E1: 578 paired clinical acquisitions, E2: healthy-stool own versus study-union FASTA); and (spectra_db_ids/) per acquisition and arm, MS2 spectra, database proteins and residues, candidates, accepted peptides, peptides per protein and sequence coverage for 701 acquisitions. COUNTS.md (the results write-up) and spectra_db_ids/SUMMARY.md give every value and definition.

Peptide rule for all tasks: Sage rank 1, target, peptide q <= 0.01, modifications stripped, I->L, distinct; 9-50 residues primary, 5-50 alongside (common.py `accepted`, identical to peptide.py `canonical`).

Code in this folder: common.py (shared readers and paths), a_accumulation.py (A), b_gmsc.py and b_extra.py (B), c_homologs.py (C), d_clustering.py (D), e_candidates.py (E1, E2), peptide.py (canonical peptide definition); spectra_db_ids/tasks.py (task list from existing results.json files), spectra_db_ids/worker.py (per acquisition and arm), spectra_db_ids/collect.py (tables and checks).
Launchers (slurm/): stage.sh, run_tasks.sh (A-E), spectra_db_ids.sh and spectra_db_ids_worker.sbatch.
Environment variables: FASTALAKE_DATA_DIR (LAB_ROOT, the data root laid out as the lab tree; read by common.py and spectra_db_ids/tasks.py), FASTALAKE_OUT_DIR (run directory; must be `<final>/discussion_counts` beside the review `final` directory, because common.py and collect.py read sibling folders such as public347/A08_347, clustering_check and len950/SFig07/work), PIPELINE_DIR (this folder), FASTALAKE_PYTHON (interpreter; the original runs used the environment in environment/).

Related pipelines: cohort_database_accumulation (source of A), cohort_candidate_competition (E1 reproduces its 578 medians), campi_candidate_clustering and campi_clustering_recount (source of D).

## Inputs

All read-only, with sha256 or size and mtime checked by the code where stated in COUNTS.md:
- LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/manifest.json and full/<cohort>/ outputs (A, B, E1, spectra_db_ids).
- Public-sequence FastaLake and matched MG searches: PREDICT LAB_ROOT/canonical/publication_data/arms/abl_nometag_noclust; CapScan final/public347/A08_347 and LAB_ROOT/MM_review_analysis/A08_4_5_capscan_public_only; healthy stool LAB_ROOT/canonical/publication_data/arms/noclust_2026-08-07 (roster LAB_ROOT/simplified_grouping/publication/analyses/simplified_grouping_2026-09-10/hstool_samples.txt). Searched FASTAs are read from each search's results.json.
- Reservoir (C): the 233,497,702-sequence public reservoir and the worked acquisition's candidate and evidence FASTAs.
- D: LAB_ROOT/review/v15_integration/analyses/clustering_90_95_2026-09-29/manifest.json (105 Sage outputs) and final/clustering_check/per_acquisition_recount.tsv.
- E2: LAB_ROOT/MM_review_analysis/A06 (healthy-stool study-union searches) and the 58-file joint search.
- spectra_db_ids: common-FASTA manifest LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/common_cache/comparison/inputs/manifest.json, A21 manifest_revision2.json, the 347 roster, and the mzML files named in each results.json.

## Outputs

A_accumulation_readout.tsv, A_per_specimen_ids.tsv, A_per_specimen_medians.tsv, A_RESULT.json; B_per_acquisition.tsv, B_summary.tsv, B_extra.txt, B_gmsc_only/<cohort>/; C_per_peptide.tsv, C_RESULT.json; D_per_acquisition.tsv, D_summary.tsv, D_RESULT.json; E1_per_acquisition.tsv, E1_summary.tsv, E1_RESULT.json, E2_hstool_per_acquisition.tsv, E2_RESULT.json; spectra_db_ids/per_acquisition.tsv, summary_medians.tsv, summary_paired_ratios.tsv, CHECKS.json. Hashes in RECEIPTS.md. Pilot outputs (*_PILOT.*) are not used.

## Resources and run time (sacct)

| task | job | cores | peak RSS | wall time |
|---|---|---|---|---|
| A | 6073021 | 2 | 5.6 GB | 2.0 min |
| B | 6073023 | 32 | 0.67 GB | 1.5 min |
| B extra | 6073075 | 2 | < 0.01 GB | 1 s |
| C | 6073031 | 64 | 0.40 GB | 2.2 min |
| D | 6073040 | 2 | < 0.01 GB | 25 s |
| E1 | 6073052 | 32 | 1.7 GB | 45 s |
| E2 | 6073053 | 16 | 2.1 GB | 2.1 min |
| spectra_db_ids worker | 6073330 (100-task array) | 1 per task | 1.01 GB (largest task) | 1.1-1.6 min per task |

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| peptide length | 9-50 residues primary, 5-50 alongside | main-text comparisons use 9-50 | 5-50 columns in every output |
| GMSC source rule | first header token; kept header is the highest-priority source (UHGG > GMGC > GMSC > human Swiss-Prot) | the reservoir keeps one header per distinct sequence | none |
| GMSC-only peptide | every searched protein containing it (I/L-equivalent substring) is GMSC-tagged | attributable to GMSC alone | B_per_acquisition.tsv |
| absent from MG | not an I/L-equivalent exact substring of any sequence in the matched MG-plus-human FASTA searched | literal sequence absence, not spectral incompatibility | none |
| homolog unit (C) | distinct sequences containing the peptide (I/L-equivalent) | counts the alternatives a peptide could explain | searched, candidate and reservoir levels |
| CapScan roster | 347 acquisitions for FastaLake-only values; 282 for matched comparisons | MG exists only for 282 | A_RESULT.json |
| candidate population | rank-one spectrum IDs reported by both arms, before q filtering | same spectra in both arms | E1_per_acquisition.tsv |

## Known limitations

- COUNTS.md names job 6073023 for task A and 6073031 for task C. sacct shows A ran as 6073021, B as 6073023 and C as 6073031; RECEIPTS.md uses the sacct values.
- PREDICT (306) and CapScan (331) study-union joint searches have no Sage output (out of memory), so E2 covers healthy stool only (COUNTS.md, section E2).
- For 4 of 282 CapScan acquisitions the A08_347 peptide set differs from the production A08 set; A08_347 is used throughout (COUNTS.md, section B).

## Notes

The code is the original code as it ran. Only the absolute LAB_ROOT path (common.py, spectra_db_ids/tasks.py, d_clustering.py) and the interpreter and log paths in the launchers were replaced by environment variables. Tasks A-E were run as interactive srun steps; slurm/run_tasks.sh reproduces them with the recorded job names, cores and memory. Nothing was rerun for this folder.
