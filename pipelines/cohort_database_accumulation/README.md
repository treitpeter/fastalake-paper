# cohort_database_accumulation

Internal ID: sequence_accumulation (review analysis of 29 September 2026)

Question. How many distinct database sequences accumulate as specimens are added, for the public-sequence FastaLake databases and the matched metagenome-plus-human databases, in PREDICT (296 specimens) and CapScan (250 specimens, 282 matched acquisitions)? Sequence identity is the SHA-256 of the exact amino-acid string (no I/L collapse). The curve is given in the fixed specimen order and as the exact mean over random specimen orders (from incidence frequencies), with the per-specimen peptide identifications (canonical, peptide q <= 0.01, 5-50 residues) alongside.

Code in this folder: worker.py (per cohort: reads every searched FASTA and Sage result in the manifest, checks peptide sets against frozen sets, writes the accumulation curve, per-specimen and per-acquisition tables, the header audit and the sequence incidence matrices), plot.py (figures and figures/SUMMARY.json), polish_figures.py (display restyling), peptide.py (canonical peptide definition imported by worker.py), PROTOCOL.md (the protocol as written before the run).
Launchers (slurm/): stage.sh, pilot.sbatch, full_PREDICT.sbatch, full_CAPSCAN.sbatch, figures.sbatch.
Environment variables: FASTALAKE_DATA_DIR (LAB_ROOT, the data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), FASTALAKE_PYTHON (interpreter). The run-time manifests list absolute input paths and are not redistributed; their sha256 is in RECEIPTS.md. worker.py finds LAB_ROOT as the fourth parent of its own directory, so the run directory must sit at FASTALAKE_DATA_DIR/review/v15_integration/analyses/<run name> (see slurm/stage.sh).

The same worker.py run produced the header audit; that analysis is documented separately in pipelines/metagenome_header_audit. The Discussion readout of these curves (values after 50 and 100 specimens, per-specimen medians, the CapScan 347 FastaLake-only extension) is pipelines/discussion_counts, task A.

## Inputs

LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/manifest.json (578 paired searches: per acquisition the FL_public and MG Sage results, searched FASTAs and frozen expected peptide sets, with size and mtime; sha256 in RECEIPTS.md). Human Swiss-Prot used to validate the MG-plus-human sets: LAB_ROOT/MM_review_analysis/A01_4_1_database_sizes/inputs/human_swissprot_canonical_archived_2026_03.fasta. Per-file input receipts (sha256, bytes, mtime) are written to full/<cohort>/INPUT_RECEIPTS.json.

## Outputs

full/<cohort>/accumulation.tsv, per_specimen.tsv, per_acquisition.tsv, MG_SEQUENCE_PRESERVATION.json, FL_public_sequence_incidence.npz, MG_sequence_incidence.npz, COMPLETE.json; figures/SUMMARY.json and figures/sequence_accumulation_*.pdf. Hashes in RECEIPTS.md.

## Resources and run time (sacct)

pilot 6025421: 1 CPU, 2.4 min. full PREDICT 6025426: 1 CPU, 48 GB request, peak RSS 27.3 GB, 26.5 min. full CAPSCAN 6025427: 1 CPU, 48 GB request, peak RSS 18.6 GB, 17.4 min. figures 6025428: 0.5 min.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| unit of accumulation | specimen (union of its acquisitions' databases) | technical replicates and aliquots of one specimen are not new biological material | acquisition-level sizes in per_acquisition.tsv |
| sequence identity | SHA-256 of the exact amino-acid string, no I/L collapse | counts database entries as searched | none |
| order | fixed specimen order and exact mean over random orders | the fixed curve depends on order | both in accumulation.tsv |
| matched rosters | PREDICT 296, CapScan 282 acquisitions (250 specimens) | MG exists only for these | CapScan 347 FastaLake-only extension in discussion_counts task A |
| peptide rule | rank 1, target, peptide q <= 0.01, modifications stripped, I->L, 5-50 residues (9-50 reproduced) | canonical definition (peptide.py) | none |

## Notes

The code is the original code as it ran. Only absolute paths and interpreter or log paths were replaced by environment variables, and module docstrings were added where missing (not to worker.py files, whose sha256 is checked by the downstream scripts and recorded in COMPLETE.json). Nothing was rerun for this folder.
