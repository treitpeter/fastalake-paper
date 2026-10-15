# metagenome_header_audit

Internal ID: sequence_accumulation header audit (review analysis of 29 September 2026)

Question. Across specimens, how often does an identical protein sequence carry different accessions, and how often does one accession name different sequences, in the public-sequence FastaLake databases and in the matched metagenome databases of PREDICT (296 specimens) and CapScan (263 specimens, 295 acquisitions; 250 specimens and 282 acquisitions before the v33 roster)? Metagenome accessions are audited three ways: as searched (hashed accessions), under the original sample-prefixed contig names, and with the sample prefix stripped (assembly-local names, from `k141_`). The audit also checks that rehashing preserved each specimen's sequence set after accounting for the appended human Swiss-Prot proteins (MG_SEQUENCE_PRESERVATION.json). Supplementary Methods 2 cites these tables in Source Data.

Code in this folder: worker.py (function `audit`, called for each namespace; the same worker.py run also writes the accumulation outputs, see pipelines/cohort_database_accumulation), plot.py (header_inconsistency_audit figure), peptide.py (imported by worker.py).
Launchers (slurm/): stage.sh, pilot.sbatch, full_PREDICT.sbatch, full_CAPSCAN.sbatch, figures.sbatch (identical to cohort_database_accumulation; one run produces both).
Environment variables: FASTALAKE_DATA_DIR (LAB_ROOT, the data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), FASTALAKE_PYTHON (interpreter). The run-time manifests list absolute input paths and are not redistributed; their sha256 is in RECEIPTS.md. worker.py finds LAB_ROOT as the fourth parent of its own directory, so the run directory must sit at FASTALAKE_DATA_DIR/review/v15_integration/analyses/<run name> (see slurm/stage.sh).

## Inputs

As pipelines/cohort_database_accumulation: LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/manifest.json (sha256 in RECEIPTS.md); the original sample-prefixed MG FASTAs named there (raw=True in `fasta`).

## Outputs

full/<cohort>/header_audit.tsv (one row per namespace), full/<cohort>/MG_SEQUENCE_PRESERVATION.json, figures/header_inconsistency_audit.pdf. Both header_audit.tsv files are in source_data/supplementary_methods_1_5/header_audit.tsv (concatenated with a cohort column, values unchanged). The CapScan rows are from the 295-acquisition re-run (below); the PREDICT rows are unchanged.

## CapScan re-run at 295 acquisitions (v33 roster)

The unchanged worker.py (sha256 331ad8e4...) was run on the production manifest plus the 13 participant-15 acquisitions of `capscan_mg_benchmark_295` (manifest built by `pipelines/capscan_mg_benchmark_295/phase3/item3_collision/build_manifest.py`; specimen id of the 13 = `15|<workbook_Main>`, 263 specimens). Control: the 282 re-run (job 6085730) reproduced the production header_audit.tsv, per_acquisition, accumulation, per_specimen and MG_SEQUENCE_PRESERVATION.json byte for byte. Negative control: a swapped raw assembly for one new row stopped at the worker gate (job 6085729). 295 run: job 6086073 (peak RSS 18.35 GiB, 16.9 min).

| CapScan, MG_assembly_local_headers unless stated | 282 | 295 |
|---|---|---|
| distinct sequences in more than one specimen under different names | 2,327,521 of 28,527,640 (8.2%) | 2,393,928 of 29,163,863 (8.2%) |
| shared assembly-local names denoting different sequences | 2,414,558 of 2,414,668 | 2,494,126 of 2,494,236 (>99.9%) |
| sequences with more than one accession, as searched (MG_as_searched) | 53 | 53 |
| accessions naming two sequences, as searched | 0 | 0 |

## Resources and run time (sacct)

Shared with cohort_database_accumulation: full PREDICT 6025426 (1 CPU, peak RSS 27.3 GB, 26.5 min), full CAPSCAN 6025427 (1 CPU, peak RSS 18.6 GB, 17.4 min).

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| header identity | first whitespace-delimited header token | the accession a search engine reports | none |
| sequence identity | SHA-256 of the exact amino-acid string | no I/L collapse, so only literal duplicates count | none |
| denominators | shared across specimens, and all distinct | rates among shared items show the inconsistency; all-distinct rates show its scale | both columns in header_audit.tsv |
| assembly-local name | header from the first `k141_` onward | strips the specimen prefix to test name reuse across assemblies | original prefixed names reported alongside |

## Known limitations

- 53 sequences per cohort carry more than one accession in the as-searched MG namespace, identically in both cohorts. They are in the as-searched files, which include the appended human Swiss-Prot; their cause was not traced (v30 Methods code check, open point 1).
- In CapScan (295 acquisitions), 2,393,930 MG sequences are shared across specimens. All 2,393,930 carry more than one name under the original prefixed headers, and 2,393,928 under assembly-local names (282 acquisitions: 2,327,523 and 2,327,521).

## Notes

The code is the original code as it ran. Only absolute paths and interpreter or log paths were replaced by environment variables, and module docstrings were added where missing (not to worker.py files, whose sha256 is checked by the downstream scripts and recorded in COMPLETE.json). Nothing was rerun for this folder apart from the CapScan 295 run described above.
