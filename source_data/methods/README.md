# Source data: Methods statistics

All tables were written by `text_statistics/methods/prepare_source_data.py` (Slurm job
`methods_prep`, 2026-10-06). Origins and their sha256 are in
`text_statistics/methods/PROVENANCE.md`. The unit also reads, without copying,
`source_data/Fig1/acquisition_sizes_726.tsv`, `source_data/Fig4/crossfit_per_acquisition.tsv`,
`source_data/Fig4/clinical_paired_statistics.tsv` and `source_data/SFig15/*.tsv`.

## reservoir_sources.tsv
sha256: 252b5daf420f520f1f4d60f8f301e3e66f78c54d2ed4b3f5343898f3397f18ed
rows: 19 (one row per reservoir and source)
produced by: cohort_database_sizes

| column | unit | meaning |
|---|---|---|
| reservoir | id | public (healthy stool and public-sequence benchmarks), PREDICT_cohort, CapScan_cohort, CAMPI |
| source | id | sequence source as named by the lake builder (UHGG, GMGC, GMSC, HUMAN_SP, PREDICT_MG, CAPSCAN_MG, UNIPROT_SP, UNIPROT_TR) |
| priority | rank | order in which sources were added; a sequence counts for the first source that contains it |
| files | count of FASTA files | input FASTA files of the source (assembly FASTAs for *_MG) |
| input_records | count of sequences | records read from the source |
| unique_sequences | count of sequences | exact sequences first added by this source; the reservoir size is their sum |
| duplicate_records | count of sequences | input_records minus unique_sequences |
| source_seconds | s | builder time for the source |
| builder_threads | count of threads | lake-builder thread setting |

## campi_comparator_records.tsv
sha256: cf32390aae243b5c1cf35cbca83e41210163698dce213cf22b0b517e419721c9
rows: 49108 (one row per FASTA record; 29,673 reference and 19,435 MG)
produced by: campi_comparator_fastas

| column | unit | meaning |
|---|---|---|
| database | id | reference or MG (supplied CAMPI comparator FASTAs, decoys removed) |
| record_index | index | 0-based position of the record in the FASTA |
| accession | id | first token of the FASTA header |
| length | residues | sequence length |
| sequence_sha256 | id | SHA-256 of the upper-case sequence |

## union_contributions.tsv
sha256: 5eaa7eb32db23fa8306754567589d3f9550b8a1a5feb5650ae8f3604039dfeb9
rows: 1042 (one row per selected FASTA and roster: 58 + 306 + 331 + 347)
produced by: cohort_selected_fasta_unions

| column | unit | meaning |
|---|---|---|
| roster | id | HSTOOL58, PREDICT306, CAPSCAN331 (v28d) or CAPSCAN347 (v28f) |
| order | index | position when FASTAs are read in acquisition-name order |
| acquisition | id | acquisition (raw-file stem) whose selected FASTA this is |
| records | count of sequences | records in the selected FASTA |
| new_unique_sequences | count of sequences | exact sequences (I/L distinct) not in any earlier FASTA of the roster; the union size is the sum |
| fasta_sha256 | id | sha256 of the selected FASTA file |

## crossfit_lda_fallback_347.tsv
sha256: 3376a7236aa524713cecb61c31f724b4308e9420e8cd54e5d71f4bf6d00408d9
rows: 591 (one row per matched acquisition)
produced by: A95_spectrum_crossfit (rerun with CapScan caches pooled over 347 acquisitions); rows 578-590 are the 13 CapScan acquisitions of the 295 matched benchmark, none with a fallback (pipelines/capscan_mg_benchmark_295)

| column | unit | meaning |
|---|---|---|
| i | index | roster index (same as crossfit_index in source_data/Fig4/crossfit_per_acquisition.tsv) |
| cohort | id | PREDICT or CAPSCAN |
| lda_fallback | id | arms whose pooled discriminant fit was singular and fell back to Sage's heuristic score (empty for none) |

## split_spectra_mg_per_acquisition.tsv
sha256: b56f3ebeb17fd2a8f4eea1cc93617d30056d9b8d9c67c631f2cf693466c04660
rows: 578 (one row per matched acquisition)
produced by: metagenome_split_spectra_comparison

Equal to the Source Data package sheet `Supplementary_Figure_06/metagenome_same_spectra/per_acquisition.tsv`
without the acquisition file name and with `group` blank for PREDICT (a specimen-visit code).
Suffix `_5` is the 5-50-residue view, `_9` the 9-50-residue view.

| column | unit | meaning |
|---|---|---|
| i | index | roster index |
| cohort | id | PREDICT or CAPSCAN |
| group | id | CapScan sampling group (stool, duodenum, jejunum, ileum, colon); blank for PREDICT |
| subject_id | id | participant: PREDICT pseudonym (PREDICT_P###), CapScan participant number, NOID_ prefix when none is recorded |
| odd_mzml_identical | boolean | odd-scan mzML identical to the one of the split-half searches |
| mg_odd_5, mg_odd_9 | count of peptides | accepted peptides, matched MG database, odd-scan spectra |
| mg_full_5, mg_full_9 | count of peptides | accepted peptides, matched MG database, all spectra (rerun) |
| mg_full_equals_frozen_5, _9 | boolean | rerun equals the original MG search |
| frozen_mg_5, _9; frozen_fl_5, _9 | count of peptides | original full-spectra MG and FastaLake searches |
| fl_odd_full_5, _9 | count of peptides | FastaLake, all-prediction selection, odd-scan spectra |
| fl_odd_self_5, _9 | count of peptides | FastaLake, odd-scan (same-half) selection |
| fl_odd_disjoint_5, _9 | count of peptides | FastaLake, even-scan (disjoint) selection |
| gain_frozen_fullspectra_5, _9 | % | 100 x (frozen_fl / frozen_mg - 1) |
| gain_odd_full_5, _9; gain_odd_self_5, _9; gain_odd_disjoint_5, _9 | % | 100 x (FastaLake / mg_odd - 1) for the three selections |

Accepted peptides: rank-one targets at Sage peptide q <= 0.01, modifications removed, I/L equivalent.

## split_spectra_mg_summary_original.tsv
sha256: 71c6dd9c41861ece03bdbfbef1bcf88ef663c3e042207620bd06bdc9d1de8ba3
rows: 40 (one row per stratum, length view and comparison)
produced by: metagenome_split_spectra_comparison (byte-identical copy of the package sheet `summary.tsv`)

| column | unit | meaning |
|---|---|---|
| stratum | id | PREDICT, CAPSCAN_all, CAPSCAN_stool, CAPSCAN_intestinal |
| min_length | residues | 5 or 9 (upper bound 50) |
| comparison | id | gain column, or disjoint_minus_full_selection_gain_points |
| n_acquisitions, n_participants | count | acquisitions and participants in the stratum |
| median_gain_pct | % (points for the contrast) | median over acquisitions, rounded to 2 decimals |
| ci95_low, ci95_high | % | percentile participant-bootstrap interval, rounded to 2 decimals |
| n_fl_gt_mg, n_fl_lt_mg, n_equal | count of acquisitions | sign of the per-acquisition value |

## predict_reconversion.tsv
sha256: 71e50432555c8288966ef6a0a63c3318ff97d2cc59ceb081d60affb2e70c909b
rows: 306 (one row per PREDICT acquisition, public-sequence search)
produced by: sage_input_completeness

| column | unit | meaning |
|---|---|---|
| acquisition | id | acquisition (raw-file stem) |
| searched_reconverted_mzml | boolean | the searched mzML was reconverted from the RAW file |
| searched_mzml_closing_tag | boolean | the searched mzML ends with its closing tag |
| nul_in_tail | boolean | NUL bytes at the end of the searched mzML |
| declared_spectra | count of spectra | spectrumList count of the searched mzML |
| max_scan | scan number | highest scan scored by Sage |
| max_scan_fraction | dimensionless | max_scan / declared_spectra |
| searched_mzml_software | id | conversion software recorded in the mzML |
| original_mzml_bytes | bytes | size of the original mzML in the deposit folder |
| original_mzml_ending | id | COMPLETE or TRUNCATED (closing tag present in the last 200 bytes) |

## compute_nodes.tsv
sha256: d52abc4938734fe212f5a4b766604626abca0667034714e93c5c70200b5df8b2
rows: 8 (one row per node)
produced by: common_resource_accounting (scontrol, 2026-10-06)

| column | unit | meaning |
|---|---|---|
| node | id | Slurm node name |
| role | id | which job ran on it |
| CPUTot | count of logical CPUs | node CPUs |
| RealMemory_MB | MB | node memory as configured in Slurm |
| AvailableFeatures | id | node features (genoa = AMD EPYC Genoa) |
| Gres | id | generic resources (GPUs) |
| queried_utc | UTC time | when the record was read |

## crossfit_interpreter.tsv
sha256: ee621676006531d4cb5c6f1b81f124ae5bab80d7c6bc8a2aafe1e6376c0723f9
rows: 1
produced by: A95_spectrum_crossfit

| column | unit | meaning |
|---|---|---|
| launcher | id | launcher that names the interpreter |
| python | version | Python version of that interpreter |
| numpy | version | NumPy version of that interpreter |
