# Code availability

The statement below is the Code availability text of manuscript v30 (paragraph 160 of `visible_v30f.txt`), quoted unchanged. Items in square brackets are still to be filled by the authors.

> FastaLake is available under the MIT license at https://github.com/MannLabs/fasta-lake. The repository's Dockerfile builds a container image (ghcr.io/mannlabs/fasta-lake) with the compiled Rust binaries and the official Sage v0.14.7 release executable, so no local Rust or Sage installation is needed. Paths and function names cited in Methods and Supplementary Methods refer to version 1.1.0rc2 (commit 9ad09c9) [Zenodo DOI to be added], in which the lookup, selection-parsimony and study-grouping functions used here are unchanged apart from added input checks. The analysis and figure code for this study, including study-level protein grouping, functional and origin profiles, and clinical, cross-omics and regional analyses, is available under the MIT license at https://github.com/treitpeter/fastalake-paper [release, commit and Zenodo DOI to be added]. Its pipelines folder holds the code of each upstream analysis, its figures, tables and text_statistics folders regenerate every panel, table and number from Source Data, and docs/REVIEWER_GUIDE.md lists every number in the manuscript with the function and Source Data table that produce it. The original EDS feature-selection screen, for which no producing code is included, is supplied as Source Data. All searches except the joint healthy-stool searches used the unmodified Sage v0.14.7 release (commit 99407db, whose executable reports version 0.14.6). The joint healthy-stool searches used a build of the same source that keeps target and decoy q values separate when pooled LFQ q values are assigned to peaks. This defect was fixed upstream after the v0.14.7 release in lazear/sage commit 0e2c7a6, and the source, patch, lockfile and checksums of the corrected build are included in the analysis release.

## Versions cited

| what | version | where |
|---|---|---|
| FastaLake software | MannLabs/fasta-lake, version 1.1.0rc2 (`fasta_lake/__init__.py`), `main` at commit 9ad09c9 (merge of PR #3); no release tag yet | https://github.com/MannLabs/fasta-lake |
| Rust crates | fasta_extractor_v2, parsimony_engine, aggregation_engine, fasta_export 1.0.0; compiled with the pinned Rust toolchain 1.91.1 | fasta-lake `rust/`, `rust-toolchain.toml` |
| Container image | ghcr.io/mannlabs/fasta-lake (compiled Rust binaries and the Sage v0.14.7 release executable) | fasta-lake `Dockerfile` |
| Sage | v0.14.7 release (commit 99407db; the executable reports 0.14.6); joint healthy-stool searches used the same source with the LFQ q-value fix identical to lazear/sage commit 0e2c7a6 | `pipelines/healthy_stool_joint_search_q/sage_build/` (patch, Cargo.lock, producing source tarball, SHA256SUMS.txt) |
| This repository | treitpeter/fastalake-paper, analysis and figure code (MIT); release, commit and DOI to be added | this repository |

## Code pointers in Methods and Supplementary Methods

The v30 text cites software functions as `MannLabs/fasta-lake, <path>, <function>` and analyses as `fastalake-paper, pipelines/<id>`. Every pipeline id the v30 text names exists in `pipelines/` (checked 2026-10-06). The software functions cited are:

| step | path | function |
|---|---|---|
| reservoir build and SHA-256 accessions | rust/fasta_extractor_v2/src/bin/lake_builder.rs | `main` |
| reservoir lookup, I/L handling, eligibility and top-30% filter | rust/fasta_extractor_v2/src/main.rs | `main`, `normalize_il_seq`, `parse_prediction_csv`, `load_peptides_top_percent` |
| chunked mode | fasta_lake/chunked.py | `extract_chunked` |
| selection parsimony and FASTA export | rust/parsimony_engine/src/main.rs | `map_peptides_parallel`, `run_razor`, `write_fasta` |
| Sage configuration | fasta_lake/search.py | `_dda_config` |
| peptide acceptance, study grouping, quantification | fasta_lake/study.py | `_passing`, `greedy_set_cover`, `study_dictionary`, `group_searches`, `_quantify` |

Production study-level groups were made with the study.py run by `pipelines/campi_simplified_grouping/run_cohort.py`; the released `fasta_lake/study.py` reproduces the production healthy-stool grouping exactly (Supplementary Methods, paragraph 464).

## Basis in the repository

| statement | where |
|---|---|
| runs on Linux, CPython 3.12.4, pinned packages | `environment/README.md` |
| every number traceable, with the v30 paragraph | `docs/REVIEWER_GUIDE.md`, `NUMBERS.tsv` in each unit |
| upstream Slurm analyses with receipts | `pipelines/<id>/RECEIPTS.md` (paths relative to `LAB_ROOT`) |
| the original EDS feature-selection screen has no producing code and is supplied as Source Data | `figures/SFig10/OPEN_ITEMS.md`; the repaired screen is `pipelines/predict_eds_screen/` |
| release checks | `python paper.py check` |

## Open

- Release, commit and Zenodo DOI placeholders of both repositories.
- The v30 text cites `campi_second_engine` for the second-engine analysis; the clinical and healthy-stool MSFragger/Percolator code is in `pipelines/clinical_second_engine`, which the text does not name.
- The corrected Sage build is in `pipelines/healthy_stool_joint_search_q/sage_build/`. Its patch makes the same q-value fix as upstream lazear/sage 0e2c7a6 (q map keyed by precursor and decoy state in `picked_precursor`) but is not byte-identical: it also adds a regression test, a ranked-score audit and the `lfq_target_decoy_audit.tsv` export. The v30 sentence "The change is identical to upstream Sage commit 0e2c7a6" should say the q-value change is identical. The binary is not redistributed.
