# Provenance: Supplementary Figure 1

`LAB_ROOT/.../` abbreviates
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/`.

## Panels

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a (v28f, 726 acquisitions) | `LAB_ROOT/.../final/suppfigs_n15/SF01/build_sf1.py` (sha256 4a12cb41...; the SF1 block of `LAB_ROOT/review/v15_integration/build_evidence.py`, sha256 1e6881b6..., lines 8-27, with input/output paths as arguments and n= tick labels derived from row counts), run by `LAB_ROOT/.../final/suppfigs_n15/SF01/run.sh` (sha256 6e705ef7...) | `LAB_ROOT/.../final/suppfigs_n15/SF01/new/Database_size_context.pdf` (= `SupplementaryFigure1_n726.pdf`) | 95c8a147dcff61366ded17d9b434e65451eb776e58d41d388b51209210bbb0be | not recorded (run.log has no job id) | CPython 3.12.12, matplotlib 3.10.8, numpy 2.1.1, pandas 2.2.2 |
| composite raster (v28f) | `LAB_ROOT/.../final/suppfigs_n15/SF01/render.py` (sha256 23b89705...; raster step of `LAB_ROOT/.../iteration3/source_data/figures/reexport.py`, sha256 cab43851...) | `LAB_ROOT/.../final/suppfigs_n15/SF01/new/SupplementaryFigure1_n726_400dpi.png` | 29930b6da0b412a612fc66d06d4f8afdd22a70801e0fe5ac315b9c64650963c8 | not recorded | CPython 3.12.4, PyMuPDF 1.27.1, Pillow 12.0.0 |
| v28d figure (710 acquisitions, for reference) | `LAB_ROOT/review/v15_integration/build_evidence.py` then `LAB_ROOT/.../iteration3/source_data/figures/reexport.py` | `LAB_ROOT/review/v15_integration/evidence/Database_size_context.pdf` (sha256 2f70eb54...), rasterised to `LAB_ROOT/.../iteration3/source_data/figures/out/SupplementaryFigure1_v21_reexport.png` | 5aca1814afd4de89fb6fb98dcd8e1040c5fb9a69b63e328c8caf37bd9e0feb3b (= v28d `word/media/image7.png`) | not recorded | as above |

Repository result: `panels.py:panel_a` with the repository environment (CPython 3.12.4,
matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3) gives a PDF that is pixel-identical at 300 dpi to
the producer PDF (`checks/panel_a.json`, 0 differing pixels; the files differ only in PDF
metadata), and `panels.py:render_composite` gives a PNG byte-identical to the producer's 400 dpi
image (`checks/SFig01_v28f.json`, 0 differing pixels, same sha256). The producer itself had
shown that `build_sf1.py` with the old 710-row table reproduces the v28d image byte for byte
(`LAB_ROOT/.../final/suppfigs_n15/SF01/run.log`). Notebook run on Slurm: 4 s wall time, 0.40 GB
peak RSS, 4 cores requested.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| v28d Supplementary Figure 1 | `word/media/image7.png`, sha256 5aca1814... | BYTE_MATCH to `LAB_ROOT/.../iteration3/source_data/figures/out/SupplementaryFigure1_v21_reexport.png`; write site `reexport.py:126`; WRITE_SITE_HASH_RECEIPT (REEXPORT_VALIDATION.json) | same raster write site; the PDF it rasterised (`ready_integration_2026-09-28/.../Database_size_context.pdf`, sha256 2f70eb54...) is byte-identical to the `build_evidence.py` output in `LAB_ROOT/review/v15_integration/evidence/`, so the original producer is `build_evidence.py` lines 8-27 | yes | Codex did not rerun the producer; this unit and the producer's own repro both reproduce the image |
| v28f candidate PNG | REBUILT_CANDIDATES.tsv, sha256 29930b6d... | REBUILT_CANDIDATE, write route `SF01/run.sh`, WRITE_ROUTE_INSPECTED_NOT_RERUN, differs from all v28d media | `render.py` via `run.sh`; reproduced byte-identically by this repository | yes | expected: it is the v28f replacement |
| v28f candidate PDF | REBUILT_CANDIDATES.tsv, sha256 95c8a147... | REBUILT_CANDIDATE, write route `SF01/run.sh:13` | `build_sf1.py` via `run.sh`; reproduced pixel-identically | yes | |

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| `source_data/Fig1/acquisition_sizes_726.tsv` (shared with Figure 1; not duplicated here) | `LAB_ROOT/.../final/fig1bd_n15/acquisition_sizes_726.tsv` | fd9ade3bb4ebd2ee6c129f69819b9697cf20abd7c10c77a12644bef25396e386 | copied unchanged |
| `source_data/SFig01/human_swissprot_canonical_2026_03.fasta.gz` | `LAB_ROOT/MM_review_analysis/A01_4_1_database_sizes/inputs/human_swissprot_canonical_archived_2026_03.fasta` | b205ed9970fee62413e2f192c6fb82daab86e4a5897ef521f907ecbfcd7b3e63 | `gzip -9 -n`; `panels.py:stat_human_reference` checks the uncompressed sha256 |

Comparison with the Source Data package (`FastaLake_Source_Data/`): the map rows SF1 and
RES-08a-d point to `Review_v16/Database_sizes/database_sizes_710.tsv` (sha256 53d90a32..., 710
rows, columns cohort, sample, entries, reservoir_entries). All 710 of its rows are in the
726-row table with equal `entries` and `reservoir_entries`; the 726-row table adds 16 CapScan
rows (participant 15). `Figure1/acquisition_sizes.tsv` (sha256 d944c2d3..., 710 rows) agrees in
the same way. The CapScan rows of the 726-row table are exactly the 347 acquisitions of the
fixed CapScan roster. The package has no 726-row table and no `Supplementary_Figure_01/` folder
(see `OPEN_ITEMS.md`).

## v30 update (2026-10-06)

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a v30 | `LAB_ROOT/.../final/v30/figures_style/SFig01/build_sfig01_style.py` (sha256 126acf893b5e5fe409f2e022d88b6879001c5933540cf386ab823fcdf7c24611) with `LAB_ROOT/.../final/v30/figures_style/_shared_sf01_03/restyle_common.py` (0c7ced18eb696d77db2639e5121fb66399debe60684097d15cd59cff816addb1) | `LAB_ROOT/.../final/v30/figures_style/SFig01/SupplementaryFigure1.pdf` and `.png`; receipts `values_check.json`, `pdf_audit.json` | PDF 70152bd0628621386aeb1c2fdb44f6c89a4fdcb560da9455f93d4539ff9dbc8c; PNG 8978bee36972f810a7fd8416b268c768ad1b643058597867e328a980c6e7168d | not recorded | paper environment, matplotlib 3.10.7, PyMuPDF 1.27.1 |

The redraw and helper code were moved into `panels.py` unchanged apart from repository-relative
paths and `fl_style` from `style/`. Notebook run on Slurm (srun, 7 s): v30 PNG byte-identical.
Source data and numbers are unchanged; the v30 legend (v30f) prints 726 and 347, as queued for v28f.
