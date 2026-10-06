# Provenance: Supplementary Figure 13

`LAB_ROOT/.../` abbreviates
`LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/`.
`SF13/` abbreviates `LAB_ROOT/.../final/suppfigs_n15/SF13/`.
`C2/` abbreviates `LAB_ROOT/.../coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6/`.

## Panels

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a-h, unredacted (v28f, 347 CapScan acquisitions) | `SF13/build_sf13.py` (sha256 8aebb13a...), run by `SF13/run.sh` (sha256 1450f668...). It holds verbatim the filters of `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/code/review_histidine_predict.py` l.16 (copy `SF13/orig/`, sha256 22f10374...) and `plot_histidine_review.py` l.40 (6063eee3...), and the HISTIDINE_CONTEXT_REPLICATION block, `dots`, `title`, `check`, `save` and rcParams of `LAB_ROOT/00_PACKAGE_2026-09-17/.../supplements/revision_2026-09-24/build_revision.py` l.13-14, 254-275 (copy `SF13/orig/build_revision.py`, e636bd85...) | `SF13/new/HISTIDINE_CONTEXT_REPLICATION.pdf` (= `expected/HISTIDINE_CONTEXT_REPLICATION.pdf`) | 5e1983adc9aa79c778d14111a28634db8993eeb6e67e448b06f8a54bbedf1972 | 6068766 (srun on the host named in `SF13/run.log`; start and end times match the output files) | CPython with matplotlib 3.10.8, pandas 2.2.2 (`run.log`) |
| a-h, redacted PDF and 300 dpi PNG (v28f) | `SF13/redact_sf13.py` (bf149838...), the logic of `LAB_ROOT/00_PACKAGE_2026-09-17/.../codex_review_2026-09-25/prepare_supplement_figures.py` (copy `SF13/orig/`, 5c52b4b8...) with the title box found by text search | `SF13/new/Supplementary_Figure_13.pdf`, `SF13/new/Supplementary_Figure_13.png` (= `expected/SFig13_v28f.pdf`, `expected/SFig13_v28f.png`) | 25974f9e4c8758f225e94d779126110214512959a3dbfe225f260593edf9fa8d; 4c26ae6bb7b1e216dbb547dc689a2817ad493fc48b807098d4c9af14e4b0dd0e | 6068766 | PyMuPDF (version not recorded in `run.log`) |
| v28d figure (331 acquisitions, for reference) | `build_revision.py` then `prepare_supplement_figures.py` (originals above) | v28d `word/media/image21.png`; reproduced by `SF13/repro/Supplementary_Figure_13.png` | d0ac4534bb5b2d02a22ff4d36fbdb122fae8e3ace08d9abc4f82c3e50730b072 | 6068766 (reproduction) | as above |
| statistics, CapScan motif family (347) | `C2/A28_fixed347/frozen/collect_revision4.py` and `histidine_helpers.py` (identical to `LAB_ROOT/MM_review_analysis/A28_4_7_histidine_inference/frozen/`), launcher `C2/jobs/A28_fixed347.sbatch` | `C2/A28_fixed347/results_revision4/full_CAPSCAN/complete_statistical_families.tsv`, `participant_motif_differences.tsv.gz` | 3c24a42eec725e0b43bdeecca3928f971925c25395dde7478d92a1db9e1c28d7; df84cc2720b0ed63448048ef0676ad33c44ef9fe46eb05414a1d98a3ecb5ae77 | 6064031_3 (completed; the earlier 6058108_3 ran out of memory and its partial output was moved aside, `C2/JOBS.txt`) | CPython 3.10.9, numpy 1.26.4, pandas 2.3.3 (anaconda 2023.03 interpreter named in the launcher, with the user-site pandas) |
| statistics and panel g,h input, PREDICT | `LAB_ROOT/MM_review_analysis/code/collect_repaired_predict_20260915.py` (96289d6e...) with the A28 frozen code, launcher `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28.sbatch` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/acquisition_enrichment.tsv.gz`, `complete_statistical_families.tsv` | 8e49c01bac8e2af1f8f2570c4718fc4bf7dea5d648af1ff62f6f5d88cc14c1ab; 85c774cb604a1d9e6c9526b04750b10d6601a07855fe1552357f66be016d28d7 | 5937597 | as above |
| v30 composite, house style (all panels) | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig13/build_sf13_style.py` (sha256 3ce1554d...) with `sf_common.py` (sha256 9a1c9929...), run by `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/run_sf456.sh` | `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/v30/figures_style/SFig13/SupplementaryFigure13.pdf` and `.png` | pdf 6b650c4281252cb75b1c234c70917f49af0d9d90c775b41248f107152ad5daa6; png 890c191b2649592435d7f9e28c1d191d44cc0d4d1b19e182cf0f4884a848d41e | not recorded (build.log has no job id) | the repository environment `fastalake-paper_v28_work/env` (run_sf456.sh): CPython 3.12.4, matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3, PyMuPDF 1.27.1 |

v30 result: `panels.py:figure_v30` in the repository environment gives a PNG byte-identical to the v30 image (`checks/SFig13_v30.json`), with every plotted value identical to the producer drawing. The v30 build also wrote `values_check.json` and `pdf_audit.json` (fonts Söhne Buch and Halbfett only, text ink #4e4f51 or white, clip check passed); they are not copied here.

Repository result (notebook run on Slurm, 4 cores, 127 s wall time, 0.72 GB peak RSS; CPython
3.12.4, matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3, PyMuPDF 1.27.1):

- `HISTIDINE_CONTEXT_REPLICATION.pdf` and `SFig13_v28f.pdf`: pixel-identical at 300 dpi to the
  producer PDFs (`checks/HISTIDINE_CONTEXT_REPLICATION.json`, `checks/SFig13_v28f_pdf.json`, 0
  differing pixels); the files differ only in PDF metadata.
- `SFig13_v28f.png`: byte-identical to the producer PNG (`checks/SFig13_v28f.json`).
- Panels a-h cut from both PNGs with the same box: 0 differing pixels each (`checks/panel_<x>.json`).
- Display table rebuilt from the family table: byte-identical to the producer's Source Data.
  The all-tests export written to `build/` is byte-identical to the producer's
  `CapScan_histidine_all_tests.tsv` (2adeeaab...).
- Motif family recomputed from participant values: 32,400 of 32,400 tests, n and P equal
  exactly, mean, q and interval bounds within 1.4 x 10^-16.
- PREDICT source family recomputed from the pseudonymised acquisition table: 1,560 of 1,560 tests,
  n equal, mean, P, q and interval bounds within 1.1 x 10^-16.

The producer's own reproduction (`SF13/CHANGES_SF13.md`, `SF13/run.log`) showed that the same
code with the 331-acquisition A28 input (sha256 90debaad...) gives the v28d image byte for byte
and the delivered 331 Source Data byte for byte.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| v28d Supplementary Figure 13 | `word/media/image21.png`, sha256 d0ac4534... | BYTE_MATCH; MATCHING_REPRODUCER; producing script `SF13/redact_sf13.py:16` via `SF13/run.sh`; selected origin `SF13/repro/Supplementary_Figure_13.png`; limit: "original historical renderer not independently located" | the historical renderer is `build_revision.py` (block l.254-275) followed by `prepare_supplement_figures.py`; in the release candidate the image was `Supplementary_Figure_09.png`, written by `prepare_supplement_figures.py`, which Codex lists among the byte-identical files | yes, with a more specific original writer | `SF13/repro` reproduces it byte for byte from the original code copied verbatim |
| v28f candidate PNG | REBUILT_CANDIDATES.tsv, sha256 4c26ae6b... | REBUILT_CANDIDATE, write route `SF13/redact_sf13.py:16`, WRITE_ROUTE_INSPECTED_NOT_RERUN | `redact_sf13.py` via `run.sh`, job 6068766; reproduced byte-identically here | yes | expected: it is the v28f replacement |
| v28f candidate PDF | REBUILT_CANDIDATES.tsv, sha256 25974f9e... | REBUILT_CANDIDATE, write route `SF13/redact_sf13.py:16` | as above; reproduced pixel-identically | yes | |

Codex number recheck (`out_recheck/REPORT.md`, section "Supplementary Fig. 13 a-f") recomputed the
motif P values exactly and q within 2.22 x 10^-16, with 25,151 significant tests, 598 displayed
dots and the same per-panel dot counts as here (109, 107, 104, 107, 86, 85). It agrees with this
unit. Codex did not recompute the PREDICT g,h family; this unit did.

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| `source_data/Fig6/capscan_histidine_statistical_families.tsv` (owned by Fig6, read here) | `C2/A28_fixed347/results_revision4/full_CAPSCAN/complete_statistical_families.tsv` | 3c24a42eec725e0b43bdeecca3928f971925c25395dde7478d92a1db9e1c28d7 | byte copy (Fig6 `prepare_source_data.py`) |
| `source_data/SFig13/capscan_histidine_motif_display.tsv` | the table above, rows selected | 3c24a42e... | `panels.derive_display_table`, run by `prepare_source_data.py`; byte-identical to `SF13/source_data/CapScan_histidine_motif_display.tsv` (269fa990...) |
| `source_data/SFig13/capscan_participant_motif_differences.tsv.gz` | `C2/A28_fixed347/results_revision4/full_CAPSCAN/participant_motif_differences.tsv.gz` | df84cc2720b0ed63448048ef0676ad33c44ef9fe46eb05414a1d98a3ecb5ae77 | byte copy |
| `source_data/SFig13/predict_histidine_acquisition_enrichment.tsv.gz` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/acquisition_enrichment.tsv.gz` (= package `supplements/revision_2026-09-24/source_data/PREDICT_histidine/acquisition_enrichment.tsv.gz`, the file the producer read) | 8e49c01bac8e2af1f8f2570c4718fc4bf7dea5d648af1ff62f6f5d88cc14c1ab | `prepare_source_data.py:pseudonymise_predict` (identifier columns replaced as text; numbers unchanged) |
| `source_data/SFig13/predict_histidine_statistical_families.tsv` | `.../full_PREDICT/complete_statistical_families.tsv` (= package `PREDICT_histidine/complete_statistical_families.tsv`) | 85c774cb604a1d9e6c9526b04750b10d6601a07855fe1552357f66be016d28d7 | byte copy |

Comparison with the Source Data package (`LAB_ROOT/.../final/figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Figure_13/`)
and `SOURCE_DATA_MAP_reanchored.tsv` (rows SF13-af, SF13-fam, SF13-gh, RES-93):

- `CapScan_histidine_motif_display.tsv` (9034b795...) and `CapScan_histidine_all_tests.tsv`
  (467cbbcc...) in the package are the 331-acquisition versions. They do not agree with v28f: 598
  rather than 597 displayed dots, 25,151 rather than 25,052 significant motif tests, cell means
  changed by up to 0.28 percentage points in proximal and distal fluid. The v28f replacements are
  `SF13/source_data/` (269fa990..., 2adeeaab...); the first is `source_data/SFig13/capscan_histidine_motif_display.tsv`,
  the second is regenerated by the notebook into `build/` from the Fig6 family table.
- `PREDICT_histidine_all_tests.tsv` (1a25fa54...): same 1,572 rows, columns and labels as
  `predict_histidine_statistical_families.tsv`, written with fewer float digits; numeric values
  agree within 4.5 x 10^-16 (checked on Slurm with the repository environment).
- `PREDICT_histidine_participant_enrichment.tsv` (6278de0c..., 11 MB, patients `PREDICT_P###`): the
  producer's patient-mean table for all routes, origins and residues (107,280 rows). The notebook
  computes the same table from the acquisition table (`panels.predict_patient_means`, 107,280
  rows). Within every route x side x origin x residue group the sorted values agree within
  1.1 x 10^-16, but the package numbers patients in a different order: none of the 93
  `PREDICT_P###` labels points to the same patient as the repository's `P###` label (matched on the
  de novo C-terminal histidine values, each a unique match). The repository's labels are the ones
  Figure 6 uses: the H rows here agree with `source_data/Fig6/predict_histidine_patient_enrichment.tsv`
  within 10^-12 on all 1,116 rows. See `OPEN_ITEMS.md`.
- Map rows SF13-af, SF13-fam and SF13-gh point to the package files above; RES-93 is LOCATED only.
  This unit recomputes RES-93: 93 of 93 patients positive at C-terminal boundaries for both routes.
