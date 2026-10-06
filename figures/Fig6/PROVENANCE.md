# Provenance: Figure 6

`CR` below stands for `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03`
and `C2F6` for `CR/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15/fig6`.

## Panels

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| a-f (builder) | `CR/final/fig6_n15/code/build_figure6.py` (sha256 `eb1e5f03acca87ff30d023581fb8594c3cd940db20d94093822e74f9e7869f4f`), run by `code/run_all.sh` (`2a2749d9...4c26`) | `CR/final/fig6_n15/build_n15/panels/Figure6.pdf`, `.png` | pdf `5506d837465e9654abb9c82386281edac0fc5fe61544aea8909e0daae5d2b3ca`; png `f466a1e8e2886a352c97a6ef1b604441e06ddac0c2d82194ea606a1bfb4e7a25` | not recorded (log `CR/final/fig6_n15/logs/build_n15.log`) | CPython 3.10.9 (Anaconda 2023.03), matplotlib 3.10.8, numpy 1.26.4, pandas 2.3.3, PyMuPDF 1.27.1; style `LAB_ROOT/figure_style/fl_style.py` (`58b64c34...64c3`) with the lab styleguide `init_plotting` and Söhne Buch |
| composite | `CR/final/fig6_n15/code/postprocess_fig6.py` (`57ae5896...ffc237ff`), `relabel_regions.py` (`781650b8...8a8154`), final re-save in `run_all.sh` | `CR/final/fig6_n15/Figure6.pdf`, `.png` (= `build_n15/Figure6.*`) | pdf `362a5cd71c47c25c8abafbff99b5bd9ae7a492ec4d310794c7ebda6978540392`; png `c822e293b35061f00124aae4bf2114d668361d29990841895f699a2c2add0b70` | not recorded | as above |
| v30 composite (panel b key "Human-assigned", flattened) | `CR/final/v29/figures/_tools/flatten.py` (sha256 `d73f26fc54d59e42e3fa439b11582f1c4158ff9e309f5ec69bab26dda2ab8cd1`) with `--relabel "Host: capsule > stool in 15/15, q < 0.001" "Human-assigned: capsule > stool in 15/15, q < 0.001"` on `CR/final/fig6_n15/Figure6.pdf`; receipt `CR/final/v29/figures/Fig6/flatten_receipt.json` | `CR/final/v29/figures/Fig6/Figure6.pdf`, `.png` (the v30 manuscript figure) | pdf `07873927e9f96c65375376cbcdfe400f05eae21c44c661c0a04a56cd2c2ae26d`; png `ca90228f7902a7528a5abf340637a9bba6c47793bb3b767d498076c88f0913d1` | not recorded | PyMuPDF 1.27.1, numpy |
| a-c inputs | A24 `C2F6/A24_fixed347/frozen/search.py` | `C2F6/A24_fixed347/results/CAPSCAN/` | see Source data | 6058149 (array task 2) | CPython 3.10.9 (Anaconda 2023.03) |
| d-f inputs, CapScan | A28 `C2F6/A28_fixed347/frozen/collect_revision4.py`, `histidine_helpers.py`; selection `CR/final/fig6_n15/code/prep_histidine_inputs.py` (`28f47288...e454`) | `C2F6/A28_fixed347/results_revision4/full_CAPSCAN/`; `CR/final/fig6_n15/inputs_n15/histidine/` | see Source data | 6064031 (array task 3; first attempt 6058108 failed, output moved aside) | CPython 3.10.9 (Anaconda 2023.03) |
| d-f inputs, PREDICT | A28 PREDICT run; `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/code/review_histidine_predict.py` | `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/histidine/` | see Source data | not recorded | not recorded |

Repository re-run (this notebook, `paper.py build Fig6` on Slurm, 4 cores, 563 MB peak RSS,
3 min 14 s): CPython 3.12.4, matplotlib 3.10.7, numpy 2.2.6, pandas 2.3.3, PyMuPDF 1.27.1.
Every panel crop, the builder PDF and PNG and the v28f PDF and PNG are IDENTICAL to the
producer outputs (0 differing pixels, `checks/*.json`).

v30 re-run (2026-10-06, Slurm jobs 6075212 and 6075295 (final), 2 cores, `tools/run_notebooks.py figures/Fig6`,
160-216 s), superseded below: the relabelled, flattened composite `panels.relabel_and_flatten` wrote was
IDENTICAL to the v30 figure at 300 dpi (PNG byte-identical, sha256 `ca90228f...d1`;
flattened PDF 0 differing pixels at 300 dpi, file sha256 `5d71b351...1933` differs from
`07873927...6d` only in PDF metadata).

Key in plot code (2026-10-06, Slurm job 6075480, 4 cores, `paper.py build Fig6`, 273 s): `panel_b` writes
"Human-assigned"; `relabel_and_flatten` only flattens. v30 PNG: 5,884 differing pixels, v30 PDF: 6,831,
all inside the panel b key (x 1266-1882, y 535-589 px at 300 dpi); panels a, c-f IDENTICAL.
`checks/Fig6_v30_png.json`, `checks/Fig6_v30_pdf.json`.

The producer's own reproduction check (`CR/final/fig6_n15/REPRODUCTION_CHECK.json`)
covered only the n = 14 control build against `CR/final/fig_v27/Figure6.*`; the checks in
this unit are the first independent re-render of the n = 15 figure.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| v28d embedded Figure 6 | Figure 6, `word/media/image6.png`, sha256 `357bbd56...a365` | BYTE_MATCH; write site `LAB_ROOT/review/v15_integration/relabel_regions.py:38` with receipt `REGION_FIGURE_CHANGES.json`; origin `LAB_ROOT/review/v15_integration/region_figures/Figure6.png` | the n = 14 build; byte-identical to `CR/final/fig_v27/Figure6.png` and `CR/final/fig6_n15/build_old/Figure6.png`; not the v28f figure | yes | v28d image is superseded by the n = 15 composite |
| v28f composite PNG | REBUILT_CANDIDATES.tsv, `CR/final/fig6_n15/Figure6.png`, sha256 `c822e293...0b70` | write route `relabel_regions.py:38` inspected, not re-run; differs from all v28d media | producer chain `run_all.sh`; re-run here, pixel-identical | yes | this unit adds the re-run Codex did not do |
| v28f composite PDF | REBUILT_CANDIDATES.tsv, `CR/final/fig6_n15/Figure6.pdf`, sha256 `362a5cd7...0392` | as above | as above | yes | |

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| `capscan_terminal_participant_region.tsv` | `C2F6/A24_fixed347/results/CAPSCAN/per_participant_region_endpoints.tsv.gz` | `6d998e903173e85cff0166be599d05fc0005d1b462bf29101a8a93fd3f746dac` | rows `min_length == 9`, unit not `assigned_bare_sequence`, decompressed (`prepare_source_data.py`) |
| `capscan_residue_participant_region.tsv` | same file | same | rows `min_length == 9`, unit `assigned_bare_sequence` |
| `capscan_terminal_statistical_families.tsv` | `C2F6/A24_fixed347/results/CAPSCAN/complete_statistical_families.tsv` | `919360392724871c2cef8e99e5a0172cdb7d8ad8a364a476778aac7c8f10f70b` | copy |
| `capscan_acquisition_design.tsv` | `C2F6/A24_fixed347/results/CAPSCAN/acquisition_design.tsv` | `4638d7894cdd3423e33a7416bfe2de11d9e761e6bbde4749d29e0ed4a03c7fc4` | copy |
| `capscan_histidine_statistical_families.tsv` | `C2F6/A28_fixed347/results_revision4/full_CAPSCAN/complete_statistical_families.tsv` | `3c24a42eec725e0b43bdeecca3928f971925c25395dde7478d92a1db9e1c28d7` | copy |
| `capscan_histidine_participant_region_enrichment.tsv` | `C2F6/A28_fixed347/results_revision4/full_CAPSCAN/participant_region_enrichment.tsv` | `a830d42712865d23f6642060c782b2dd00b49db19b551dc3575cdb359f619722` | copy |
| `capscan_histidine_participant_reciprocal.tsv` | `C2F6/A28_fixed347/results_revision4/full_CAPSCAN/participant_reciprocal_controls.tsv` | `74613f6da70994b5307cf95660f40f1eb3a22cc9757175c48d5b28161980263c` | copy |
| `capscan_histidine_primary_source_enrichment.tsv` | `CR/final/fig6_n15/inputs_n15/histidine/CAPSCAN_primary_source_enrichment.tsv` | `aa6587cd2d51c19a1ed0d38d510708d420ade2552387435a4d7a31f432f9845b` | `panels.derive_histidine_inputs` from the A28 family table; byte-identical to the producer input |
| `capscan_histidine_primary_reciprocal_tests.tsv` | `CR/final/fig6_n15/inputs_n15/histidine/CAPSCAN_primary_reciprocal_tests.tsv` | `6b12f31b2e7360dd3675ab37c6e2e2c33689ce97c582fc84728d5ca356fe4109` | as above |
| `predict_histidine_primary_source_enrichment.tsv` | `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/histidine/PREDICT_primary_source_enrichment.tsv` | `8f2cb3185f336832a59d0ba0be519e76eda1134275ff8febb98c604c887de3cf` | copy |
| `predict_histidine_primary_reciprocal_tests.tsv` | `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/histidine/PREDICT_primary_reciprocal_tests.tsv` | `61e8490827b12dcd29720e3750706e4900184b2e88fda92ecdafd243db258ee8` | copy |
| `predict_histidine_statistical_families.tsv` | `LAB_ROOT/MM_review_analysis/recovery_2026-09-15/cohort_recovery/A28/results_revision4/full_PREDICT/complete_statistical_families.tsv` | `85c774cb604a1d9e6c9526b04750b10d6601a07855fe1552357f66be016d28d7` | copy |
| `predict_histidine_patient_enrichment.tsv` | `LAB_ROOT/MM_review_analysis/end_of_day_2026-09-15/histidine/PREDICT_patient_source_enrichment.tsv` | `d4f768b9a0a41b57a042e3e3888e7fe7615c733de086b28047ecd9635d783c27` | raw patient id replaced by pseudonym P001-P093 in ascending raw-id order, column renamed `patient`; map kept outside the repository |

### Comparison with the Source Data package (`CR/final/figure_sd/source_data/package/FastaLake_Source_Data/Figure6/`)

| this unit's table | package sheet | values agree? |
|---|---|---|
| `capscan_terminal_participant_region.tsv`, `capscan_residue_participant_region.tsv` | `per_participant_region_endpoints.tsv.gz` (`a89b485c...`) | no: the package sheet is the n = 14 / 331-acquisition build; this unit uses the n = 15 / 347 replacement |
| `capscan_terminal_statistical_families.tsv` | `CAPSCAN_A24_terminal_complete_statistical_families.tsv` (`05ab4715...`) | no: package is n = 14 (q 0.000571 for a, b; here 0.000290) |
| `capscan_histidine_statistical_families.tsv` | `CAPSCAN_A28_histidine_complete_statistical_families.tsv` (`90debaad...`) | no: package is n = 14 |
| `capscan_histidine_primary_source_enrichment.tsv`, `capscan_histidine_primary_reciprocal_tests.tsv` | `CAPSCAN_primary_source_enrichment.tsv` (`46cfdbfb...`), `CAPSCAN_primary_reciprocal_tests.tsv` (`d84b102d...`) | no: package is n = 14 |
| `capscan_histidine_participant_reciprocal.tsv` | `CAPSCAN_A28_histidine_matched_reciprocal_controls.tsv` (`9f3191dc...`) | not compared row by row (different table level; package is n = 14) |
| `predict_histidine_primary_source_enrichment.tsv` | `PREDICT_primary_source_enrichment.tsv` | yes, byte-identical (`8f2cb318...`) |
| `predict_histidine_primary_reciprocal_tests.tsv` | `PREDICT_primary_reciprocal_tests.tsv` | yes, byte-identical (`61e84908...`) |
| `predict_histidine_statistical_families.tsv` | `PREDICT_A28_histidine_complete_statistical_families.tsv` | yes, byte-identical (`85c774cb...`) |
| `predict_histidine_patient_enrichment.tsv` | none in `Figure6/` (map row RES-93 points to `Supplementary_Figure_13/PREDICT_histidine_participant_enrichment.tsv`) | not compared |

`SOURCE_DATA_MAP_reanchored.tsv` rows RES-89, RES-90, RES-91, RES-92, FIG6-f, FIG6-q and
FIG6-c point to the n = 14 package sheets and quote the v28d (n = 14) values; FIG6-q's
family size (4,536) is unchanged.

### Comparison with Codex number recheck (`CODEX_TASKS_2026-10-05/out_recheck/REPORT.md`)

Codex rows N0189-N0194 and N0201 agree with this unit on every v28f value (5.1-5.6%,
3.1%, 0.34%, fifteen, 15.2 and 17.6 pp, about 40% and 2.4%). Codex computed the histidine
shares as 39.864/40.290/24.352% from acquisition-level tables; this unit, using the
figure's own aggregation (mean over participants of the region-group mean), gets
39.873/40.297/24.370%. Both round to the same printed values (39.9%, 40.3%, 24.4%).
Codex's 39.880% and 2.386% for line 175 equal this unit's mean over regions 1-4 of the
non-specific route.

### Comparison with Codex statistics register (`CODEX_TASKS_2026-10-05/out_statsreport/STATISTICS_REGISTER.tsv`)

| Codex row | Codex value | this unit | agree? | evidence |
|---|---|---|---|---|
| S113 (host H-share, fluid vs stool) | BH family residue-composition, 13,608 tests | BH family `terminal_fraction__length_9`, 4,536 tests | no | in `source_data/Fig6/capscan_terminal_statistical_families.tsv` the test (route unspecific, unit canonical_peptide, origin Human, metric `H_fraction_nonKR`, contrast capsule minus stool, length 9) sits in `terminal_fraction__length_9`; `panels.stat_terminal_families` assigns it there; STATS.tsv row 4 |
| S114 (microbial H-share, fluid vs stool) | as S113 | as for S113, origin Microbial | no | STATS.tsv row 5 |

The 13,608-test family `all18_nonKR_composition__length_9` holds the panel c residue tests
(unit assigned_bare_sequence, metric `residue_H_fraction_nonKR`, regions 1-2 or 3-4 minus
stool; STATS.tsv rows 6-7), not the panel b tests. Codex also records n = 14 for S113 and
S114 (v28d build); the v28f build has n = 15.
