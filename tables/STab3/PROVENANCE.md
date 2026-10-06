# Provenance: Supplementary Table 3

`F` = `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01`
`R` = `F/biology_refit/claude_review_2026-10-03`

## Panels

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| statistics, all sheets | `R/iteration1/fix_rifaximin/rifaximin_sensitivity.py` (sha256 99889423a575c0852683c85f20e753ad9cfa95c869a26a31b12d03887a15ac7e), launcher `run.sbatch` | `R/iteration1/fix_rifaximin/out/paired_endpoints_by_scenario.tsv`; `.../out/paired_composition_by_scenario.tsv`; `.../out/counts.json` | f33383cad5fcbaca59aa1586de4435b225931c6d5129d50c12845965132e3f7d; b3dc72b682754e1a2c47f0ed5b193bcefb8ab173fedef1167c6f5abde371bf58; 447677655f48c01c25898b3f2fea6f95f425c54a6721847142c0ea89228aa2f9 | 6048159 (main, 7 s, MaxRSS 3.5 MB reported), 6048158 (negative control) | Python 3.12.4, numpy 2.2.6, pandas 2.3.3, scipy 1.16.3 |
| Endpoints sheet input | not found (see OPEN_ITEMS.md, UNTRACED) | `R/iteration1/fix_rifaximin/SupplementaryTable6_rifaximin_sensitivity.tsv` | d4273f1aefad9ca37bafdf8357c1cd67a923af47c87a1ea3f02e4520af9dbf3f | none recorded | none recorded |
| workbook | `R/final/figure_sd/source_data/build/specials.py::st6` (sha256 4b5d700bd04089af4682f851a310f77b882a5d1ed166be19703ecfa9191a36d4) called by `build_package.py` (sha256 730ded5c10f3fc6964d554e249478711d45586470696b77bdf053166d0e74a04); renumbered Supplementary Table 6 to 3 in `selection.py` | `R/final/figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Tables/Supplementary_Table_3.xlsx` (written 2026-10-04 06:00 CEST; row in `work/full/MANIFEST_build.tsv`) | 5d9815c575a8cc15f9f840994b2bf255fc95f1c0b3ddbb4aaf218c49a39a0f3d | not recorded in the build logs | openpyxl 3.1.5 |

Upstream inputs of the producer (controlled, not in this repository):
`LAB_ROOT/canonical/publication/metadata/PREDICT_metadata_2026-08-03.tsv` (fc5dbf07df2fd25c7def3e9d583b2745c9f237b12a6a9c6e78ecefb249e90238) and
`F/clinical_refit/acquisition_accounting.tsv` (3880642bef93ccb390ff7be86b30a4aec715134b4f96c63ce0d5c6a4f5deaa23), written by `F/clinical_refit.py`
(b496050501da5075efdf4aff75127bd357b2d4aceb277cb2209399f5e3735072, job 6042666), whose
antibiotic rule is that of `LAB_ROOT/FIGURATOR/Figure4/v5_2026-09-07/derive_clinical.py`
(924a04c9ee07e05166400b803c9d2f5a4bbb328ec362a5bbc2a635ead346c051) lines 21-27. `panels.derive_state` re-derives
the state from the metadata and asserts equality with the state in `acquisition_accounting.tsv` for all 306 acquisitions.

What was ported. `panels.py` holds `rifaximin_sensitivity.py` with paths made relative and
constants named; `derive_state`, `paired_inference`, `bh` and the scenario and count logic
are unchanged. Two additions: the sign-flip extreme count and draw count are returned
(for STATS.tsv), and the result tables are returned instead of written. Not ported: the
block that counts state changes among the 85 EDS display samples (it needs a third
controlled file and does not enter this table). `build_workbook` is `specials.py::st6`
with the Read_me text as a constant and the composition header restored from `state`.

Recomputation in this repository (Slurm, 4 cores, 4 GB, notebook run 2 s,
`python paper.py build STab3` with the controlled inputs):

- Endpoints and composition: bit-identical to the producer outputs (all 12 numeric
  columns, 8 rows; all 5 composition columns, 8 rows).
- Counts: identical to `counts.json` (all keys except EDS_85).
- Against the published Endpoints sheet: maximum relative difference 4.2 x 10^-13
  (the sheet's input TSV carries p and q values written with fewer digits).
- Definition A reproduces `primary_paired_statistics.tsv` (rtol 1e-10).
- Negative control (every HE flag inverted): A unchanged (42 patients), B collapses to 7
  patients, as in the original job 6048158.
- Rebuilt workbook against `expected/Supplementary_Table_3.xlsx`: IDENTICAL, 213 cells,
  3 sheets (`checks/workbook.json`).

Earlier independent check (not reused here): `R/iteration1/fix_rifaximin/independent_check.py`
(job 6048159) shares no code with the producer and gave identical n, means and CIs for
A-D; its P values from 10^6 sign flips with another seed agreed within Monte Carlo error.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| Supplementary Table 3 (workbook) | none | not covered: the census lists the 23 images embedded in the v28d docx; tables are not embedded | workbook from `specials.py::st6`, statistics from `rifaximin_sensitivity.py` | NA | the Codex number recheck (`out_recheck/REPORT.md`) has no row for this table |

Codex statistics register (`out_statsreport/STATISTICS_REGISTER.tsv`):

| register row | definition, endpoint | register gaps | filled here (STATS.tsv) |
|---|---|---|---|
| S055, S056 | A, human share; A, depth | software | software, statistic, seeds per endpoint (primary rows also belong to Fig4) |
| S057, S058 | B, human share; B, depth | software | software, statistic, seeds per endpoint |
| S059, S060 | C, human share; C, depth | effect size and CI, exact P or q, software | all from the code; values are in the table, not in the text |
| S061, S062 | D, human share; D, depth | effect size and CI, exact P or q, software | as above |

Register values that agree with the code: test, two-sided, unit patient, n (42, 40, 35, 38),
BH family of two endpoints within each definition, 200,000 sign flips with plus-one,
50,000 bootstrap resamples, seed 20260907. The register does not say that the depth
endpoint's sign flips use seed 20260908 (OPEN_ITEMS.md, TEXT-GAP).

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| `endpoints_by_definition.tsv` | `R/iteration1/fix_rifaximin/SupplementaryTable6_rifaximin_sensitivity.tsv`; package sheet `Supplementary_Table_3.xlsx/Endpoints` | d4273f1aefad9ca37bafdf8357c1cd67a923af47c87a1ea3f02e4520af9dbf3f | byte copy (`prepare_source_data.py`); rebuilt sheet equals the package sheet in every cell |
| `endpoints_producer_output.tsv` | `R/iteration1/fix_rifaximin/out/paired_endpoints_by_scenario.tsv` | f33383cad5fcbaca59aa1586de4435b225931c6d5129d50c12845965132e3f7d | byte copy; agrees with the package sheet to 4.2 x 10^-13 relative |
| `composition_by_definition.tsv` | `R/iteration1/fix_rifaximin/out/paired_composition_by_scenario.tsv`; package sheet `Composition` | b3dc72b682754e1a2c47f0ed5b193bcefb8ab173fedef1167c6f5abde371bf58 | header `antibiotic_state` written `state`; values unchanged; rebuilt sheet equals the package sheet |
| `primary_paired_statistics.tsv` | `F/clinical_refit/clinical_paired_statistics.tsv` | 9518b0117267eb3198b1aafe785dab47b53b020bd4deb0d9d2fb24f67ccc3bc7 | byte copy (same file as `source_data/Fig4/clinical_paired_statistics.tsv`) |
| `definition_counts.tsv` | `R/iteration1/fix_rifaximin/out/counts.json` | 447677655f48c01c25898b3f2fea6f95f425c54a6721847142c0ea89228aa2f9 | nested keys flattened with '.', `EDS_85` dropped |

`SOURCE_DATA_MAP_reanchored.tsv` rows ST3-a to ST3-l (40 patients; 11.4697; 5.74566;
17.6549; 0.000479998; -2378.06; -3651.47; -1140.04; 0.000539997; A equals primary; 35; 38)
all agree with the recomputed values in NUMBERS.tsv.
