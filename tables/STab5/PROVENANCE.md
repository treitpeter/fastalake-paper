# Provenance: Supplementary Table 5

## Panels

| panel | producer script (LAB_ROOT/...) | producer output (LAB_ROOT/...) | sha256 of producer output | Slurm job id(s) | interpreter and key package versions |
|---|---|---|---|---|---|
| tests, eligibility, summary | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/individual_functions.py (sha256 cc110778f31b9768e82a2dda533fa1deddb1c39ec07e7571732fd1f3ff8f7105; byte-identical copy run from .../claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/) | .../iteration1/fix_eggnog_native/downstream/verified50/B/individual_refit/crossomics_same_function_correlations.tsv; eligibility.tsv; same_function_summary.tsv | bbab94eb3dcbc52e24b841b070316cc5c5614080bcec0138a112659f9c90c18f; ccba7b357e04b450e61d586a8bf8934a719a96b0162a2324c9c9cd941d8d2b15; a7f59dbe3703c0f3d9be55d6d451aca737a813bbbfe04e78977ef713bd106a82 | 6048399 (refit receipt COMPLETE.json); profiles 6048301 | CPython 3.12.4, scipy 1.16.3, statsmodels 0.14.5, pandas 2.3.3, numpy 2.2.6 (the producer interpreter has these versions today; versions at run time not recorded in the receipt) |
| workbook (producer) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration2/cr003_followup/scripts/tables_cr003.py, variant verified50 (sha256 7453d63a10a4cb6f53f52d850c14fd51bcb009f4f263d5316b57b6eb43439f94) | .../iteration2/cr003_followup/tables/verified50/Supplementary_Table_4.xlsx (numbered 4 before the v21 renumbering) | c7dd4195c2404825cb89e2378e00959885cc467e67e87cc48b479c2899bae65e | not recorded | openpyxl 3.1.5 |
| workbook (as published) | LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/figure_sd/source_data/build/build_package.py with specials.py (st4, elig, elig_sum) (sha256 730ded5c10f3fc6964d554e249478711d45586470696b77bdf053166d0e74a04; 4b5d700bd04089af4682f851a310f77b882a5d1ed166be19703ecfa9191a36d4) | .../final/figure_sd/source_data/package/FastaLake_Source_Data/Supplementary_Tables/Supplementary_Table_5.xlsx | cef0a5eb4c7774120fc05e021660210a25b4220101c7b4fc72595ca4b63c5641 | build.sbatch (job id not recorded) | openpyxl 3.1.5 |
| support files | same build | .../FastaLake_Source_Data/Supplementary_Tables/Table5_support/eligibility_matched_correlations.tsv; common_eligibility_summary.tsv | d1ad7b4296d71f20fc160a568eaec576631021f87ce4afca0b98f6351fa719aa; cbbfa654762645f458be66a90768d6af2371caab597fff3f4a7d77e42050db00 | as above | as above |

Results of this unit (`checks/table.json`, `checks/recompute.json`): the rebuilt workbook
equals the published one in all 1,351,138 cells (sheet names and order, shape, values; no
tolerance used); both support files are byte-identical; the recomputed 30,581 tests,
184,386 eligibility rows and nine-family summary equal the produced tables (keys, patient
counts and eligibility flags identical; maximum absolute difference 1.0e-16 for rho and
1.1e-16 for P and q; the recomputed test TSV is byte-identical to the produced one).

The workbook bytes differ from the published file only by the timestamps that openpyxl
writes; this is why the comparison is by cell value.

## Codex census reconciliation

| panel or image | Codex CENSUS.tsv row (image id) | Codex attribution | this unit's attribution | agree? | note |
|---|---|---|---|---|---|
| Supplementary Table 5 | none | not covered (the census covers embedded figure images only) | producer chain above | NA | no image in this unit |

Codex number recheck (`out_recheck/REPORT.md`): no row for the numbers of v28d lines 141
or 829 that this unit covers. Codex statistics register: rows S093-S096 cover this table;
the gaps they list (family sizes, n per test, implementation, software) are filled in
`STATS.tsv` from the code and listed as TEXT-GAP in `OPEN_ITEMS.md`.

The Source Data map row (`SOURCE_DATA_MAP_reanchored.tsv`) and the earlier inventory name
`Biology_refit/individual_refit/same_function_summary.tsv` for line 141; that package file
is byte-identical to `source_data/STab5/same_function_summary.tsv` (sha256 a7f59dbe...).

The earlier paper repository carries a different Supplementary_Table_5 workbook (v12, old
allocation, 17,914 tests). It was not used.

## Source data

| source_data file | derived from (LAB_ROOT/... or Source Data package sheet) | sha256 of origin | transformation (code path) |
|---|---|---|---|
| crossomics_same_function_correlations.tsv | LAB_ROOT/.../biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native/downstream/verified50/B/individual_refit/crossomics_same_function_correlations.tsv; equals package Biology_refit/individual_refit/ copy and sheet All_30581_tests | bbab94eb3dcbc52e24b841b070316cc5c5614080bcec0138a112659f9c90c18f | byte copy (tables/STab5/prepare_source_data.py) |
| same_function_summary.tsv | .../verified50/B/individual_refit/same_function_summary.tsv; equals package copy and sheet Nine_family_summary | a7f59dbe3703c0f3d9be55d6d451aca737a813bbbfe04e78977ef713bd106a82 | byte copy |
| eligibility.tsv | .../verified50/B/individual_refit/eligibility.tsv; equals package copy and sheet All_eligibility | ccba7b357e04b450e61d586a8bf8934a719a96b0162a2324c9c9cd941d8d2b15 | byte copy |
| previous_allocation_same_function_correlations.tsv | LAB_ROOT/.../biology_refit/individual_baseline/crossomics_same_function_correlations.tsv (the file tables_cr003.py reads) | 3cf4e8e0b22812d156d201ec926c95e53023a7467653363b53ce3dfc657ef8a8 | byte copy |
| profiles_{OG,KO,Pathway}_metaG.tsv.gz, profiles_{OG,KO,Pathway}_metaT.tsv.gz | LAB_ROOT/analyses/enterotypes_2026-09-07/metaG_OG_TPM.tsv (c39df737891e118b92c3251aca87be433fae886fe16ceca2ad4e8eb55ebe1f73), metaT_OG_TPM.tsv (351a6cbeead4b21b7b5aa40e4c27918170f8a6093b2fa76aab8dc3083f792266), PREDICT_metaG_KO.tsv (077b3f48920f8041e579acaf80d0ff285de4dc4188e8161593ebaf00ebcef034), PREDICT_metaT_KO.tsv (4cb6c0000a46daefeb7073944377d5779efc085b6b5d617851247e9fe7015a51), PREDICT_metaG_Pathway.tsv (cc692926acf327b0ccc9f9e71489d4a13a2138432e4d4cd6c1b6d4cab89d3d52), PREDICT_metaT_Pathway.tsv (d8d267bc41a0610247e4c8fbb46e5c2fdf4dc131a31072cdfc306ccfc5fedb85) | as listed | rows of the 77 patients, columns of the layer-pair vocabularies, non-zero values in long format, patients pseudonymised (prepare_source_data.py) |
| profiles_{OG,KO,Pathway}_microbial_metaP.tsv.gz | LAB_ROOT/.../iteration1/fix_eggnog_native/profiles/verified50/PREDICT_noclust/{OG,KO,Pathway}_other_unresolved.npz and _axes.json (KO npz dda86f8256cb1325a414961317adff27fa1c886becf560a004b5acfc37f41179; OG 2e9334518e4cdcbf90348e7a4a21d3b1aeb39f26771b2ce273618b865c1c8b26; Pathway 47ba8caa332e00aaab17c14a5367138bf1e13b769fe59d4ddfb42cc5ea367790) | as listed | acquisitions collapsed to patients by the original `collapse` (mean), then as above |
| profile_features.tsv, patients.tsv | the same inputs | NA | feature vocabulary per level and layer; the 77 pseudonyms |

Pseudonyms are the Source Data package codes (`PREDICT_Pnnn`, key
`LAB_ROOT/.../final/figure_sd/source_data/internal_only/PREDICT_PARTICIPANT_KEY.tsv`,
sha256 f26ee8475b34233dce0b3af3b5365302f32069dc24e59a00f4b1ec2c4adc158e), mapped from
the study sample code through `LAB_ROOT/analyses/enterotypes_2026-09-07/PREDICT_fullstack_primary_samples.tsv`
(sha256 5e967b0aa182f2dfd6b2d902d2d1b393dd96caf43ef619575ee5c5b972aa09a7). Both are
controlled and are not in this repository. Preparation ran on Slurm (4 cores, 24 GB
requested, 40 s wall time).
