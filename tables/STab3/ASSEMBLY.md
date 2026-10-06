# Assembly: Supplementary Table 3

The table is a workbook, not a composite figure. It is in v28d as published and needs no
change for v28f.

1. Statistics. Input: controlled `PREDICT_metadata_2026-08-03.tsv` and
   `acquisition_accounting.tsv`. Tool: `LAB_ROOT/.../iteration1/fix_rifaximin/rifaximin_sensitivity.py`
   (here `tables/STab3/panels.py:stat_definitions`), Slurm job 6048159. Output:
   `paired_endpoints_by_scenario.tsv`, `paired_composition_by_scenario.tsv`, `counts.json`
   (sha256 in PROVENANCE.md).
2. Endpoints table for the workbook. Input: `paired_endpoints_by_scenario.tsv`. Tool: not
   recorded (OPEN_ITEMS.md, UNTRACED). Changes: columns renamed
   (`n_patients_increase` to `n_patients_on_greater_than_off`, `q_two_clinical_endpoints`
   to `q_BH_two_endpoints_within_definition`, `acq_off`/`acq_on` to
   `n_acquisitions_off`/`n_acquisitions_on`), `definition` and `units` columns added,
   metric names written in words, p and q values written with fewer digits. Output:
   `SupplementaryTable6_rifaximin_sensitivity.tsv` (d4273f1a...), here
   `source_data/STab3/endpoints_by_definition.tsv`.
3. Workbook. Input: step 2 output and `paired_composition_by_scenario.tsv`. Tool:
   `LAB_ROOT/.../final/figure_sd/source_data/build/specials.py::st6` via `build_package.py`
   (here `tables/STab3/panels.py:build_workbook`). Changes: Read_me sheet written, numeric
   strings converted to numbers, creator set to FastaLake; the package step neutralised
   document properties. Built as Supplementary Table 6 and renumbered to 3 in
   `selection.py`. Output: `Supplementary_Tables/Supplementary_Table_3.xlsx`.

Published workbook: `expected/Supplementary_Table_3.xlsx`, sha256
5d9815c575a8cc15f9f840994b2bf255fc95f1c0b3ddbb4aaf218c49a39a0f3d. The notebook's rebuild
(`build/STab3/Supplementary_Table_3.xlsx`) differs in bytes (zip timestamps) and is
identical in every cell (`checks/workbook.json`). There is no PNG or PDF composite.
