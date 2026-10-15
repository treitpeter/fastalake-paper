# Controlled inputs

Patient-level clinical covariates of the PREDICT cohort (for example age, sex, MELD, CRP,
antibiotic exposure, visit dates, tube identifiers and raw patient identifiers) are not in
this repository. Code that needs them reads them from the directory named by the
environment variable `FASTALAKE_CONTROLLED_DIR`. When the variable is unset or a file is
missing, the code prints which file it needs and skips that step; it never fails silently
and never substitutes other data.

Public tables use pseudonymised patient identifiers only. The mapping from pseudonym to
clinical identifier is itself a controlled input.

Access: requests go through the PREDICT consortium data hub. Contact Jonel Trebicka
(contact details in the article), who initiates the required data transfer agreement.
Access is subject to the ethics approval and data protection terms of the PREDICT
study. The manuscript reports approval 35281-2/2017/EKU for the MUCOSA-PREDICT
sub-study at the University of Debrecen and written informed consent.

## Format

One row per file. sha256 is of the exact file the published results were computed from.

| file | columns | sha256 | needed by (notebook or pipeline, cell or function) | purpose |
|---|---|---|---|---|
| `paired_antibiotic_profiles.tsv` | patient_id, antibiotic_state, Human, Microbial_species, Microbial_above_species, Other_annotated_origins, Unresolved_origin, pg_parsimony, n_acquisitions | `b99dd07a5cf75150f7cb6520eedcc816cefcd6e562dea9b457c309dfa95f27ec` | figures/Fig4: figures/Fig4/panels.py:load_paired_controlled (Panels e, g; stat_eg_paired); text_statistics/methods: text_statistics/methods/stats.py:clinical_paired_tests (cell "Line 275") | per-patient off/on values for Fig. 4e,g and the paired statistics (copy of LAB_ROOT/.../manuscript_completion_2026-10-01/clinical_refit/paired_antibiotic_profiles.tsv) |
| `heatmap_85_patient_order.tsv` | microbiome_id, index, sample_id, patient_id, clinical_group, antibiotic_state, Human, pg_parsimony, dysbiosis_score | `43dc72c4b31e4bcf44c2c6bfbcc16a0a713ea47a8f377e86c87280fc5327b9ba` | figures/Fig4: figures/Fig4/panels.py:load_antibiotic_track (Panel h) | antibiotic track of Fig. 4h, aligned to heatmap_tracks.tsv by EDS score (copy of LAB_ROOT/.../verified50/B/eds_refit/input_data/heatmap_85_patient_order.tsv) |
| `predict_pseudonyms.tsv` | raw_patient_id, pseudonym | `cbc1ae4a900e541c3422f0c7ce91742e19a34679cc5b7fb24094573c51e57f27` | figures/Fig6: figures/Fig6/prepare_source_data.py (writes it; no notebook cell reads it) | map between PREDICT raw patient ids and the pseudonyms P001-P093 in `source_data/Fig6/predict_histidine_patient_enrichment.tsv`; other units that pseudonymise PREDICT patients (for example SFig13) should use the same map. |
| `patient_order_and_covariates.tsv` | antibiotic_state, clinical_group (rows matched by position and checked against dysbiosis_score, Human, pg_parsimony; the file also holds sample_id, microbiome_id, index, patient_id, which the code does not read) | `8a4844354b3dfcdc8087e19def097e25dcfebfc59fabf1f8bfbad918f06baa68` | figures/SFig10: SFig10.ipynb Panel a (panels.load_controlled), Panel b (panels.render), Statistics for panel a (panels.stat_a_adjustment_grid, panels.stat_a_joint_adjusted); tables/STab4: STab4.ipynb setup cell (panels.load_controlled), Checks and numbers (panels.build_sheets with ctrl), Statistics: the three pairs (panels.stat_joint_adjusted), Statistics: 192 adjustment estimates (panels.stat_adjustment_grid) | antibiotic track of panel a, and the antibiotic and clinical-group covariates of the adjusted correlations. Request this file, not the lab copy 6d8b8835... with raw identifiers; the notebook reproduces every value and both pixel checks with it (fix run 2026-10-06) |
| `SFig10_v30.png` | image (v30 figure, single panel, includes the per-patient antibiotic track) | `9b85a672a7c39d87371ec10e5d5a8c5e7361590e150a33303d8c13bebb5a6fd8` | figures/SFig10: SFig10.ipynb final cell (pixel check) | reference image of the PNG pixel check; replaces SFig10_v28f.png (7cc12746...) |
| `SFig10_v30.pdf` | vector image (same figure) | `7093d6897d946f253380df3a0f7431dc92bde1e85c3237783a1da427433d909b` | figures/SFig10: SFig10.ipynb final cell (pixel check) | reference image of the PDF pixel check; replaces SFig10_v28f.pdf (9e7c863d...) |
| `acquisition_accounting.tsv` | per-acquisition class fractions with patient_id, visit, antibiotic_state, clinical_group | `3880642bef93ccb390ff7be86b30a4aec715134b4f96c63ce0d5c6a4f5deaa23` | pipelines/predict_clinical_refit: clinical_refit.py main; tables/STab3: tables/STab3/STab3.ipynb, section "Statistics (controlled inputs)", panels.load_controlled / stat_definitions | paired antibiotic and clinical-group inference |
| `SData1_PREDICT_withheld_cells.tsv` | sheet, row_index, column_index, column, value (cells of august_visit_id, august_group, later_visit_number in PREDICT metadata, PREDICT mapping changes and PREDICT repeated tubes; recorded_visit_number in PREDICT stack eligibility; all cells of the PREDICT_ANTIBIOTIC42 rows of PREDICT membership (306 rows) and PREDICT patients (93 rows)); 5,254 cells | `1d52de25c42605a8d8ae90135050921e838131a8ac5cd312f5e2302bc869cfff` | tables/SData1: tables/SData1/panels.py:controlled_cells, notebook cells Table v28d and Numbers | complete v28d/v28f workbooks, the full-cell check against the published workbook, and the count "42 patients observed in both antibiotic states" (v28d line 501) |
| `PREDICT_metadata_2026-08-03.tsv` | sample_id, visit_rifaximin_for_he, the nine core antibiotic fields (visit_antibiotic_prophylaxis_for_sbp, visit_antibiotic_prophylaxis_for_variceal_bleeding, visit_antibiotic_prophylaxis_for_intervention, visit_other_type_of_antibiotic, visit_antibiotic_for_infection_nonprophylactic, visit_rifaximin, visit_norfloxacin, visit_cephalosporines, visit_penicillines), all visit_antibiotic_treatment_* fields except visit_antibiotic_treatment_other_specify | `fc5dbf07df2fd25c7def3e9d583b2745c9f237b12a6a9c6e78ecefb249e90238` | tables/STab3: tables/STab3/STab3.ipynb, section "Statistics (controlled inputs)", panels.load_controlled / stat_definitions | antibiotic state under the four exposure definitions |
| `legacy_259_heatmap_columns.tsv` | sample_id, patient_id, clinical_group, antibiotic_state, Human, pg_parsimony, microbiome_id, dysbiosis_score (259 rows; pseudonymised Source Data copy PKG/EDS_focused_20261003/source_data/) | `a5c5bf922943e00937d77bf3fbd8b5cf2973d68fdb1d846cc2ceeceecc134740` | tables/STab4: panels.load_controlled, panels.build_sheets (Observation_metadata, matrix headers) | complete published layout of the 259-observation sheets |
| `Supplementary_Table_4.xlsx` | the published workbook (contains the withheld columns) | `4fe89a55c7607b0e5932377cbd488c031bc1f0c7f99c255f6c9eb7d77747e403` | tables/STab4: Checks and numbers cell (panels.compare_workbooks) | cell-for-cell check of the complete rebuilt layout against the published file |
| `PREDICT_fullstack_primary_samples.tsv` | MicrobPred_ID, idpatient, visit_number, Tube_ID | `5e967b0aa182f2dfd6b2d902d2d1b393dd96caf43ef619575ee5c5b972aa09a7` | tables/STab5: tables/STab5/prepare_source_data.py:main (not the notebook); text_statistics/supplementary_methods_1_5: prepare_source_data.py:predict_distances (not the notebook) | maps the 77 study sample codes to patients for pseudonymisation of the shipped profiles |
| `PREDICT_PARTICIPANT_KEY.tsv` | study_code, package_code | `f26ee8475b34233dce0b3af3b5365302f32069dc24e59a00f4b1ec2c4adc158e` | tables/STab5: tables/STab5/prepare_source_data.py:main (not the notebook) | Source Data package pseudonyms PREDICT_Pnnn used in source_data/STab5/ |

## Withheld from public tables

- The PREDICT study visit code (column `group` of the Source Data package sheet `Figure4/benchmark_0_per_acquisition.tsv`: admission, readmission and scheduled-visit codes) is a clinical event field. It is removed from `source_data/SFig06/` (figures/SFig06 OPEN_ITEMS). Whether the package sheet is republished without it is an author decision.
- Supplementary Data 1 (`tables/SData1`) withholds the PREDICT visit, group and repeated-tube cells listed under `SData1_PREDICT_withheld_cells.tsv`.

## How code reads a controlled file

    import os
    from pathlib import Path
    CONTROLLED = os.environ.get("FASTALAKE_CONTROLLED_DIR")
    path = Path(CONTROLLED) / "<file>" if CONTROLLED else None
    if path is None or not path.is_file():
        print("SKIPPED: needs controlled input <file> (see controlled/CONTROLLED_INPUTS.md)")
    else:
        ...

## Published Fig. 4 composites

The published Fig. 4 composites show the per-patient antibiotic panels (e, g, h) and are held with the controlled inputs, not in the repository: `Fig4_expected/Fig4_v28f.{png,pdf}`, `Fig4_v30.{png,pdf}` and, for v33 (b,c at 295 CapScan acquisitions), `Fig4_v33.{png,pdf}` (sha256 of the PNG `3946a6c9ce18d7f8837faca3291a15179859618da4ac2943ec0bda2a60733974`), and, for v36 (b labels give the median paired change), `Fig4_v36.{png,pdf}` (PNG sha256 `f7c0eb116c5a93087ec0b22e124db248017447d50846fa220f2e2102f81a1b1b`, PDF sha256 `822ab39ec8a8d87bab763d1738fc108082ff9f97463fa3490a4cf776bec15ebd`; the PNG is the image embedded in `FastaLake_v36.docx`), plus `Fig4_artwork/panel_a.pdf`. `figures/Fig4/panels.py:check_against_v30` uses the newest of `Fig4_v36.png`, `Fig4_v33.png`, `Fig4_v30.png` that is present.
