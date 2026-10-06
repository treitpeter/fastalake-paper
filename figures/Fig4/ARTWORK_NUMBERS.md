# Artwork numbers: Fig4

`artwork/panel_a.pdf` is the panel-a band (top 77 pt) of the v28f composite `expected/Fig4_v28f.pdf`, clipped as vector (ASSEMBLY.md, step A). Its typed text, extracted from the PDF: "Metagenomics", "Metatranscriptomics", "Metaproteomics", "Metabolomics", "Liver disease / ACLF", "Patient stool", "Matched sampling", "Companion dataset" (three times), "306 LC–MS acquisitions", "FastaLake · this study", "Metagenome benchmark", "296 matched acquisitions". Two of these carry numbers.

| artwork file | typed text | value | table it comes from (repo path, row, column) | status |
|---|---|---|---|---|
| artwork/panel_a.pdf | 306 LC–MS acquisitions | 306 | source_data/Fig4/predict_acquisition_roster.tsv, all rows (count of acquisition_id) | MATCH |
| artwork/panel_a.pdf | 296 matched acquisitions | 296 | source_data/Fig4/predict_acquisition_roster.tsv, rows with matched_metagenome = True (count); equal to the PREDICT rows of crossfit_per_acquisition.tsv and spectrum_audit_per_acquisition.tsv | MATCH |

Both counts are recomputed in the notebook (section "Panel a") and listed in NUMBERS.tsv (v28d lines 121 and 123).
