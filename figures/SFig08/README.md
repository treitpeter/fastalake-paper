# Supplementary Figure 8

Clinical database comparators and evidence ablation. Panel a shows whether the PREDICT
peptide gain of public-sequence FastaLake over the specimen metagenome depends on how the
metagenome database is built: it adds a curated specimen MG arm (FastaLake selection
applied to the specimen metagenome) and a full FastaLake arm that includes the cohort
metagenomes. Panels b-d show how much of a fixed FastaLake reporting dictionary depends on
peptides that the matched metagenome search did not accept. Panel a counts peptides of
9-50 residues. The v30 figure is regenerated from `source_data/SFig08/` and is
pixel-identical to the v30 render (`checks/`).

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | accepted 9-50-residue peptides per acquisition in four PREDICT arms; bars are medians | v30: curated specimen MG a median 1.5% fewer than specimen MG (fewer in 177 of 296); public FastaLake +10.6% over curated, higher in 292 of 296 (123); legend 296 acquisitions (441) | 296 acquisitions per arm (PREDICT) | none; median of per-acquisition gains, sign counts | Panel a | PREDICT_curation_counts_by_length.tsv (PREDICT_curation_counts.tsv for the 5-50 control) | predict_curated_mg_comparator |
| b | LFQ intensity removed (%) when FastaLake-only features are dropped | legend (441); no value printed | 296 PREDICT, 282 CapScan acquisitions | none; medians | Panel b | clinical_fixed_dictionary_ablation.tsv | fixed_database_evidence_ablation |
| c | KOs losing all support (%) | legend (441); no value printed | 296 PREDICT, 282 CapScan acquisitions | none; medians | Panel c | clinical_fixed_dictionary_ablation.tsv | fixed_database_evidence_ablation |
| d | Bray-Curtis distance, normalized all-evidence vs shared-only KO profile | legend (441); no value printed | 296 PREDICT, 282 CapScan acquisitions | none; medians | Panel d | clinical_fixed_dictionary_ablation.tsv | fixed_database_evidence_ablation |

`STATS.tsv` has no rows: the figure and the sentences it supports report medians and
counts only, with no test or interval.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `Q_CUTOFF` | 0.01 (Sage peptide q, per search) | production acceptance threshold (Methods, v28d line 229) | none for this figure |
| `LENGTH_WINDOW` | 5-50 residues | accepted set from which the 9-50 counts are taken; the 5-50 column must equal the shipped table | none |
| `LENGTH_WINDOW_V30` | 9-50 residues | v30 legend a: "Accepted peptides of 9-50 residues" | 5-50 counts in `PREDICT_curation_counts.tsv` |
| gain estimator | median over acquisitions of 100 x (test/reference - 1) | Methods, v28d line 229 | ratio of cohort medians is a different quantity and is not used |
| `N_PREDICT`, `N_CAPSCAN` | 296, 282 | acquisitions with a specimen-matched metagenome database | CapScan participant 15 has a stool metagenome only (v28d line 129) |
| `PRODUCER_RC`, `COL`, `FIGSIZE` | as in the producer | used for the hidden producer render that `render_v30` checks against | none |
| `TITLE_BBOX`, `COMPOSITE_DPI` | producer-recorded title rectangle; 300 dpi | v28d title removal (`redact_title`, no longer used for v30); 300 dpi output | none |
| `V30_W_MM`, `V30_H_MM` | 160 x 134 mm | v30 placed size (the v29 empty title band removed) | none |
| `ARM_COLOUR`, `COHORT_COLOUR` | FastaLake arms fl_style FL blue, metagenome arms METAG orange; b-d both cohorts FL blue | v30 colours by meaning, as main Fig. 4 | producer `COL` per column |
| `V30_TICK_PT`, `V30_LAB_PT`, `V30_DOT_S`, `V30_DOT_A` | 6.5 pt, 7.0 pt, 9, 0.6 | v30 type and marker sizes | producer defaults |

## Known limitations

- Acquisitions are the unit. Patients (PREDICT) and participants (CapScan) contribute
  several acquisitions; the medians are not patient-balanced.
- Panels b-d are conditional on the FastaLake dictionary. They do not reconstruct a
  metagenome-only workflow, and absence from the metagenome identification set does not
  mean absence from its database.
- The curated specimen MG arm was searched without LFQ (Methods, v28d line 229).
- The v30 figure is restyled (Söhne, colours by meaning, no panel titles) but keeps the
  internal arm names on the panel a axis (FL cohortMG, FL public, MG, MG own curated) and
  the cohort order CapScan before PREDICT in b-d (`OPEN_ITEMS.md`).
- Panel a at 9-50 residues reverses the v28d curated-MG result (5-50: +9.1%, higher in 286 of
  296; 9-50: a median 1.5% fewer, lower in 177 of 296).
- The searches and the ablation that produced the two source tables are not re-run here;
  their code is referenced in `UPSTREAM.json`.
