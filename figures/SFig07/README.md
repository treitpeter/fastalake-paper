# Supplementary Figure 7

Supplementary Figure 7 asks what changes when every healthy-stool acquisition is searched
against one common study FASTA instead of its own selected FASTA: how many of the original
peptide identities are kept (a), how peptide and protein-group counts change per participant
in the two clinical cohorts (b), and whether microbial protein-group quantities agree better
between technical (c) and independently prepared (d) replicate pairs. The original producer
drew all four panels on one canvas; the notebook redraws that canvas with the same code and
cuts each panel from it.

## Panels

Claims are cited by v28d line; the printed values are those of the v30 text (v30f), which
moves panels a and b to peptides of 9-50 residues. The figure is the v30 house-style render.

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | original peptide identities (9-50 residues) retained per acquisition: joint search with study-level q (corrected and original Sage build), joint search with per-acquisition q, repeat search, separate common-FASTA searches | 115 (v30: 11.4%, 11.6%, about 94%, 34.2%; 10,493 to 13,149) | 58 acquisitions per column | none (medians, descriptive) | Panel a | `identity_comparison.tsv`, `identity_per_acquisition.tsv`, `completeness_by_length.tsv` | healthy_stool_common_fasta, healthy_stool_joint_search_q, healthy_stool_clinical_length_recount_9_50 |
| b | participant-level percentage change in accepted 9-50-residue peptides and directLFQ protein groups, common versus own FASTA | 437 (PREDICT 93, CapScan 15); 269 (v30: 4.4%, 2.4%, 3.6%) | 93 PREDICT patients, 15 CapScan participants | none | Panel b | `common_search_participants.tsv`, `identity_summary_by_length.tsv` | clinical_common_fasta, healthy_stool_clinical_length_recount_9_50 |
| c | median microbial protein-group pair CV, technical replicate pairs, own and matched rows | 769 (20.1 to 15.2%, r 0.83 to 0.87; study level 20.3 to 14.5%, 0.83 to 0.87) | 66 acquisition pairs (dependent) | none (medians, descriptive) | Panel c | `microbial_group_pairs.tsv` | healthy_stool_joint_search_q |
| d | as c, preparation replicate pairs | 769 (32.9 to 35.1%, r 0.76 to 0.72; study level 33.1 to 33.6%, 0.75 to 0.74) | 105 acquisition pairs (dependent) | none (medians, descriptive) | Panel d | `microbial_group_pairs.tsv` | healthy_stool_joint_search_q |

The leading-species agreement on line 769 is read from `technical_agreement_summary.tsv`
(see `NUMBERS.tsv`). Panels c and d and line 769 use protein-group intensities gated by
5-50-residue identities; they were not rebuilt at 9-50.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `MIN_LENGTH_PANEL_A` | 9 | v30: panels a, b and the Results numbers use peptides of 9-50 residues | 5-50 rows are in `identity_comparison.tsv` and the `*_by_length.tsv` tables (min_length 5) |
| `PEPTIDE_Q_CUTOFF` | 0.01 | rank-one target peptides at peptide q <= 0.01 (v28d line 765); applied upstream | none |
| `N_ACQUISITIONS` | 58 | the healthy-stool study (v28d lines 113, 115) | none |
| `TIERS` | TechnicalReplicates, WorkflowReplicates | panels c and d (v28d line 437) | other tiers are not paired in the upstream table |
| `TIER_PAIRS` | 66, 105 | all pairs of 12 and 15 acquisitions (v28d line 437) | none |
| `PAIR_CV_DDOF` | 0 | population SD / mean for a pair (v28d lines 437, 769); applied upstream | all-replicate sample-SD CV, upstream A21 `replicate_cv.tsv` |
| `PAIR_SUMMARY` | median | median over rows within a pair, then over pairs (v28d line 769) | none |
| `MICROBIAL_CLASS` | Microbial | panels c, d use microbial protein groups | Human, Mixed, Unknown classes exist upstream |
| `FIG_SIZE_IN`, `PNG_DPI` | 7.1 x 6.4 in, 300 dpi | producer canvas | none |
| `FONT_FILE` | soehne-buch.ttf | producer font file (TrueType, not the .otf) | `.otf` fallback renders different text pixels |
| `COL`, `RC` | producer palette and rcParams | used only for the hidden producer render that `render_v30` checks against | none |
| `V30_W_MM`, `V30_H_MM` | 160 x 144 mm | placed size in the v30 docx | none |
| `C_ORIG`, `C_STUDY`, `C_ACQ`, `C_JOINT_ORIG`, `C_UNION` | fl_style REF, FL, SKY, PURPLE, GREEN | v30 colours by meaning (reference grey; joint study q blue; joint acquisition q sky blue; original joint build purple; common-FASTA green) | producer palette `COL` |
| `V30_TICK_PT`, `V30_CAT_PT`, `V30_LAB_PT`, `V30_LEG_PT` | 6.5, 5.6, 7.0, 6.5 pt | v30 type sizes at placed size | none |
| `V30_DOT_S`, `V30_DOT_A` | 9, 0.6 | v30 dot size and alpha | producer 8, 0.48 |

## Known limitations

- The v30 figure is one canvas, so per-panel images in the notebook are crops at each axes'
  tight box; the pixel check is made on the whole composite (PNG and PDF), which is IDENTICAL
  to the v30 render.
- Panels c, d and the line 769 values use 5-50-residue identities to gate protein groups; the
  upstream comparison has no 9-50 counterpart.
- Panel b at 9-50 changes direction relative to v28d (common FASTA now raises peptide counts in
  68 of 93 PREDICT and 13 of 15 CapScan participants); `gained_peptides` and `lost_peptides` were
  not recounted at 9-50 and are empty.
- Per-pair CVs, completeness counts, distinct-peptide counts and leading-species agreement
  were computed upstream from search outputs that are not in this repository; this unit
  recomputes the medians over acquisitions, participants and pairs and reads the upstream counts.
- Pairs in panels c and d share acquisitions and are not independent; no test is made.
