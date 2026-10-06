# Figure 3

Figure 3 uses the 58 healthy-stool acquisitions in five sampling series. It asks how much of an
acquisition's peptide identification is kept when its search database is selected from another
acquisition, what pooling the databases of 12 technical replicates adds, and how closely
protein-group and KEGG-ortholog (KO) intensities agree between acquisition pairs within each
series. Panels b-e and f are regenerated from `source_data/Fig3/` by the code of their original
producers and are pixel-identical to those producers' renders. Panel a is a schematic. In v30,
panels b-e use 9-50-residue peptides and the b legend reads "Same series" / "Different donor"; the
v30 composite is rebuilt by `panels.assemble_v30` and is pixel-identical to the v30 figure PNG.

## Panels

| panel | what it shows | claim (v30 paragraph; v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | sampling design of the five series | 58 acquisitions in five series (v30 48; v28d 99) | 58 acquisitions: 12, 15, 10, 11, 10 | artwork | artwork | artwork (`ARTWORK_NUMBERS.md`) | none |
| b | 9-50-residue peptide count with a same-series (T1) or different-donor (T4) database, % of own | 96.0% technical, 55.8% between individuals (v30 48; v28d 99 printed 96.4%, 55.1% at 5-50) | 58 target acquisitions, one T1 and one T4 search each | none (medians) | Panel b | transfer_pool_comparisons.tsv | healthy_stool_database_transfer |
| c | own-database 9-50-residue peptide identities recovered, % | 89.6% technical, 50.2% between individuals (v30 48; v28d 99 printed 89.2%, 49.1%) | as b | none (medians) | Panel c | transfer_pool_comparisons.tsv | healthy_stool_database_transfer |
| d | 9-50-residue peptide count with own, all-12 and other-11 databases, % of own | +10.0% in every acquisition (other 11); +10.7% (all 12) (v30 50; v28d 103 printed +9.1%, 4,450 gained, 1,589 lost) | 12 technical-replicate acquisitions | none (medians, count) | Panel d | transfer_pool_comparisons.tsv | healthy_stool_database_transfer |
| e | own-database identities recovered with the unions, % | 96.8% (other 11), 97.5% (all 12) retained (v30 50; v28d 103 printed 94.9%) | 12 acquisitions | none (median) | Panel e | transfer_pool_comparisons.tsv | healthy_stool_database_transfer |
| f | pairwise r and CV heatmaps, protein groups and KOs | CV 21% to 52%, r 0.82 to 0.48; KO r 0.90 and 0.81, above protein groups in every design (v30 54; v28d 109; unchanged) | 66, 105, 45, 55, 45 acquisition pairs (dependent) | none (medians) | Panel f | protein_group_pairwise.tsv, ko_pairwise.tsv, design.tsv | peptide_centric_downstream, healthy_stool_ko_precision_refit |

The donor-assignment count in Methods line 739 / v30 Supplementary Methods paragraph 379 (nine differ, one coincides) is computed in the
final cell from `transfer_pool_comparisons.tsv`.

The 60 numbers printed inside panel f (median r, CV and feature count for 20 class labels) are
checked against the producer's receipt (`expected/panel_f_producer_receipt.json`, 20 of 20
records equal, 10 of 10 heatmap matrices equal by sha256) and by the pixel check; they are not
manuscript-text numbers and are not repeated in `NUMBERS.tsv`.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `MIN_LENGTH` | 9 | v30 reports b-e at 9-50 residues: the accepted 5-50 set filtered to >= 9 residues, q not re-estimated (v30 legend paragraph 52, Supplementary Methods paragraph 379) | 5-50 residues: rows with `min_length` 5 in `transfer_pool_comparisons.tsv` (v28d values; Source Data) |
| `PEPTIDE_Q` | 0.01 | peptide q threshold applied upstream in A15 (Methods 739) | none in this figure |
| `TRANSFER_ARMS` | T1_within_group, T4_different_donor | same-series donor and inter-individual-series donor (Methods 739) | none |
| `POOL_ARMS` | own, pool_all12, pool_other11 | other-11 excludes the target's own database (Results 103) | all-12 shown alongside in d, e |
| `DESIGNS` | five healthy-stool series in figure order | design of the cohort (Results 99) | none |
| `PG_CLASSES` | host, microbial | group-level classes (Methods 247, 651) | `all` class in `protein_group_pairwise.tsv` (not plotted) |
| `KO_PROFILES` | new_all, new_other_unresolved | revised allocation, all origins and human-assigned excluded (legend 107, Methods 747) | previous allocation `old_all`, `old_microbial` in `ko_pairwise.tsv`; Source Data HSTOOL_KO_refit |
| `HEATMAP_SCALE` | 0 to 1 | producer setting; all plotted r lie in [0, 1] (asserted) | none |
| `HEATMAP_COLOURS` | white to #0072B2 | producer setting | none |
| `JITTER_SEED_BC`, `JITTER_SEED_DE` | 3, 9 (plus series or arm index) | producer's horizontal jitter of points; affects display only | none |
| `JITTER_HALFWIDTH_BC`, `JITTER_HALFWIDTH_DE` | 0.05, 0.08 | producer's jitter widths | none |
| `SIZE_BCDE_MM` | 180 x 97 mm canvas | producer setting (saved with a tight bounding box) | none |
| `SIZE_F_PT` | 518.74 x 211.87 pt | producer setting: page width and f band of the FAR21 composite | none |
| `V30_*` | band rows 333-1520 at 300 dpi; legend box and baselines; 8 / 6.5 / 7 pt | v30 composite step (`build_fig3_len950.py`, `relabel.py`) | none |
| `FONT_BCDE`, `FONT_F_REGULAR`, `FONT_F_BOLD` | soehne-buch.otf; soehne-buch.ttf, soehne-halbfett.ttf | the font files each producer loaded | without them text falls back to the default face and the pixel checks are skipped |

## Known limitations

- `NUMBERS.tsv` carries the v28d columns and the four v30 columns (`v30_para`, `v30_quote`,
  `v30_printed_value`, `v30_status`); v30-only numbers have `NA` in the v28d columns.
- Panel f is unchanged in v30 and stays at 5-50 residues (its producers have no 9-50 counterpart).

- Panels b-e were drawn on one canvas by the original producer and are regenerated together
  (`panel_bcde`). The composite step re-typesets all b-e text and relabels the b legend
  (`ASSEMBLY.md`), so the regenerated b-e match the producer render pixel for pixel but differ
  from the composite band in text (`checks/composite_band_bcde.*`). The full v30 composite
  rebuild (`checks/composite_v30.json`) is IDENTICAL.
- Panel b-e use the figure style as it stood on 2026-08-02 (0.6 pt lines, tight bounding box),
  moved into `panels.py` (`_ProducerStyle`), not the later `style/fl_style.py`.
- The 12 technical replicates are injections of one digest, and pairs in f share acquisitions.
  All statistics are descriptive; no test is made in this figure.
- The per-pair r and CV values are inputs here; they were computed upstream (A33, A90) and are
  not recomputed from intensity matrices in this repository.
- Panel a is artwork extracted from an earlier manuscript figure; no script that draws it exists.
