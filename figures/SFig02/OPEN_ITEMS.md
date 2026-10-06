# Open items: SFig02

- MISMATCH: v28d line 59 prints "86,000 query–seed evaluations" at 13 or more residues; the count is 85,835 (`stat_ab_lookup`, NUMBERS.tsv row for line 59). v28f should print 85,835, or keep 86,000 with "about". NUMBERS status `v28f_change:85,835`.
- UNTRACED: v28d line 561, "13,512 peptide–acquisition–seed observations through the mass route across 60 acquisitions and two seeds" belongs to Supplementary Table 6 and is not in this unit's source data; left to the STab6 unit (NUMBERS rows OPEN).
- UNTRACED: v28d line 665, "233,497,702-sequence reservoir" is a count of reservoir FASTA entries (stated in the A02 upstream report) that cannot be recomputed from repository tables; left to the unit that owns the reservoir size (NUMBERS row OPEN).
- AUTHOR: v28d line 679 says "19 settings: minimum evidence lengths of 6–12 at top 30%, and score fractions of 5–100% at minimum length nine". The 19 arms behind the 1,102 observations are 18 distinct settings plus the older 97-pool reference arm at nine residues and 30%, which is in Source Data but not drawn. The v28d text no longer says that the 1,102 records include this reference arm; v28f could add it.
- AUTHOR: whether Supplementary Fig. 2e,f is updated for v28f is undecided. The unit inventory reports that a lab status note (V27_STATUS, not re-read here) listed the native and shuffled-control bag reruns as due for v28; v28d text and legend use the 48-acquisition bounded comparison, and no newer e,f render was found. This unit reproduces the v28d e,f.
- TEXT-GAP: the quartile implementation for Supplementary Fig. 2c,d is not stated (Python `statistics.quantiles(n=4, method='exclusive')`; Codex register S012, software). Filled in STATS.tsv from the code.
- TEXT-GAP: the software that matched queries to the reservoir for Supplementary Fig. 2a,b and line 59 is not named in the sentence (the production Rust extractor; Codex register S131-S133, software). Filled in STATS.tsv.
- UNTRACED: the script that wrote the selection-sensitivity workbook (`SupplementaryTable_selection_sensitivity.xlsx`) was not identified. Its 1,102 rows equal the upstream sweep table `per_sample_by_arm.tsv` in every compared column (PROVENANCE.md, Source data).
- UNTRACED: the interpreters that ran `build_revision.py` and `integrate_v18.py` are not recorded; their PDFs name Matplotlib 3.10.8, the repository uses 3.10.7. The regenerated panels are pixel-identical, so this has no visible effect.
- PIXEL: panels a-d and panels e,f are each one canvas in the original producers, so they are regenerated and checked per canvas (`panel_abcd`, `panel_ef`, notebook sections "Panels a-d" and "Panels e,f"), not per letter. All four checks are IDENTICAL (0 differing pixels); the two manuscript images are byte-identical to v28d `image8.png` and `image9.png`.
- PIPELINE: database_selection_specificity: upstream scramble and lookup code not in `pipelines/` (pipeline_dir null).
- PIPELINE: healthy_stool_selection_sweep: upstream sweep code not in `pipelines/`; the id is provisional and chosen not to collide with ids used by other units.
- PIPELINE: bag_mass_lookup_search: upstream bag/mass search and fragment-audit code not in `pipelines/`.
- UNTRACED: no Slurm job id is recorded for the v30 house-style build (`build_sfig02_style.py`); outputs are identified by sha256 (ASSEMBLY.md, v30).
- AUTHOR: the v30 PNGs are not yet among the media of the v30f docx; placement in the next build is pending (ASSEMBLY.md, v30).
- MISMATCH: v30 text (v30f) still prints "86,000 query–seed evaluations"; recomputed 85,835 (NUMBERS.tsv row for line 59).

## Shared requests

- `source_data/SFig02/selection_length_lookup.tsv` is byte-identical to Source Data package `Figure1/length_by_arm.tsv`; if the Fig1 unit needs it, it should reference this file rather than copy it.
- The fragment audit contains one PREDICT acquisition, renamed `PREDICT_acq_01` in `source_data/SFig02/bags_fragment_audit_summary.tsv`. A repository-wide pseudonym scheme for PREDICT acquisitions would keep names consistent across units; the mapping to the raw file name is not in the repository.
