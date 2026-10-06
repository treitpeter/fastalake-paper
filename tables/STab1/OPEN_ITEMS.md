# Open items: STab1

- AUTHOR: The v28d legend at line 475 is the title only; for v28f it should be replaced by `expected/LEGEND.txt` (first sentence unchanged, adds row definition, thresholds, 80% bidirectional coverage, peptide acceptance and the KO recovery definition) (legend, line 475).
- AUTHOR: The Source Data package still carries the superseded 45-row, 9-sheet `Supplementary_Tables/Supplementary_Table_1.xlsx` (sha256 9731b751...ba61); swap in the 105-row workbook and add the TSV as `Table1_support/Supplementary_Table_1_105_rows.tsv`, keeping the eight `Table1_support/*.tsv` (PROVENANCE.md, Source data).
- AUTHOR: `SOURCE_DATA_MAP_reanchored.tsv` row ST1 still says "45 search rows" and RES-41 maps line 89 to the 97/99%-only `paired_count_ratios.tsv`; the 95/90/75/50% medians and the 96.4-99.0% KO range have no map row (NUMBERS.tsv rows for line 89).
- TEXT-GAP: Line 89 does not say that 96.4-99.0% is the range of the six per-threshold medians (median over 15 acquisitions of per-acquisition KO recovery, lowest at 50%, highest at 99%); a reader could take it as the range over acquisitions, which is wider (NUMBERS.tsv, line 89; STATS.tsv last row).
- TEXT-GAP: Line 89 reports medians without any dispersion; the code computes none. Per-acquisition values are in the table (STATS.tsv, uncertainty column).
- TEXT-GAP: Line 89 and the legend call the unit "acquisitions" without stating that, by their names, S11_Fraction1 to S11_Fraction4 are four fractions of one CAMPI sample (STATS.tsv, distinct_or_repeated).
- UNTRACED: The interpreter and package versions of the original table build (Slurm job 6065679) and of the independent recount (6064573-6064575) are not recorded; the table regenerates byte-identically with CPython 3.12.4, pandas 2.3.3, openpyxl 3.1.5 (PROVENANCE.md).
- UNTRACED: `LEGEND.txt` is listed as an output in the builder's docstring but the builder does not write it; how it was written is not recorded (ASSEMBLY.md step 2).
- PIPELINE: campi_candidate_clustering: the clustering, Sage searches and KO collection (`clustering_90_95_2026-09-29/prepare.py`, `frozen/selection.py`, `frozen/search.py`, `frozen/collect.py`, launchers) and the earlier 97/99% H04/H08 runs are not in `pipelines/`; their manifest and launchers hold absolute lab paths, so moving them needs path parametrisation by the orchestrator (UPSTREAM.json).
- PIPELINE: campi_clustering_recount: `clustering_check/recount.py` is not in `pipelines/` (absolute paths in its header) (UPSTREAM.json).
- ZENODO: The per-peptide inputs behind the counts (`fixed_peptide_KO_mapping.tsv.gz`, 0.8 MB; `KO_fates.tsv.gz`, 0.4 MB; `member_specific_peptides.tsv`) and the Sage result files are not in the repository; the two per-acquisition tables can be checked only against the independent recount here, not recounted from peptides.

## v30 reference (added 2026-10-06)

- UNTRACED: the v28d quote of NUMBERS.tsv rows 2 (15), 3 (0.9%), 4 (2.2%), 5 (2.5%), 6 (4.1%), 7 (7.2%), 8 (8.8%) was not located in the v30 text (`visible_v30f.txt`); the wording or the value changed or the sentence was removed. v30_status NOT_IN_V30; the v28d columns are unchanged.
