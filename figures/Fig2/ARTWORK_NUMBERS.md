# Artwork numbers: Figure 2

Panel a contains one schematic, the SIHUMIx community artwork `artwork/bottle_61.pdf` (reused from
an earlier figure version; no generating script). Its text is "CAMPI SIHUMIx", "8-member gut model ·
15 acquisitions" and the eight species names. The second line was removed from the v28d figure by
`ASSEMBLY.md` step 6, so v28d shows no typed number in the artwork. The rows below trace the numbers
the artwork file itself carries. The database table in panel a is drawn by code from
`source_data/Fig2/metrics_per_acquisition.tsv` and is covered by the pixel check, not listed here.

| artwork file | typed text | value | table it comes from (repo path, row, column) | status |
|---|---|---|---|---|
| artwork/bottle_61.pdf | 8-member gut model | 8 | source_data/Fig2/reservoir_protein_presence_by_proteome.tsv, the 8 rows whose label does not end in _SHUFFLED, column label | MATCH |
| artwork/bottle_61.pdf | 15 acquisitions | 15 | source_data/Fig2/paired_reference_recovery.tsv, distinct values of column sample (15) | MATCH |
| artwork/bottle_61.pdf | species names (B. thetaiotaomicron ... L. plantarum) | 8 names | source_data/Fig2/member_summary.tsv, column member | MATCH |
