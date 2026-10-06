# Source data: Supplementary Figure 1

The per-acquisition size table this figure reads is shared with Figure 1b and lives under
`source_data/Fig1/acquisition_sizes_726.tsv` (sha256
`fd9ade3bb4ebd2ee6c129f69819b9697cf20abd7c10c77a12644bef25396e386`, 726 rows, one row per
acquisition); its data dictionary is in `source_data/Fig1/README.md`. The notebook checks that
sha256 before use.

## human_swissprot_canonical_2026_03.fasta.gz
sha256: ec416a6d4316597a1ee37cc2a7314dee29809ae1c3443530d1339b91e40c0492 (gzip file); uncompressed sha256 b205ed9970fee62413e2f192c6fb82daab86e4a5897ef521f907ecbfcd7b3e63
rows: 20,431 FASTA entries (one entry per protein sequence)
produced by: cohort_database_sizes (input file of the size audit; UniProtKB/Swiss-Prot reviewed human canonical proteome, release 2026_03, isoforms excluded), compressed with `gzip -9 -n`

| column | unit | meaning |
|---|---|---|
| header line (`>sp|accession|entry_name description`) | id | UniProtKB/Swiss-Prot accession, entry name and description |
| sequence lines | amino-acid residues | protein sequence, 60 residues per line |

UniProt data are distributed under CC BY 4.0.
