# Open items: SFig03

- AUTHOR: The legend (v28d line 417) calls MG "the multi-omic co-assembly proteins"; the original CAMPI Methods built the SIHUMIx metagenome database from DNA only (Codex triage T01-T05). Suggested v28f legend wording: "the metagenome-derived proteins (MG)". The figure itself does not change. (legend, line 417)
- AUTHOR: Decision A2 (second engine for all cohorts; PREDICT 1.051, CapScan 1.364, healthy stool 0.998 in final/second_engine_all) is pending for the v28f text on line 89. This unit covers the CAMPI-only sentence as printed in v28d; if A2 extends the sentence, the new numbers belong to the unit that owns second_engine_all. (NUMBERS.tsv, line 89)
- TEXT-GAP: The v30 sentence for line 89, "recovering a median 85.3% of comparator peptide identities", does not say whether the median is per comparator or over both. All three versions print as 85.3% at 9-50 residues (REF 85.331%, MG 85.324%, both 85.327%); NUMBERS.tsv uses the median over both. (NUMBERS.tsv row 16, STATS.tsv rows 1-3)
- TEXT-GAP: Line 83 "3- to 17-fold fewer member-specific peptides" holds at 5-50 residues (3.37 to 16.75); at 9-50 it would be 3.38 to 22.6. The sentence does not name the length; suggested addition "(5–50 residues)". (NUMBERS.tsv rows 11-12)
- TEXT-GAP: Line 89 medians are over acquisitions (n = 15 per comparator) at 9-50 residues (v30) and peptide q <= 0.01 per acquisition and engine; the sentence does not state n, the length population or that the medians are of per-acquisition ratios. The legend (417) says "Dots are acquisitions". (STATS.tsv rows 1-3)
- AUTHOR: Taxon names in panels a and b use the reference-set name "Lactobacillus plantarum" (abbreviated L. plantarum); the text on line 83 uses "Lactiplantibacillus plantarum". The abbreviation is the same, so the figure need not change. (panels a, b)
- MISMATCH: NUMBERS.tsv keeps v28d line numbers and quotes as CONVENTIONS.md requires; v30 values are in the status column (line 89: 1.00 and 85.3% at 9-50 residues). The two legend rows for "c,d use peptides of 9–50 residues" quote v30 legend text that is not on v28d line 417, so a quote check against v28d fails for them. (NUMBERS.tsv rows 26-27)
- MISMATCH: The v30 Results (v30f) no longer print the line-83 values "twelve", "seven" and "3- to 17-fold"; their NUMBERS.tsv rows are kept against v28d. (NUMBERS.tsv rows 2-12)
- UNTRACED: no Slurm job ids are recorded for the len950 and v30 producer runs; outputs are identified by sha256 (PROVENANCE.md).

## Shared requests

- ZENODO: The accepted MSFragger peptide lists of the 45 searches (`results/<sample>__<arm>/accepted_5_50.txt` and `accepted_9_50.txt` of A07, 90 files, 31 MB in total) would let a reader recount panel d from peptides, as panels a-c are recounted here. Currently panel d is checked only against the per-search accepted counts. (panel d)
- PIPELINE: A14_campi_member_gap and campi_second_engine are referenced in UPSTREAM.json but not ported to pipelines/ (pipeline_dir null); their search code contains absolute lab paths. If ported, Fig2 and SFig04 (CAMPI) may share them. (UPSTREAM.json)
