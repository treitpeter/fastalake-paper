# Source data: Supplementary Figure 8

`PREDICT_curation_counts.tsv` and `clinical_fixed_dictionary_ablation.tsv` are the exact files
the v28d producer read (byte-identical copies). `PREDICT_curation_counts_by_length.tsv` is the
9-50-residue recount the v30 panel a is drawn from.
They contain no clinical covariates. Samples are identified by their LC-MS/MS
acquisition file name; the run number in a PREDICT name links to a patient only through
a controlled file that is not in this repository.

## PREDICT_curation_counts.tsv
sha256: 33502dc224471a7e13c77302e65f3ac610640e19417381f0f06260fba6120661
rows: 1,184 (one row per acquisition and database arm; 296 acquisitions x 4 arms)
produced by: predict_curated_mg_comparator

Identical to Source Data package sheet `Supplementary_Figure_08/PREDICT_curation_counts.tsv`.

| column | unit | meaning |
|---|---|---|
| sample | id | PREDICT LC-MS/MS acquisition |
| arm | id | `FL_cohortMG` full FastaLake including cohort metagenomes; `MG` specimen-matched metagenome plus human Swiss-Prot (specimen MG); `FL_public` public-sequence FastaLake; `MG_own_curated` FastaLake selection applied to the specimen metagenome (curated specimen MG) |
| peptides | count of peptides | distinct accepted peptides: rank-1 target PSMs with Sage peptide q <= 0.01, 5-50 residues, modifications removed, I and L merged |

## PREDICT_curation_counts_by_length.tsv
sha256: 01b7f7c288fdb0dc799ab00c299a7d80b70886f37e7b47fc191ad5e0c82ecc3a
rows: 1,184 (296 acquisitions x 4 arms)
produced by: predict_curated_mg_comparator (9-50 recount, `recount_curation.py`)

| column | unit | meaning |
|---|---|---|
| sample | id | PREDICT LC-MS/MS acquisition |
| arm | id | as in `PREDICT_curation_counts.tsv` |
| peptides_5_50 | count of peptides | distinct accepted peptides of 5-50 residues (equals `peptides` in `PREDICT_curation_counts.tsv`, all rows) |
| peptides_9_50 | count of peptides | the same accepted set restricted to 9-50 residues; plotted in panel a (v30) |

## clinical_fixed_dictionary_ablation.tsv
sha256: e8f214c89a6670df364906ea34cf653932388d33e35a95342f2c521cccf70893
rows: 578 (one row per acquisition; 296 PREDICT, 282 CapScan)
produced by: fixed_database_evidence_ablation

Same values as Source Data package sheet `Supplementary_Figure_08/clinical_fixed_dictionary_ablation.tsv`;
that sheet was re-written by pandas and differs in the last binary digit of
`intensity_removed_percent` and `KO_support_lost_percent` in 29 rows (largest relative
difference 1.9e-16). This file is the one the figure was drawn from.

"Shared" features have a peptide accepted by the same-acquisition matched metagenome
search; "FastaLake-only" features do not. Group intensity is split equally across the
group's KO terms; intensity of groups without a KO stays in the total.

| column | unit | meaning |
|---|---|---|
| sample | id | LC-MS/MS acquisition |
| all_LFQ_intensity | LFQ intensity (arbitrary units) | summed positive feature intensity, all features |
| shared_LFQ_intensity | LFQ intensity (arbitrary units) | summed intensity of shared features |
| FL_only_LFQ_intensity | LFQ intensity (arbitrary units) | summed intensity of FastaLake-only features |
| fraction_intensity_removed | dimensionless | FL_only_LFQ_intensity / all_LFQ_intensity |
| all_KOs | count of KOs | KOs with any support |
| shared_KOs | count of KOs | KOs with support from shared features |
| FL_only_supported_KOs | count of KOs | KOs supported only by FastaLake-only features |
| normalized_KO_profile_BrayCurtis_all_vs_shared | dimensionless | Bray-Curtis distance between the all-evidence and shared-only KO profiles, each normalized to sum 1 |
| FL_only | count of features | quantified features classed FastaLake-only |
| outside_local_representative | count of features | features whose group representative is not in the acquisition's own FASTA |
| shared | count of features | quantified features classed shared |
| cohort | id | `PREDICT` or `CapScan` |
| intensity_removed_percent | % | 100 x fraction_intensity_removed (panel b) |
| KO_support_lost_percent | % | 100 x (1 - shared_KOs / all_KOs) (panel c) |
| BrayCurtis | dimensionless | copy of normalized_KO_profile_BrayCurtis_all_vs_shared (panel d) |
