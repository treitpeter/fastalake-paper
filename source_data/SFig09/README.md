# Source data: Supplementary Figure 9

CapScan matched acquisitions: FastaLake (public-sequence reservoir, databases selected from
the scored spectra) against the specimen-matched metagenome (MG) search, with the size of
each MG database. Peptide counts are unique modification-stripped peptides of 5-50 residues
at Sage peptide q <= 0.01. MG database size is the number of predicted-protein records kept
by the assembly sanitizer (`n_kept`), without the added human Swiss-Prot entries.

Participant identifiers (`subject`, 1 to 15) are the CapScan study's participant numbers.
No clinical covariates are included.

## capscan_assembly_dna_yield_gain_by_acquisition.tsv
sha256: 351bdf31dc30f108e11796f15614cb05299ba8d37a8991d6ca4f10f7fa1a9c4d
rows: 282 (one row per matched CapScan acquisition; 226 intestinal fluid, 56 stool; 250 specimens; 15 participants)
produced by: capsule metagenome join (see `figures/SFig09/UPSTREAM.json`; join code not saved upstream)

| column | unit | meaning |
|---|---|---|
| sample | id | acquisition (raw file stem) |
| region | id | sampling region: duodenum, jejunum, ileum, colon (Regions 1-4 in the figure) or stool |
| subject | id | CapScan participant number |
| biosample | id | NCBI BioSample of the specimen whose assembly forms the MG database |
| assembly_proteins | count of predicted-protein records | MG database size, sanitizer `n_kept`, human reference excluded |
| dna_ng | ng | recorded DNA quantity from the specimen metadata; empty for 183 acquisitions |
| volume_ul | µl | recorded specimen volume from the specimen metadata; empty where not recorded; not used by the figure |
| FL_peptides | count of peptides | FastaLake peptides, 5-50 residues, q <= 0.01 |
| MG_peptides | count of peptides | MG-search peptides, 5-50 residues, q <= 0.01 |
| gain_percent | % | 100 x (FL_peptides / MG_peptides - 1), rounded to 3 decimals |
| MG_identity_recovery | fraction | share of MG-search peptides also identified by FastaLake, rounded to 4 decimals; not used by the figure |

## capscan_acquisitions_9_50.tsv (v30 figure input)
sha256: 5bb9ec898849e853e93ef0d01da9b7eaebab64db2b1a075b89f8d9d44fe5ad39
rows: 282 (as `capscan_assembly_dna_yield_gain_by_acquisition.tsv`)
produced by: `figures/SFig09/panels.py:acquisitions_9_50` (lab: len950 build.py) from that table and `a08_benchmark282_old.tsv`

Columns and meanings as in `capscan_assembly_dna_yield_gain_by_acquisition.tsv`, except:

| column | unit | meaning |
|---|---|---|
| FL_peptides | count of peptides | FastaLake peptides, 9-50 residues (`n_fl_9` of `a08_benchmark282_old.tsv`) |
| MG_peptides | count of peptides | MG-search peptides, 9-50 residues (`n_mg_9`) |
| gain_percent | % | 100 x (FL_peptides / MG_peptides - 1) at 9-50, rounded to 3 decimals |
| MG_identity_recovery | fraction | unchanged 5-50 value; not used by the figure |

This is the table the v30 figure is drawn from. When the CapScan MG benchmark is extended to
295 acquisitions, a replacement with the same columns is passed by `SFIG09_ACQ_TSV`.

## capsule_assembly_gain_tertiles_9_50.tsv (expected panel b values, v30)
sha256: 761b3e620045ad96c4cdcaefb49ca4eb3bc6b8eeefb48c3e4e4b2719b8f2638f
rows: 3

| column | unit | meaning |
|---|---|---|
| assembly_tertile | id | Smallest, Middle or Largest third of intestinal-fluid MG database size |
| median_gain | % | median acquisition-level gain at 9-50 residues |
| n | count of acquisitions | acquisitions in the third |
| min_assembly, max_assembly | count of predicted-protein records | bounds of the third |

## capsule_assembly_gain_per_specimen.tsv
sha256: 800c0e370c61f8557539c1ec149d981965c49c893797910301108ca49cc144e9
rows: 250 (one row per specimen)
produced by: aggregation in `figures/SFig09/panels.py:derive_tables` (ported from the v15 figure builder)

| column | unit | meaning |
|---|---|---|
| biosample | id | specimen |
| subject | id | CapScan participant number |
| region | id | sampling region or stool |
| assembly_proteins | count of predicted-protein records | median MG database size over the specimen's acquisitions |
| gain_percent | % | median gain over the specimen's acquisitions (unrounded) |
| acquisitions | count of acquisitions | acquisitions of this specimen |

## capsule_assembly_gain_acquisitions.tsv
sha256: 89956d96f81cc50656da184c078c6f14b50e44381b4c1fe3131eb7eca4ddc1c9
rows: 226 (one row per intestinal-fluid acquisition)
produced by: aggregation in `figures/SFig09/panels.py:derive_tables`

Columns as in `capscan_assembly_dna_yield_gain_by_acquisition.tsv`, plus

| column | unit | meaning |
|---|---|---|
| gain_exact | % | 100 x (FL_peptides / MG_peptides - 1), unrounded |
| assembly_tertile | id | Smallest, Middle or Largest third of assembly_proteins (pandas.qcut, 3 equal-count bins) |

## capsule_assembly_gain_tertiles.tsv
sha256: 9d9b3e2bf774245e74406e5111e69ad02654eebe4eea746641c48b4c71144796
rows: 3 (one row per third of MG database size)
produced by: aggregation in `figures/SFig09/panels.py:derive_tables`

| column | unit | meaning |
|---|---|---|
| assembly_tertile | id | Smallest, Middle, Largest |
| median_gain | % | median gain_exact over the acquisitions in the bin |
| n | count of acquisitions | acquisitions in the bin |
| min_assembly | count of predicted-protein records | smallest MG database size in the bin |
| max_assembly | count of predicted-protein records | largest MG database size in the bin |

## CAPSCAN_sanitize_report.tsv
sha256: 0120dd064a2b20ab39c59177519cb637aeb701bbcd392ca691260b630a9f7205
rows: 335 (one row per CapScan specimen assembly)
produced by: assembly sanitizer of the FastaLake chain (see `figures/SFig09/UPSTREAM.json`)

| column | unit | meaning |
|---|---|---|
| sample | id | NCBI BioSample of the assembly |
| n_in | count of records | predicted-protein records in the input FASTA |
| n_kept | count of records | records kept (the MG database size used here) |
| n_malformed | count of records | malformed records removed |
| n_dup_identical | count of records | duplicate accessions with identical sequence removed |
| n_dup_conflicting | count of records | duplicate accessions with differing sequence removed (first kept) |

## a08_benchmark282_old.tsv
sha256: 4c732679b085a01245d39f324d0884d2eff1986f73914bdd9f37c7ddec2eee44
rows: 282 (one row per matched CapScan acquisition)
produced by: capscan_public_only_search (recount from the per-acquisition peptide sets of the production arm)

| column | unit | meaning |
|---|---|---|
| sample | id | acquisition |
| group | id | sampling region or stool |
| n_fl_5, n_mg_5 | count of peptides | FastaLake and MG peptides, 5-50 residues, q <= 0.01 |
| rec_5 | fraction | MG peptides (5-50) also identified by FastaLake |
| n_fl_9, n_mg_9 | count of peptides | as above, 9-50 residues |
| rec_9 | fraction | as rec_5, 9-50 residues |

## a08_benchmark282_new.tsv
sha256: fa6eaf0fddd45be5a5dc989c5c64397421a82713beb0fd084a7dc9aea1700265
rows: 282 (one row per matched CapScan acquisition)
produced by: capscan_public_only_search, rerun on the 347-acquisition roster (v28f)

Columns as in `a08_benchmark282_old.tsv`.

## capscan_roster_347.tsv
sha256: 473ff6afcf7348163f1fbaa8d624a8078eed9348a15019952363bc86e70f7aa5
rows: 347 (one row per CapScan acquisition in the v28f roster; 15 participants)
produced by: v28f CapScan roster (participant-15 addition); used only to cross-check region and participant

| column | unit | meaning |
|---|---|---|
| sample_id | id | acquisition (raw file stem) |
| region | id | sampling region or stool |
| region_ordinal | dimensionless | 1-4 for duodenum to colon, 5 for stool as given by the roster |
| subject_id | id | CapScan participant number |
