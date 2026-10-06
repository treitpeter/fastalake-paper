# clinical_group_ko_homogeneity

Internal ID: homogeneity_v30 (working-directory name; no lab code assigned)

Status: analysis running; code and receipts to be added.

Question. Are the KEGG orthology (KO) annotations of the proteins that can explain an
accepted peptide more homogeneous with public-sequence FastaLake than with the matched
metagenome database, in the 296 PREDICT and 282 CapScan matched acquisitions?

Definitions (identical for both arms): accepted peptide = Sage rank 1 target, peptide
q <= 0.01, 9-50 residues after removing modifications, I/L-equivalent; peptide group = the
distinct I/L-normalised sequences of that search's FASTA that contain the peptide; KO =
eggNOG-mapper 2.1.12 (eggNOG 5.0.2, DIAMOND mode, production settings) looked up by sha256
of the I/L-normalised sequence, one annotation source for both arms.

Planned contents (from the running analysis):

- `hom.py` with subcommands roster, extract (per acquisition and arm), index (KO index and
  the sequences still lacking annotation) and metrics (per-acquisition and summary tables);
- Slurm launchers for extraction (578-task array), indexing, eggNOG-mapper on the
  unannotated sequences (12 shards) and metrics, parametrised by environment variables;
- README Choices table and RECEIPTS.md with job ids and output sha256.

Producer location: `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/homogeneity_v30/`.
