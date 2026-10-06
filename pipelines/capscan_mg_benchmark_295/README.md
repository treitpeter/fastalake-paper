# capscan_mg_benchmark_295

Internal ID: mg295 (working-directory name; no lab code assigned)

Status: analysis running; code and receipts to be added.

Question. Does the public-sequence FastaLake versus matched-metagenome comparison in CapScan
hold when the 13 participant-15 acquisitions whose specimen has a completed public shotgun
run are added, giving 295 matched acquisitions on the 347-acquisition roster instead of 282
on the 331-acquisition roster? See the DNA count reconciliation in
`pipelines/capscan_participant_roster/README.md`.

Planned contents (from the running analysis):

- matched MG FASTA for the 13 specimens, built with the production rehash and human
  Swiss-Prot concatenation and gated by sha256 against a rebuilt existing pair (control) and
  a one-residue change (negative control);
- same-spectrum Sage MG searches for the 13;
- spectrum-split (disjoint) and 10-fold cross-fit comparisons with the production worker
  code and an extended 591-row roster;
- MSFragger and Percolator second-engine searches for the 13;
- comparison against the 282-acquisition results;
- Slurm launchers parametrised by environment variables, README Choices table, RECEIPTS.md
  with job ids and output sha256.

Producer location: `LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/final/mg295/`.
