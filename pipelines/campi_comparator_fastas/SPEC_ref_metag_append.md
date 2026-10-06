# CAMPI final-FASTA ablation: REF_MG

All 15 acquisitions from the canonical CAMPI mock manifest: S01–S10, S11_Fraction1–4, and S12.

Sequence flow: broad UniProtKB-expanded canonical noclust lake → evidence selection → sample selection → razor/hash-acc → complete CAMPI reference addition → complete CAMPI metagenome addition → unchanged Sage 0.14.6.

This arm adds the full frozen sequence set in `inputs/REFERENCE.json` to each acquisition's final FASTA from `LAB_ROOT/canonical/publication_data/arms/uniprotkb_2026_03_noclust_campi_ref_append_2026-09-09`. The REF+MG layer therefore also retains the entire REF layer. Exact normalized sequences are deduplicated using the committed Rust lake_builder; original base accessions have priority. Distinct I/L sequences are retained. Every added source accession maps to its final search representative in `4_reference_append/<sample>/reference_to_search_accession.tsv`. Namespace is CAMPI_MG__. No selection, razor, clustering or species filtering runs after addition. Sage generates fresh decoys under its unchanged settings.

The explicit `4_reference_append` output is symlinked into the canonical `4_parsimony/razor_hash-acc` contract solely for unchanged search, annotation and dual-key consumers. This path does not imply that the augmented database underwent another parsimony step. Cohort outputs live separately from all existing arms. Upstream base lake is reused, with min de novo length 9 and top fraction 0.30.

The complete REF comparator target file includes its contaminant targets; the MG file includes its target contaminants too. Full target sets are intentional; existing decoys are excluded. Broad UniProt contains potential foreign-entrapment species, so old entrapment interpretations are not reused.

Preflight 5895451 used archived S01 only to exercise the full append. These preflight counts are not experimental results for the new UniProt run. Production inputs depend on current base-arm completion. Production code commit: 8bc8dae965fc74dea58c346780b66576f5309875. Submission commands and dependencies: `jobs.json`. The final validation job checks all 15 searched FASTAs, nonempty results/LFQ and paired Sage settings/raw paths against the blind arm. Scientific comparison awaits successful outputs and canonical aggregation/grouping review.

For the cumulative ablation compare blind vs REF vs REF+MG with paired acquisitions. Keep database size and identified accessions distinct from harmonized groups. Return Supplementary Figure 8 a–d to main Figure 2 using recomputed data; omit e, as recorded in the Figure 2 panel decision.
