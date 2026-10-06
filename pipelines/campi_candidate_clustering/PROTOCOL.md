# CAMPI clustering extension: 90%, 95%, 75% and 50% identity

This extends the completed H04/H08 97%/99% clustering sensitivity to 90%, 95%, 75% and 50% identity, with functional coverage. No other output is modified.

All 15 CAMPI acquisitions use the exact pinned broad-UniProt candidate FASTAs, original predictions, razor/hash-acc selector, MMseqs/Linclust binary and Sage0.14.6 search configuration from the prior comparison. The intervention is clustering candidate proteins AFTER broad-reservoir lookup and BEFORE parsimony. It does not test clustering the whole reservoir. Coverage0.8, cov-mode0, cluster-mode2, sorted accession input, createdb shuffle0 and two threads stay fixed.

Each acquisition must reproduce the unclustered selected sequence set. S01 and the largest candidate input S06 are serial pilots; both also reproduce the existing99% representative and selected sequence sets. Only90%,95%,75% and50% receive new Sage searches. Exact source checks precede computation. Failure blocks the dependent main array; retain all failures.

Primary denominator is the original fixed contemporaneous H08 unclustered search, shared with97%/99%. S02 uses the saved baseline from the documented recovery. Its one-gain/one-loss confidence-fit difference from the archive is retained; no new baseline is chosen. Export archive/S02-repeat sensitivities. Recount frozen old accepted sets, preserve every zero and gain/loss, compare peptide identities, exact selected-sequence availability and member-specific REF evidence. Group counts from native Sage outputs are not substituted for peptide endpoints.

Report acquisition-paired peptide count ratios and100*(1-ratio), database entries, REF identity recovery and gains/losses. Combine90/95/97/99 against the same baseline. Confidence remains native peptide q<=0.01 and canonical5–50 residues, following the source comparison. Check eligibility against the stored baseline sets. Observations are conditional search yields, not equal-error sensitivity, strain truth or whole-workflow FDR calibration.

Execution: low priority (nice10000), 2CPUs/32G, six-hour cap per task, one task at a time. Pilot (two acquisitions sequentially) -> full remaining13 array capped1 -> collector. Submission is not completion. Full receipt, source checks and independent summary reconstruction are required before numbers or figures are used. Existing97/99 results and adverse findings remain preserved.

## Functional-coverage extension
75% and 50% serve as a community benchmark. Use one frozen existing exact-sequence annotation dictionary from CAMPI FL+REF+MG. Report compatible KEGG Orthology (KO) coverage and the stricter single-KO peptide evidence separately. Keep unmapped sequences, mapped peptides without KO and multi-KO ambiguity explicit. Do not transfer a cluster representative annotation to its members. This estimates retained annotated function evidence within the fixed dictionary, not all true functions, activity, abundance or functional equivalence of homologs. Compare unclustered,99%,97%,95%,90%,75%,50%; export per-acquisition retained/lost/gained KO identities and paired summaries.
