# CAMPI fusion control on the reported broad reservoir — 22 September 2026

The historical-reservoir control was judged insufficient support for the current benchmark, so the matched experiment was run on Slurm. This work is restricted to CAMPI.

## What changes

The historical experiment used 233,477,360 sequences. This rerun uses the exact 375,333,654-entry stage-0 FASTA from `uniprotkb_2026_03_noclust_2026-09-09`, the same 15 prediction CSVs and mzML files, that arm's acquisition-specific search settings, and production extractor/parsimony/Sage binary hashes. No CAMPI REF or MG is added to the reservoir. Two prespecified seeds, 20260907 and 20260908, give 30 searches. No favourable-value gate or adaptive seed choice is used.

The original control construction and estimator are retained so the reservoir mismatch is addressed directly. Every reservoir protein is linked to its seeded composition-preserving, cleavage-compatible control before pooled evidence lookup. Both partners can supply selection evidence and pass together through acquisition lookup and razor/hash-accession parsimony. They become separate target/control records only after selection. Controls are not appended to the archived final FASTAs. Protein pairs, ambiguous sequence overlaps and non-bijective peptide mappings remain visible. The historical pairing audit includes its declared initiator-methionine-clipping sensitivity; final Sage search chemistry itself is copied exactly from production. See `HISTORICAL_SPEC.md` for the original construction and its statistical assumptions.

## Provenance and gates

`inputs/PROVENANCE.json` records paths, sizes and nanosecond mtimes for all frozen inputs, hashes of original search configurations and selected FASTAs, and the verified production binary hashes. Configurations and binaries are copied locally. `SOURCE_SHA256.json` freezes the runner, independent auditor, driver, tests, configurations and protocol.

1. Preflight runs the existing transformation, collision/accounting and control-only-evidence-retention fixtures, reproduces S01's unmodified production selected-protein sequence multiset, then computes full reservoir/prediction/mzML SHA-256 hashes on the compute node.
2. Two preparation tasks transform and search the complete broad reservoir, one per seed. Each confirms all 375,333,654 linked entries were searched.
3. Two S01 pilot searches, one per seed, use the frozen production search configuration. Independent final-FASTA pairing reconstruction and raw-result recount must agree exactly with the original estimator. Gates check construction and accounting, not whether the diagnostic is below 1%.
4. The remaining 28 searches are dependency-gated on both pilots, with at most six concurrent tasks. Each receives the same independent audit.
5. Collection verifies 30 validated searches and all per-search output hashes. It reports acquisition results, pooled counts within seed, paired/unpaired components and acquisition medians. Length populations are 5–50, nested 9–50 and 5–8; q thresholds are .001/.005/.01/.02/.05. The shorter-length sensitivity does not recalculate q values.

## Scheduling and endpoints

Low priority (`--nice=10000`) on p.hpcl92/p.hpcl93/p.hpcl8. Preflight: 8 CPUs, 96 GB, 2 h. Preparation: two tasks maximum, 16 CPUs/96 GB each, 6 h. Searches: 8 CPUs/96 GB, 1 h each; two pilots, then six concurrent remaining tasks. Collector: 2 CPUs/8 GB, 30 min. Phases use `afterok`; failed pilots cannot release the full array. Existing jobs are preserved.

The actual queue was checked before submission; no CAMPI fusion rerun was active. Preparation scans the full 136-GB reservoir and writes complete paired reservoirs, requiring roughly 0.6 TB for both seeds plus selected databases and search outputs. The checked filesystem had 196 TB free. No original output is overwritten or removed.

Completion means the matched experiment and independent numerical audits pass. Scientific interpretation remains a separate endpoint. This corrects the reservoir/configuration mismatch but does not make an assumption-dependent fusion diagnostic a certified error rate for ordinary target-only searches or protein groups. All outcomes and the historical results are preserved. No manuscript result is replaced until the new results are reviewed.
