# native_crossomics_profiles

Internal ID: A86 (previous folder name A86_native_crossomics_profiles)

How well do DNA and RNA gene-TPM profiles, host and non-host protein intensity and metabolites agree in sample geometry without summing to OG, KO or pathway? Genes are keyed by the SHA-256 of their full encoded amino-acid sequence so that sample-local assembly names are not joined by name. Used by SFig11 and the native-profile statements of supplementary methods 8-11 (retained signal fraction of the 10% prevalence filter, native feature counts, CapScan DNA run membership).

Status: code copied as it ran (path edits only), not re-run in this repository. See RECEIPTS.md for original job ids and hashes.

## Inputs

Frozen input manifest `native/INPUT_MANIFEST.json` (sha256 `437c76019adcf128d7f2bdf5d7caaf97222b76c15bc1d0c1df8bc925d2df3fb8`) lists every FASTA, TPM table, Salmon checkpoint and protein source with its hash; it holds absolute lab paths and is not in the repository.

## Outputs

`native/results/{PREDICT,CAPSCAN}/` (`correlations.tsv`, `feature_quality.json`, `DNA_native_cells.npz`, `RNA_native_cells.npz`, `gene_sequence_features.tsv.gz`, rosters) and `SOURCE_DATA_retained_signal_by_profile.tsv` (hashes in RECEIPTS.md).

## Slurm resources and run time (original logs, `sacct`)

| cohort | job | cores | requested mem | peak RSS | wall |
|---|---|---|---|---|---|
| CapScan | 5990406 | 1 | 32G | 6.19 GB | 20 min 16 s |
| PREDICT | 5990407 | 1 | 32G | 4.04 GB | 13 min 7 s |

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| gene coordinate | SHA-256 of the full amino-acid sequence, uppercase, stop removed; no I/L folding, no clustering | sample-local gene names cannot be joined | none |
| duplicate sequences within a specimen | summed | identical proteins | none |
| unresolved PREDICT TPM genes | exported and counted separately; resolved plus unresolved TPM conserves input | 2,327 genes of one sample are absent from its archived FASTA; not guessed | strict gate kept for CapScan |
| prevalence filter | at least 10% of units (primary); all observed features (sensitivity) | existing framework | both reported |
| transform | relative abundance, square root; metabolites positive log10, at least 70% finite, median imputation, constants removed, standardized | existing framework | none |
| distance | Euclidean between samples; Spearman between upper triangles | existing framework | none |
| permutations | 1999, seed 20260907, two-sided plus-one p | protocol | none |
| CapScan region | per-feature region mean removed; whole participant blocks permuted | repeated measures | none |
| BH | within each new primary pyramid: 10 tests PREDICT, 6 CapScan | complete families | sensitivity families separate |

## Notes

Code edits relative to the original: lab absolute paths replaced by `FASTALAKE_DATA_DIR`; lab-tree roots located by a marker file replaced by `FASTALAKE_DATA_DIR`. No statistics changed. The scripts keep their original integrity gates (`verify()`, `gate()`, hash checks against `*_CODE.json` and `*_INPUTS.json` manifests of the frozen run). Those manifests contain absolute lab paths and are not in the repository; the hashes they held are listed in RECEIPTS.md. Because the code here carries the path edits described below, its own sha256 differs from the original (both are in RECEIPTS.md), so a re-run needs the manifests regenerated for this copy.Commands in the original run: `prepare`, `pilot`, `run --cohort PREDICT`, `run --cohort CAPSCAN`. `native.py` is the expanded-cohort runner (16-participant extension, `prepare|pilot|full`); its outputs are not mapped to a specific published number here (see OPEN_ITEMS.md). `crossomics.py` of A84 and `distance_sensitivity.py` of A85 import `profile_distance`, `bh` and `correlations` from `native_pyramids.py` by parsing.
