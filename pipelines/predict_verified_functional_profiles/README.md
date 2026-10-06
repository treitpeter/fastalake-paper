# predict_verified_functional_profiles

Internal ID: A88 (previous folder name A88_predict_verified_functional_profiles)

Re-derive the PREDICT peptide function evidence, quantity matrices and OG, KO and pathway profiles when a native MGnify UHGG eggNOG row is kept only if its own seed orthologue aligns to the attached UHGG sequence (DIAMOND bitscore at least 50% of the bitscore MGnify reported, variant `verified50`). Used by STab5 (306 acquisitions, non-host profiles `*_other_unresolved`) and by A70, A84, A87, A108.

Status: code copied as it ran (path edits only), not re-run in this repository. See RECEIPTS.md for original job ids and hashes.

## Inputs

Frozen evidence tasks of the grouping step (`evidence_all_cohorts/INPUTS.json`, task `PREDICT_noclust`: feature ledger, peptide table, acquisition summary and annotation inputs, each with sha256) and the annotated-origin consensus policy. Task entries hold absolute lab paths; the manifest hash is 458fba5f52166d97336f5cb350e94a071e2e4f2a6f9db830b54e518e23287b52.

## Outputs

`profiles/verified50/PREDICT_noclust/`: `KO_all.npz`, `KO_other_unresolved.npz`, `OG_*`, `Pathway_*` with axes files, `COMPLETE.json` (replay error, change counts). Hashes in RECEIPTS.md.

## Slurm resources and run time (original logs, `sacct`)

Job 6048301: 4 cores, 40G requested, peak RSS 16.09 GB, wall 29 min 11 s. Command `corrected_reproject.py PREDICT_noclust verified50`. The matching replay variant (`replay`, keep every native row) reproduces the frozen outputs exactly before the verified variant is used.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| verification threshold | DIAMOND bitscore at least 50/100 of the MGnify bitscore (`verified50`) | removes circular native orthologue rows | `verifiedNN` for other NN; only 50 is used by the paper |
| rows kept | native rows pass: 8,085,318; fail: 286,165 (of 8,371,483 native accessions) | from the job log | none |
| allocation | equal share to distinct eligible annotations; human-assigned signal excluded; unresolved origin retained | protocol of the parent run | none |
| frozen | taxonomy, parent graph, peptides, acquisition roster, LFQ ledger re-checked | no regrouping | none |

## Notes

Code edits relative to the original: lab absolute paths replaced by `FASTALAKE_DATA_DIR`; lab-tree roots located by a marker file replaced by `FASTALAKE_DATA_DIR`. No statistics changed. The scripts keep their original integrity gates (`verify()`, `gate()`, hash checks against `*_CODE.json` and `*_INPUTS.json` manifests of the frozen run). Those manifests contain absolute lab paths and are not in the repository; the hashes they held are listed in RECEIPTS.md. Because the code here carries the path edits described below, its own sha256 differs from the original (both are in RECEIPTS.md), so a re-run needs the manifests regenerated for this copy. The projection modules `run_projection.py`, `evidence_rules.py`, `policy.py` and `profiles.py` are imported by `corrected_reproject.py` from lab folders and are copied here because they are on the call chain. `run_projection.py` needs the `pyahocorasick` package. The original `run.sh` is only `exec python -u "$@"` and is replaced by the launcher.
