# Corrected Sage build for the joint healthy-stool searches

The joint searches of the 58 healthy-stool acquisitions used Sage built from the upstream v0.14.7 source with one correction: pooled LFQ q values are assigned back to peaks keyed by (precursor, decoy), so target and decoy values no longer overwrite each other. All other searches in the paper used the unmodified Sage v0.14.7 release.

## Identity

| item | value |
|---|---|
| upstream | https://github.com/lazear/sage, tag `v0.14.7`, commit `99407db6e3754b31a9b88b7316a0aee67293c93f` (crates report version 0.14.6) |
| patch | `sage-v0.14.7-fastalake-lfq-identity-producing.patch` (3 files: `crates/sage/src/fdr.rs`, `crates/sage-cli/src/main.rs`, `crates/sage-cli/src/output.rs`) |
| full producing source | `sage-v0.14.7-fastalake-lfq-identity-producing-source.tar.gz` (upstream tree with the patch applied, upstream MIT `LICENSE` inside; also `LICENSE.sage`) |
| lockfile | `Cargo.lock` (original upstream lockfile, used with `--locked`) |
| source manifest | `PRODUCING_SOURCE_MANIFEST.json` (sha256 of the 41 source files); `verify_source.py <tree>` checks a tree against it |
| toolchain | Rust/Cargo 1.75.0 (the pinned `time` dependency does not build with 1.91.1) |
| build job | Slurm 5922884 (build and regression test `mm_lfq_identity_regression`, passed); first attempt 5922849 |
| searches | Slurm 5923007 (pilot), 5923008 (joint search) |

## Relation to upstream commit 0e2c7a6

Upstream fixed the same defect after v0.14.7 in lazear/sage commit `0e2c7a6` ("fix: bug in MS1 picked-precursor FDR (#150)"), which changes only `crates/sage/src/fdr.rs`: the q-value map is keyed by `(score.ix, score.decoy)` and peaks are looked up by their full `(PrecursorId, decoy)` key. The patch here makes the same change to the same function (`picked_precursor`), written slightly differently (the closure binds the full key as `key` instead of `(ix)`), so the q-value assignment is identical. The patch is not byte-identical to 0e2c7a6: it also adds (1) a regression test with separate target and decoy q-value expectations, (2) an optional ranked-score audit written when `MM_LFQ_RANK_AUDIT` is set, and (3) an `lfq_target_decoy_audit.tsv` export. None of these changes scoring, peak integration or FDR estimation. Checked 2026-10-06 against a clone of the upstream repository: `git apply` of the patch on 99407db succeeds and the result matches all 41 manifest hashes.

## Build

    git clone https://github.com/lazear/sage && cd sage && git checkout 99407db6e3754b31a9b88b7316a0aee67293c93f
    git apply --check <this folder>/sage-v0.14.7-fastalake-lfq-identity-producing.patch
    git apply <this folder>/sage-v0.14.7-fastalake-lfq-identity-producing.patch
    cp <this folder>/Cargo.lock .
    python3 <this folder>/verify_source.py .        # expect PASS, 41 files
    cargo +1.75.0 test --locked -p sage-core mm_lfq_identity_regression -- --nocapture
    cargo +1.75.0 build --locked --release -p sage-cli
    target/release/sage --version                    # sage 0.14.6

Or unpack the source archive and run the last three commands in it. A rebuild can give a binary with different bytes (paths, compiler, linker); the hash below identifies the executable that ran.

## Checksums

| file | sha256 |
|---|---|
| producing executable `target/release/sage` (LAB_ROOT/MM_review_analysis/A12_7_quantification/source/target/release/sage, 21,591,120 bytes; not redistributed) | `21241c22a7e886789e3abfb3f129c299ec06080b8ccb504d54906e17689a3420` |
| `sage-v0.14.7-fastalake-lfq-identity-producing.patch` | `303c49ec41a82c4b8ca65e945a3faf103059509d65ade92c77b70dd7a0022558` |
| `Cargo.lock` | `57528de307f3574b8fd0615b5455fc1afe6d1cc25c5fbac05f1a56740eb36f70` |
| `sage-v0.14.7-fastalake-lfq-identity-producing-source.tar.gz` | `4b9e555a8610e416a1867c24f271f7ff226fdf1ed77a73c7e1c7bbbdb17fef0d` |
| `PRODUCING_SOURCE_MANIFEST.json` (this copy; the `source` field was changed to a LAB_ROOT path, original sha256 `a07182c711d63433253f64521d25d2f2ad15824012485088b51120c57ae96934`) | `108f01c766cd55af26ef9244969ac493105a7eaef01aef5846b7573ac239d818` |

`SHA256SUMS.txt` lists every file in this folder (`sha256sum -c SHA256SUMS.txt`).

Origin of these files: LAB_ROOT/release_prep/fasta-lake-paper_2026-10-06/current/v15/evidence/sage_lfq_fixed/ (review bundle) and LAB_ROOT/MM_review_analysis/A12_7_quantification/ (build and search jobs).
