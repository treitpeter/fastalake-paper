# Code availability

FastaLake is released under the MIT licence at
[MannLabs/fasta-lake](https://github.com/MannLabs/fasta-lake). The manuscript cites
**v1.1.0rc2, commit 9ad09c9**, archived at
[doi:10.5281/zenodo.23215189](https://doi.org/10.5281/zenodo.23215189).

The study code is available under the MIT licence at
[treitpeter/fastalake-paper](https://github.com/treitpeter/fastalake-paper).
The repository includes the [improvements listed here](../CHANGELOG.md), which
will be included in the next bioRxiv revision. The current manuscript reference
is recorded in [MANUSCRIPT.json](MANUSCRIPT.json). Cite the analysis commit used;
the study-code repository does not yet have its own archival DOI.

## Reproduction

`pipelines/` contains upstream analysis code; `figures/`, `tables/` and
`text_statistics/` contain the notebooks and their input tables. The
[reviewer guide](REVIEWER_GUIDE.md) maps registered claims to code and inputs.
An entry in that guide is not evidence of complete manuscript coverage: unresolved
entries remain explicit, and a public run reports computations requiring unavailable
large or controlled inputs as skipped. See the [check documentation](../README.md).

The Docker image `ghcr.io/mannlabs/fasta-lake:v1.1.0rc2` includes the Rust programs
and official Sage v0.14.7 executable. All searches except the joint healthy-stool
searches used that Sage release (commit 99407db; the executable reports 0.14.6).
Joint healthy-stool searches used a build of the same source with a correction
separating target and decoy q values when assigning pooled LFQ peak q values.
The source, patch, lockfile and checksums are in
`pipelines/healthy_stool_joint_search_q/sage_build/`. The corresponding upstream
fix is [lazear/sage commit 0e2c7a6](https://github.com/lazear/sage/commit/0e2c7a6).
The local patch also includes validation code and is not byte-identical to that
upstream commit.

For the historical Fig. 4h display, the selected-feature list and display code are
preserved. The original feature-selection code and complete selection-screen output
were not recovered; the later screen is in `pipelines/predict_eds_screen/`.
