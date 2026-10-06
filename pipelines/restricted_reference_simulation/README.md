# restricted_reference_simulation

Internal ID: A71 (previous folder name A71_restricted_simulation)

## Question

This pipeline produces: known-truth simulation outcomes (simulation_all_outcomes.tsv, simulation_mixture_outcomes.tsv), bootstrap summaries (simulation_uncertainty.tsv) and design constants, line 733

Used by: figures/SFig05, text_statistics/supplementary_methods_6_7.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `generate.py` | Construct matched mixtures and independent spectral perturbations; no search. |
| `evaluate.py` | Truth evaluation after all searches, with a separately implemented recount. |
| `collect.py` | Collect all planned units, preserving mixture-level replication and uncertainty. |
| `common.py` | Shared paths, provenance and small format helpers for this experiment. |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `simulation_all_outcomes.tsv` | not found on disk |
| `simulation_mixture_outcomes.tsv` | not found on disk |
| `simulation_uncertainty.tsv` | not found on disk |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `PROJECT` | `ROOT.parents[1]` | no comment in original code | not stated in code; see common.py |
| `PREVIOUS` | `PROJECT / "analyses/fdr_validation_2026-09-08"` | no comment in original code | not stated in code; see common.py |
| `FUSION` | `PROJECT / "analyses/fusion_entrapment_2026-09-07"` | no comment in original code | not stated in code; see common.py |
| `RELEASE` | `PROJECT / "release_prep/github_biorxiv_2026-09-07"` | no comment in original code | not stated in code; see common.py |
| `RUNTIME` | `RELEASE / "logic_review_2026-09-10/campi_clean/_runtime/bin"` | no comment in original code | not stated in code; see common.py |
| `REFERENCE` | `PREVIOUS / "results/campi_reference_v2/REF_targets.fasta"` | no comment in original code | not stated in code; see common.py |
| `NOISE` | `PREVIOUS / (` | no comment in original code | not stated in code; see common.py |
| `CONFIG` | `PROJECT / (` | no comment in original code | not stated in code; see common.py |
| `THRESHOLDS` | `(0.001, 0.005, 0.01, 0.02, 0.05)` | no comment in original code | not stated in code; see common.py |
| `N_SPECTRA` | `2048` | no comment in original code | not stated in code; see common.py |
| `PART_SIZE` | `128` | no comment in original code | not stated in code; see common.py |
| `ACQUISITIONS` | `("donor", "test")` | no comment in original code | not stated in code; see common.py |
