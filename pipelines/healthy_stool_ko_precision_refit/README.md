# healthy_stool_ko_precision_refit

Internal ID: A90 (previous folder name A90_hstool_ko_precision_refit)

## Question

This pipeline produces: verified50/hstool_precision_refit/pairwise_native_features.tsv (ko_pairwise.tsv): Fig3 f KO row; measure()/pairs() used by A107 for healthy-stool KO correlations

Used by: figures/Fig3, text_statistics/supplementary_methods_1_5.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `hstool_precision.py` | Reproduce the frozen KO precision estimator, then compare revised annotations. |
| `run_hstool_precision.sh` | launcher |
| `profiles.py` | Conserved, distinct-term reporting profiles on the frozen expanded reference. |
| `run_downstream.sh` | Run the frozen biology_refit drivers (byte-identical copies) in the shadow for one variant. |
| `make_shadow.sh` | Build a read-only shadow of biology_refit (B) and the hstool driver (C) whose profile inputs point at |
| `prep_native.py` | Prepare the non-circular native UHGG eggNOG verification. |
| `verify_shard.py` | Non-circular check of each native MGnify eggNOG row: does the protein sequence the row is attached to |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `pairwise_native_features.tsv` | `2cf176c091f92f89162acfb091c3e3da6ca04886380f564c158040ac521a01f3` |
| `pairwise_native_features.tsv` | `0d53beafd1e1d0625b8de8269f33474ef2d414f440510cdc9a69ff4e31e2e828` |
| `pairwise_native_features.tsv` | `2cf176c091f92f89162acfb091c3e3da6ca04886380f564c158040ac521a01f3` |
| `pairwise_native_features.tsv` | `0d53beafd1e1d0625b8de8269f33474ef2d414f440510cdc9a69ff4e31e2e828` |
| `pairwise_native_features.tsv` | `0d53beafd1e1d0625b8de8269f33474ef2d414f440510cdc9a69ff4e31e2e828` |
| `pairwise_native_features.tsv` | `2cf176c091f92f89162acfb091c3e3da6ca04886380f564c158040ac521a01f3` |
| `pairwise_native_features.tsv` | `0d53beafd1e1d0625b8de8269f33474ef2d414f440510cdc9a69ff4e31e2e828` |
| `pairwise_native_features.tsv` | `52f9a14957c72db844a3689930e9571a611c6c6f8c39e71ade2e77737e76d3c8` |
| `pairwise_native_features.tsv` | `2cf176c091f92f89162acfb091c3e3da6ca04886380f564c158040ac521a01f3` |
| `pairwise_native_features.tsv` | `2cf176c091f92f89162acfb091c3e3da6ca04886380f564c158040ac521a01f3` |
| `pairwise_native_features.tsv` | `2cf176c091f92f89162acfb091c3e3da6ca04886380f564c158040ac521a01f3` |
| `pairwise_native_features.tsv` | `2cf176c091f92f89162acfb091c3e3da6ca04886380f564c158040ac521a01f3` |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `A33` | `R/'MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19'` | no comment in original code | not stated in code; see hstool_precision.py |
| `OLD` | `A33/'results/HSTOOL__simplified/source_data/Figure3'` | no comment in original code | not stated in code; see hstool_precision.py |
| `NEW` | `H/'profiles/HSTOOL_noclust'` | no comment in original code | not stated in code; see hstool_precision.py |
| `UNASSIGNED` | `'Unassigned_at_this_resolution'` | no comment in original code | not stated in code; see profiles.py |
| `HUMAN` | `'Human_assigned_signal_withheld'` | no comment in original code | not stated in code; see profiles.py |
| `LEVELS` | `['KO','OG','Pathway']` | no comment in original code | not stated in code; see profiles.py |
| `INP` | `json.loads((E/'INPUTS.json').read_text())` | no comment in original code | not stated in code; see prep_native.py |
