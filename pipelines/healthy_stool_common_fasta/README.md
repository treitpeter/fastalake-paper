# healthy_stool_common_fasta

Internal ID: A21 (previous folder name A21_hstool_quant_comparison)

## Question

This pipeline produces: results/identity_comparison.tsv (SFig07 a study-level, repeat and common-FASTA columns); results/leaders/technical_agreement_summary.tsv (line 769); common FASTA union summary 153,267 (line 751)

Used by: figures/SFig07, text_statistics/supplementary_methods_6_7.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `search_revision2.py` | Reconstruct identified, quantified and transferred HSTOOL evidence separately. |
| `compare_revision2.py` | Matched HSTOOL identity, feature/group precision and missingness comparisons. |
| `leaders_revision2.py` | Reproduce baseline leaders, then audit fixed-dictionary and native-set alternatives. |
| `slurm/compare.sbatch` | launcher |
| `slurm/compare_revision2.sbatch` | launcher |
| `slurm/extract.sbatch` | launcher |
| `slurm/extract_revision2.sbatch` | launcher |
| `slurm/leader_full.sbatch` | launcher |
| `slurm/leader_full_after_retry.sbatch` | launcher |
| `slurm/leader_full_revision2.sbatch` | launcher |
| `slurm/leader_pilot.sbatch` | launcher |
| `slurm/leader_pilot_after_retry.sbatch` | launcher |
| `slurm/leader_pilot_revision2.sbatch` | launcher |
| `slurm/native_full.sbatch` | launcher |
| `slurm/native_pilot.sbatch` | launcher |
| `slurm/native_same_file_128G.sbatch` | launcher |
| `slurm/pilot.sbatch` | launcher |
| `slurm/pilot_revision2.sbatch` | launcher |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `identity_comparison.tsv` | `c5a863acb448d698d0a89b440e5c0a42064e55a5475b0de3d793f6f836287090` |
| `technical_agreement_summary.tsv` | `b04ab6e490e0b40c08c3e5fb920c383f6a3995b7a0ce32b96115bc13c53dc401` |
| `technical_agreement_summary.tsv` | `c9d59a24319a1106ab0e1e7895b50133040d78a28faca01918c1895bf7632751` |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

Requests in the launchers (`#SBATCH` lines):

- `slurm/compare.sbatch`: no resource lines
- `slurm/compare_revision2.sbatch`: no resource lines
- `slurm/extract.sbatch`: no resource lines
- `slurm/extract_revision2.sbatch`: no resource lines
- `slurm/leader_full.sbatch`: no resource lines
- `slurm/leader_full_after_retry.sbatch`: no resource lines
- `slurm/leader_full_revision2.sbatch`: no resource lines
- `slurm/leader_pilot.sbatch`: no resource lines
- `slurm/leader_pilot_after_retry.sbatch`: no resource lines
- `slurm/leader_pilot_revision2.sbatch`: no resource lines
- `slurm/native_full.sbatch`: no resource lines
- `slurm/native_pilot.sbatch`: no resource lines
- `slurm/native_same_file_128G.sbatch`: no resource lines
- `slurm/pilot.sbatch`: no resource lines
- `slurm/pilot_revision2.sbatch`: no resource lines

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| none found as top-level constants | | | |
