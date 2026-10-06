# healthy_stool_joint_search_q

Internal ID: A22 (previous folder name A22_joint_search_acquisition_q)

## Question

This pipeline produces: out/compare/identity_per_acquisition.tsv (SFig07 a acquisition-q column), identity_summary.tsv, completeness_all.tsv (line 115), microbial_group_pairs.tsv; joint_fileq coverage and iteration-1 counts (32.0%, 32.9%), lines 765-769

Used by: figures/SFig07, text_statistics/supplementary_methods_6_7.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `extract_fileq.py` | D3 / CR-030 follow-up: re-extract HSTOOL LFQ features gated by accepted same-acquisition |
| `compare_fileq.py` | D3 / CR-030 follow-up: completeness, identity and microbial protein-group precision with |
| `slurm/compare.sbatch` | launcher |
| `slurm/full.sbatch` | launcher |
| `slurm/tiny.sbatch` | launcher |
| `recount_joint_q.py` | CR-030: per-acquisition re-estimation of peptide-level target/decoy q inside the joint HSTOOL search. |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `identity_per_acquisition.tsv` | `56c7271c5092a44bf665f2b6ba441da79e463694149a51d1bfb877ec6e28c23a` |
| `identity_summary.tsv` | `f4b9f08e71920b6e08aad21436cc3ed487ebf2ae9049048aef1f168602e0aa7a` |
| `completeness_all.tsv` | `0d73c3788ea060b63c16a068c0b4ba5a74def2f78dc0ef094e8e0e0b3d451f60` |
| `microbial_group_pairs.tsv` | `e8d58907fb5564149344d76b09e6dbe8b3dd3062e83084ed14da083143471329` |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

| job | name | cores | requested memory | peak RSS (GB) | elapsed | state |
|---|---|---|---|---|---|---|
| 6048635 | FL_D3_tiny | 2 | 16G | 1.2 | 00:02:39 | COMPLETED |
| 6048638 | FL_D3_compare | 2 | 16G | 1.0 | 00:01:18 | COMPLETED |
| 6048637_0 | FL_D3_full | 2 | 16G | 1.6 | 00:03:31 | COMPLETED |
| 6048637_1 | FL_D3_full | 2 | 16G | 1.7 | 00:03:29 | COMPLETED |
| 6048637_2 | FL_D3_full | 2 | 16G | 2.4 | 00:03:03 | COMPLETED |

Requests in the launchers (`#SBATCH` lines):

- `slurm/compare.sbatch`: #SBATCH --partition=p.hpcl8,p.hpcl93,p.hpcl92; #SBATCH --cpus-per-task=1; #SBATCH --mem=16G; #SBATCH --time=01:00:00
- `slurm/full.sbatch`: #SBATCH --partition=p.hpcl8,p.hpcl93,p.hpcl92; #SBATCH --cpus-per-task=1; #SBATCH --mem=16G; #SBATCH --time=02:00:00; #SBATCH --array=0-2
- `slurm/tiny.sbatch`: #SBATCH --partition=p.hpcl8,p.hpcl93,p.hpcl92; #SBATCH --cpus-per-task=1; #SBATCH --mem=16G; #SBATCH --time=01:00:00

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `MOD` | `re.compile(r'\[[^\]]*\]')` | no comment in original code | not stated in code; see recount_joint_q.py |
| `TOK` | `re.compile(r'[A-Z](?:\[[^\]]*\])?/\[[^\]]*\]-?')` | no comment in original code | not stated in code; see recount_joint_q.py |
| `PICKED` | `False` | no comment in original code | not stated in code; see recount_joint_q.py |

## Corrected Sage build

The joint searches whose q values this pipeline re-estimates were run with a corrected build of Sage v0.14.7 (LFQ q values assigned by the full target/decoy key; the same change as upstream lazear/sage commit 0e2c7a6). Source pin, patch, lockfile, producing source archive, build instructions and checksums (producing executable sha256 `21241c22a7e886789e3abfb3f129c299ec06080b8ccb504d54906e17689a3420`, build job 5922884, joint search job 5923008) are in `sage_build/` (see `sage_build/README.md`).
