# minlen7_search_sensitivity

Internal ID: A09 (previous folder name A09_minlen7_sensitivity)

## Question

This pipeline produces: Per-acquisition peptide counts by length population for all repeated min-length-7 searches (SFig04 a-c), CAMPI member counts (d), fusion-control summaries (e), grouping; minlen7_benchmark_per_acquisition.tsv for 14.5%/45.0% gains and 85.4%/84.0% recovery, line 777

Used by: figures/SFig04, text_statistics/supplementary_methods_6_7.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `prepare.py` | Build isolated min_len=7 Sage configs for M09 wave 1. Writes configs/, tasks.tsv, PREPARE.json. |
| `prepare_addendum1.py` | Addendum 1: five LFQ-on repaired PREDICT FL_cohortMG searches (indices 1985-1989), appended; existing rows untouched. |
| `code__run_task.py` | Run one M09 search. Usage: run_task.py <index> [--pilot]. Writes a receipt JSON; exits non-zero on any failed check. |
| `slurm/code__task.sbatch` | launcher |
| `validate_pilot.py` | Gate: every pilot receipt passes. Also records baseline comparison (not a gate) and Slurm MaxRSS. |
| `collect_fusion.py` | no module docstring in original |
| `collect_healthy_claude.py` | Healthy-stool own / common-FASTA / joint strategies at 5-50 (published) and fresh >=7. |
| `collect_peptides.py` | no module docstring in original |
| `dispatch.py` | no module docstring in original |
| `joint_lfq_audit.py` | no module docstring in original |
| `prepare_wave2.py` | no module docstring in original |
| `revalidate_fieldlimit.py` | Revalidate wave-2 tasks whose receipt FAILED only because csv's default field limit (131072) was exceeded while |
| `revalidate_joint_pilot.py` | no module docstring in original |
| `run_search.py` | no module docstring in original |
| `summarize_peptides.py` | no module docstring in original |
| `validate_pilots.py` | no module docstring in original |
| `member_counts.py` | CAMPI member-specific peptide counts for baseline 5-50, nested >=7, fresh >=7, nested >=9 and fresh7_nested9. |
| `run_fusion.py` | M09 fusion control on the CURRENT broad-reservoir CAMPI fusion (campi_fusion_rerun, 375,333,654 entries; manuscript P37). |
| `group_fresh.py` | M09 grouping: run the EXECUTED production grouping (study.py sha256 5f69c457…) unchanged on fresh min_len-7 searches. |
| `bench.sh` | Memory benchmark = real pilot searches (frozen run_task.py --pilot, writes run/pilot/) under /usr/bin/time -v. |
| `collect.py` | A2 collector: writes RESULTS.md (old 331 vs new 347 for every printed number of the min-length-7 sensitivity). |
| `slurm/collect.sbatch` | launcher |
| `slurm/grouping.sbatch` | usage: sbatch ... grouping.sbatch <driver.py>   (run331/grouping/group_fresh.py or run/grouping/group_fresh_347.py) |
| `preflight.py` | A2 preflight. usage: preflight.py <searches|grouping331|grouping347> [--negative] |
| `preflight_all.sh` | A2 preflight on a compute node: control, build, positive checks, negative controls (must exit 1). |
| `collect16.py` | Supp. Fig. 4c 'CapScan FL' (FL_public) row, 331 vs 347, and the search count. |
| `collect_peptides.M09.py` | no module docstring in original |
| `prepare16.py` | min_len=7 FL_public configs for the 16 participant-15 A08_347 searches, by the M09 prepare.py transformation |
| `2_code__run_task.py` | Run one M09 search. Usage: run_task.py <index> [--pilot]. Writes a receipt JSON; exits non-zero on any failed check. |
| `slurm/2_code__task.sbatch` | launcher |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

See RECEIPTS.md.

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

Requests in the launchers (`#SBATCH` lines):

- `slurm/code__task.sbatch`: #SBATCH --cpus-per-task=16
- `slurm/collect.sbatch`: #SBATCH --cpus-per-task=1
- `slurm/grouping.sbatch`: #SBATCH --cpus-per-task=2
- `slurm/2_code__task.sbatch`: #SBATCH --cpus-per-task=16

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `OUT` | `F + "/canonical/publication_data/arms/M09_minlen7_2026-09-24"` | no comment in original code | not stated in code; see prepare.py |
| `NEW_MIN` | `7` | no comment in original code | not stated in code; see prepare.py |
| `OUT` | `F + "/canonical/publication_data/arms/M09_minlen7_2026-09-24"` | no comment in original code | not stated in code; see prepare_addendum1.py |
| `SAGE_SHA` | `"d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625"` | no comment in original code | not stated in code; see code__run_task.py |
| `ROOT` | `Path.cwd().resolve();O=Path(__file__).resolve().parent.parent` | no comment in original code | not stated in code; see collect_fusion.py |
| `POPS` | `["baseline5", "nested7", "fresh7", "nested9", "fresh7_nested9"]` | no comment in original code | not stated in code; see member_counts.py |
| `BASE` | `F + "/simplified_grouping_results/manuscript_cohorts_2026-09-10_2200"` | no comment in original code | not stated in code; see group_fresh.py |
| `STUDY` | `BASE + "/code/fasta_lake/study.py"` | no comment in original code | not stated in code; see group_fresh.py |
| `SAGE` | `F + "/MM_review_analysis/A13_4_7_predict_unspecific/bin/sage"` | no comment in original code | not stated in code; see preflight.py |
| `SAGE_SHA` | `"d3820543a31bfa2a556e04f91719204ce0ac6f34dcca6670314a13df87064625"` | no comment in original code | not stated in code; see 2_code__run_task.py |
