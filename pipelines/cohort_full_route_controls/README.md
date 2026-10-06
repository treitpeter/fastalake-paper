# cohort_full_route_controls

Internal ID: A91 (previous folder name A91_full_route_controls)

## Question

This pipeline produces: tables/<cohort>/{per_acquisition,stage_control_counts,baseline_comparison,summary,control_discoveries}.tsv for STab2 Cohort_overview, Per_acquisition, Stage_controls sheets

Used by: tables/STab2.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `full_route_entrapment.py` | Execute a reservoir-entry control through the unclustered production recipe. |
| `prepare_experiment.py` | Freeze a full-route experiment after successful control generation. |
| `make_controls.cpp` | include <algorithm> |
| `scan_overlap.cpp` | include <algorithm> |
| `slurm/run_stage.sbatch` | launcher |
| `collect_cohort.py` | Require a complete cohort, collect control traces and scan target overlap. |
| `slurm/collect.sbatch` | launcher |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `summary.tsv` | `233f8509033ec03ecebfc29b3bbcd77aeafb98274637395c1747cb4b5a11cc56` |
| `summary.tsv` | `a27905b90a74f0619b91d0a8395247c87b8caa6508bfd84bcb95106217d74ef3` |
| `summary.tsv` | `fcfe131ec3c57d5a8d08bf3e10b59969e88ee79ef944a40c4806b3518203a92e` |
| `summary.tsv` | `90aa4037a641b85bc578de640fb21033fa1adb04424a5f1fbbef9a703310b2d8` |
| `summary.tsv` | `e93983b04916fdbd3a4cc16fbc765b6dba014246f9534df1981a4b26d628be37` |
| `summary.tsv` | `784e836457068051aba9b730866de91750c8e79f1bf8feed5ae7fc2c3817d2fe` |
| `summary.tsv` | `90e15ff46cba3d7dbb22b36d862c5fee74069aeb00248e089996133657d1795c` |
| `summary.tsv` | `4880a00f99950cfb6eb4c0a47ae39ee90f86304cb008c22172e55be66208ee9a` |
| `summary.tsv` | `03e30fe0b9c70ddc97f400611b9b062d014cd30e5edeafbc86660f870d449888` |
| `summary.tsv` | `f69775e53d8c5558177aa5b7148868a319eab1fb84e6f8904049ed0ef8627c55` |
| `summary.tsv` | `11a1f99733ab026c6aa41d7b85e406bd18b376c8337f607cf3a250cfc917f990` |
| `summary.tsv` | `5e975735277956ca53905ef370b0ed9fdacd279ecb2df99179a79cb2f199e529` |
| `per_acquisition.tsv` | `7fbbd87610c1a5298f80a7e3a004470e2cda563da260d4f6d8e0e2a8916f761d` |
| `per_acquisition.tsv` | `c182e2b1cce6748a8479762f7737fc10eb52f850533d91f27ae6d46d20e434c9` |
| `per_acquisition.tsv` | `378e2e585174abb1de0d0de4f41679c7276ab6bb47e6090ca5ba90f6677a1af6` |
| `per_acquisition.tsv` | `8bd413c252992de908822dceae425229eb18488dac01e4b89939e524cd4541d5` |
| `stage_control_counts.tsv` | `fbce9906dd3d2d6570b9f29adf3a602bdd7f11eb43637b7c4a7aba622b60f925` |
| `stage_control_counts.tsv` | `893029d1a66256745d68af3288f9bcc9cd7ceda0cebe78fa1f970f6d9d0e2922` |
| `stage_control_counts.tsv` | `49ab806664aa642aac80ca1df5bb18a9d374a696c17396dc599158b07072b203` |
| `stage_control_counts.tsv` | `61453e81dfee1380eaf1294494524dbaf804fbe328fd92af2e71295895922a7f` |
| `baseline_comparison.tsv` | `f963ed20745e7564117349001b7e5990f3e9b282218f1e9dfdad6fc87ac855f7` |
| `baseline_comparison.tsv` | `eacc2fa2950f5671898588699f56cbf06d002fd1da19cacdefe828016238dbc8` |
| `baseline_comparison.tsv` | `ad00f14eabdfbe7f7a9b9b1d57bc44747f092ec90b946632b66b2330a7f4a14d` |
| `baseline_comparison.tsv` | `73a0d370fd8a55e2547a39f2b15d6c5506fec934f87d2466c5e7ff47032a0d91` |
| `control_discoveries.tsv` | `155c43af8f931cb0d15b98c9c385a4ab3cd9200e981eebb3739dacf0650cbe42` |
| `control_discoveries.tsv` | `8987ef72a5173d40abf3321e6bc30d0eceb1bb6602f93e83ba02231603138c21` |
| `control_discoveries.tsv` | `7dafc42ccc5d25e6740c294da9fcce13b97207844ffd040dde8d40256f8b36dc` |
| `control_discoveries.tsv` | `62110ad8b152b282fdbb0cfca3a74cf0a470d6622f6db2d5860265fe360a9209` |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

No job id could be traced for this analysis; see OPEN_ITEMS.md.

Requests in the launchers (`#SBATCH` lines):

- `slurm/run_stage.sbatch`: #SBATCH --time=02:00:00
- `slurm/collect.sbatch`: #SBATCH --partition=p.hpcl93; #SBATCH --cpus-per-task=16; #SBATCH --mem=16G; #SBATCH --time=02:00:00

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `ARCH` | `F/'FIGURATOR/Supplementary/ENTRAPMENT/rerun_2026-09-03'` | no comment in original code | not stated in code; see prepare_experiment.py |
