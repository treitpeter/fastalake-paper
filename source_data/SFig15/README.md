# Source data: Supplementary Figure 15

All three tables were made by `figures/SFig15/prepare_source_data.py`. The reservoir sizes used
for v28d line 193 are read from `source_data/Fig1/acquisition_sizes_726.tsv` (owned by Fig1;
dictionary in `source_data/Fig1/README.md`).

## compute_stage_summary.tsv
sha256: 37c1d94d821d5ba9defb5ff5c8495905ce93209a136952dec8b600b74816498f
rows: 28 (one row per cohort x workflow stage)
produced by: common_resource_accounting (Computing block of the producer, 347-acquisition CapScan roster); recomputed byte for byte from `compute_task_records.tsv` by `figures/SFig15/panels.py:stage_summary`

| column | unit | meaning |
|---|---|---|
| cohort | id | CAMPI, HSTOOL (healthy stool), PREDICT or CAPSCAN |
| stage | id | workflow stage: 02_evidence (shared lookup), 05_persample (acquisition lookup), 06_parsimony, 07_sage, 08_aggregate, 09_entrap_shuffle (shuffle search), 10_fdp_shuffle (shuffle diagnostic) |
| n_tasks | count of Slurm tasks | tasks in this cell; for per-acquisition stages one task per acquisition |
| n_memory_sampled | count of Slurm tasks | tasks with a usable sampled peak memory |
| wall_median_min | minutes | median job execution elapsed time; empty if no task |
| memory_median_GiB | GiB (2^30 bytes) | median sampled peak memory over tasks with a sample; empty if none |
| allocated_CPU_hours | CPU-hours | sum over tasks of allocated CPU time (elapsed x allocated CPUs); not consumed CPU time |
| population_scope | text | which acquisitions the tasks cover: the manuscript subset, or "shared stage as run: <n> acquisitions" |

## compute_task_records.tsv
sha256: fb248054124262d7deff9a9e7653b85d6c4f055c08081cc8b1c29c15ff458c2a
rows: 2,905 (one row per Slurm job or array task)
produced by: common_resource_accounting (sacct records of the August production chain for HSTOOL, PREDICT and CAPSCAN, restricted by array index to the manuscript acquisitions; September broad-reservoir run for CAMPI). The raw acquisition file names were replaced by pseudonymous ids; every other field is as produced.

| column | unit | meaning |
|---|---|---|
| cohort | id | CAMPI, HSTOOL, PREDICT or CAPSCAN |
| stage | id | workflow stage as in the summary, plus 01_gate, 11_eggnog and lake (CAMPI reservoir build), which are not drawn |
| jobid | id | Slurm job id, `<array job>_<index>` for array tasks |
| name | id | Slurm job name |
| state | id | Slurm state (all COMPLETED) |
| elapsed_s | seconds | job execution elapsed time |
| cpu_alloc_s | CPU-seconds | allocated CPU time (elapsed x allocated CPUs); non-CAMPI rows |
| cpu_used_s | CPU-seconds | consumed CPU time; non-CAMPI rows |
| ncpus | count of CPUs | allocated CPUs; non-CAMPI rows |
| reqmem_gb | GB | requested memory; non-CAMPI rows |
| submit, start, end | ISO 8601 time | Slurm submit, start and end times |
| node | id | compute node |
| partition | id | Slurm partition |
| acquisition | id | pseudonymous acquisition id `<cohort>_A<rank>` for per-acquisition tasks of HSTOOL, PREDICT, CAPSCAN; empty for shared stages and for CAMPI (CAMPI tasks are matched by array index) |
| population_scope | text | as in the summary |
| memory_GiB | GiB (2^30 bytes) | sampled peak memory used for the summary; empty if no usable sample |
| wall_minutes | minutes | elapsed_s / 60 |
| allocated_CPU_hours | CPU-hours | cpu_alloc_s / 3600 (non-CAMPI) or elapsed_s x allocated_logical_cpus / 3600 (CAMPI) |
| sampled_peak_gib | GiB | CAMPI only: sampled peak memory from sacct MaxRSS |
| raw_MaxRSS | text | CAMPI only: sacct MaxRSS as reported |
| program_threads | count of threads | CAMPI only: threads given to the program (`auto` for Sage) |
| allocated_logical_cpus | count of CPUs | CAMPI only: allocated logical CPUs |
| cpu_consumed | text (time) | CAMPI only: sacct TotalCPU |

## campi_alphanovo_inference_sacct.tsv
sha256: d007807d2819bc8b6213c9ee05aebd6d02428beb9ab640197a3b75d2eb0d6029
rows: 51 (17 array tasks x 3 rows: job, batch step, extern step)
produced by: common_resource_accounting (H10 compute provenance: `sacct --format=JobID,JobName,State,ElapsedRaw,NCPUS,TotalCPU,MaxRSS,Start,End,NodeList --parsable2 --noheader` for the CAMPI AlphaNovo inference arrays); header added, `|` changed to tab

| column | unit | meaning |
|---|---|---|
| JobID | id | `5450698_<i>` (2 CAMPI faecal acquisitions, not in the 15) or `5450700_<i>` (the 15 CAMPI acquisitions); `.batch`/`.extern` rows are job steps |
| JobName | id | Slurm job or step name |
| State | id | Slurm state |
| ElapsedRaw | seconds | job wall time, including model loading |
| NCPUS | count of CPUs | allocated CPUs |
| TotalCPU | text (time) | consumed CPU time |
| MaxRSS | text | peak resident memory of the step (job rows empty) |
| Start, End | ISO 8601 time | start and end times |
| NodeList | id | GPU node |
