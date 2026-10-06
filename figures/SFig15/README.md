# Supplementary Figure 15

Shows the measured computational cost of the CPU production workflow: how long each stage
ran, how much memory it used and how many CPU-hours it was allocated, per cohort. The four
panels are heatmaps of seven workflow stages (rows) by four cohorts (columns), computed from
the Slurm accounting records of the production runs. The figure supports the statement that
the computational requirements are modest (v28d lines 193 and 299).

## Panels

The producer draws all four panels in one matplotlib figure; the notebook draws them together
in the section `Panels a-d: stage summary and drawing` and lists each panel's printed values in
its own section.

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| a | median task execution minutes per stage and cohort | lookup 8-12 min per cohort (193; v30f no longer prints the selection 1-3 min and Sage under 1 min of v28d); CAMPI 11, 0.48, 0.63, 0.23 min (299) | per-acquisition stages: 15 CAMPI, 58 healthy-stool, 306 PREDICT, 347 CapScan acquisitions (one Slurm task each); shared stages: one task per cohort | none (descriptive) | Panel a | `source_data/SFig15/compute_stage_summary.tsv` (recomputed from `compute_task_records.tsv`) | common_resource_accounting (not in `pipelines/`) |
| b | median usable sampled peak memory (GiB) | legend (469) | as a; memory-sampled tasks per cell in panel d | none (descriptive) | Panel b | as a | as a |
| c | sum of allocated CPU-hours | legend (469) | as a | none (descriptive) | Panel c | as a | as a |
| d | tasks with a usable memory sample / tasks | legend (469); scopes in Methods (903) | as a | none (descriptive) | Panel d | as a | as a |

Text numbers that are not printed in the figure but come from the same records: the CAMPI
reservoir build time and threads (`compute_task_records.tsv`, stage `lake`), the CAMPI
AlphaNovo inference median (`campi_alphanovo_inference_sacct.tsv`) and the reservoir size range
(`source_data/Fig1/acquisition_sizes_726.tsv`, owned by Fig1). The unit has no statistical
tests; `STATS.tsv` has its header only.

## v30 composite

From manuscript v30 on, the figure is the house-style redraw of the same values (`panels.py:figure_v30`, notebook section `v30 composite (house style)`; ASSEMBLY.md, second part). Panels and lettering are unchanged; only type, colours, size and the placement of panel titles changed. The function compares every plotted array with the producer drawing and raises on any difference. The regenerated PNG is byte-identical to the v30 image (`checks/SFig15_v30.json`).

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `STAGES` | 02_evidence, 05_persample, 06_parsimony, 07_sage, 08_aggregate, 09_entrap_shuffle, 10_fdp_shuffle | the CPU production stages (Methods, v28d line 903); gate (01), eggNOG (11) and reservoir build (`lake`) are not drawn | none in the figure; the CAMPI reservoir build is quoted in the text (line 299) |
| `PER_ACQUISITION_STAGES` | 05, 06, 07, 09 | stages run as one array task per acquisition, restricted to the manuscript acquisitions by array index (Methods, line 903) | shared stages keep their 97- and 368-acquisition scopes (line 903) |
| `N_STAGE` | 07_sage | its task count gives the acquisition n in the x tick labels | any per-acquisition stage gives the same n (checked by `stat_acquisition_counts`) |
| `PANEL_COLUMNS` | wall_median_min, memory_median_GiB, allocated_CPU_hours, n_memory_sampled | the four summaries in the legend (line 469) | none |
| shading | log10(1 + value), `CMAP` = Blues, white text above shade `SHADE_TEXT_THRESHOLD` = 0.62 | producer | none |
| `SELECTION_STAGES` | 05_persample + 06_parsimony, summed per acquisition, then the median | reading of "database-selection time" (line 193) that gives the printed 1-3 min; no original code computed it | sum of the two stage medians gives 1.12-3.14 min, the same printed range (see `OPEN_ITEMS.md`) |
| `CAMPI_INFERENCE_ARRAY` | 5450700 (15 tasks) | the CAMPI inference array of the 15 CAMPI acquisitions; job wall time, so it includes model loading | the per-acquisition execution spans in the prediction logs give the same median to 1 s |
| `V30_SIZE_MM`, `V30_CELL_PT`, `V30_TICK_LEN`, `V30_H_PAD`, `V30_W_PAD`, `V30_SPACE` | 160 x 155.5 mm; 7 pt; 2.6 pt; constrained-layout pads | v30 house style at the placed size; shading white to #0072B2 on log10(1 + value), cell text white or ink by higher contrast | v28f producer: Blues, 0.62 shade threshold |
| `V30_LETTER_PT`, `V30_DPI` | 7 pt, 300 dpi | house panel letters; PNG = PDF rasterised at 300 dpi | none |
| `RC`, `FIGSIZE_IN`, `RASTER_DPI` | producer rcParams (DejaVu Sans, 8 pt), 7.1 x 6.9 in, 300 dpi | producer | none |

## Known limitations

- CPU-hours are allocated CPU-hours (elapsed time x allocated CPUs, Slurm `CPUTimeRAW`), not
  consumed CPU time (`cpu_used_s` in the task records), as the legend and Methods state (v28d
  lines 469, 903). Memory is the sampled peak RSS in GiB (2^30 bytes).
- All timings are Slurm job execution times. They exclude queueing, raw-file conversion, de
  novo inference (except the separate CAMPI AlphaNovo median) and downstream analyses.
- Memory is sampled by Slurm accounting; short tasks often have no usable sample (panel d), and
  CAMPI aggregation has none.
- CAMPI rows come from the September broad-reservoir production run; the other cohorts come
  from the August production chain. Shared healthy-stool and CapScan stages ran over 97 and 368
  acquisitions, more than the manuscript subsets.
- The `acquisition` column of the task records is a pseudonymous id; raw acquisition file names
  are not shipped (`OPEN_ITEMS.md`).
- Node hardware (512 GB, AMD EPYC, NVIDIA L40S) is not in any shipped table.
- Definitional constants in the legend (2^30 bytes, log10(1 + value)) have no NUMBERS.tsv row.
