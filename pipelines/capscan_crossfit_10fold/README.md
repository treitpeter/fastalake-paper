# capscan_crossfit_10fold

Internal ID: A95 (previous folder name A95_crossfit_10fold)

Does the FastaLake gain hold when evidence and search use disjoint spectra (10-fold spectrum cross-fit)?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `analyze_final.py` | call chain of the used outputs |
| `common.py` | call chain of the used outputs |
| `recovery.py` | call chain of the used outputs |
| `run__analyze_final.py` | call chain of the used outputs |
| `sage_lda.py` | call chain of the used outputs |
| `sage_q.py` | call chain of the used outputs |
| `worker.py` | call chain of the used outputs |

Launchers: `slurm/final.sbatch`, `slurm/run.sbatch`.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: Merged A95_spectrum_crossfit (same cross-fit code). crossfit_per_acquisition.tsv: pooled 10-fold cross-fitted FL and MG peptide counts (5-50 and 9-50 residues), | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/Fig4, text_statistics/abstract_discussion, text_statistics/methods, text_statistics/supplementary_methods_6_7. What is used: Merged A95_spectrum_crossfit (same cross-fit code). crossfit_per_acquisition.tsv: pooled 10-fold cross-fitted FL and MG peptide counts (5-50 and 9-50 residues), frozen same-spectrum counts, LDA fallbacks, cross-fitted recovery; 578-acquisition roster plus 347-roster CapScan rerun; Abstract gains (line 35), Methods line 267, validation, one-fold control and pilot fold records

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

| jobs | cores | requested memory | peak RSS | longest task | summed task time | states |
|---|---|---|---|---|---|---|
| 129 tasks | 16, 2, 32, 64, 8 | 16G, 2G, 86G, 90G | 78.9 GB | 1 h 05 min | 55 h 40 min | COMPLETED 129 |

Partitions: p.hpcl8, p.hpcl92, p.hpcl93. Job ids: 6051943, 6065335, 6065336, 6065344.

## Choices

Approach: per-acquisition cross-fit counts for PREDICT and CapScan.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| folds | 10 | deterministic fold = SHA-256 of seed:scan, first 8 bytes, mod 10 | 2-fold odd/even split in A96 |
| MS1 | every MS1 spectrum goes to every output | only MS2 are partitioned | none |

## Known limitations

See `OPEN_ITEMS.md`.

## Cross-fitted PEP sensitivity (added 2026-10-06)

Question. Does the cross-fitted FastaLake gain over the matched metagenome search hold when a PSM posterior error probability gate (PEP <= 0.01) is added to pooled peptide q <= 0.01 (v30 Methods paragraph 125, Supplementary Methods paragraph 326)? Read-only on the 578 x 30 fold searches of this pipeline (iteration4/crossfit/outputs).

Code: `cf_pep.py` (modes `validate <i>` and `run <chunk>`; 20 chunks of 29 roster rows), `cf_analyze.py` (medians and participant-cluster bootstrap), `PEP_MODEL.json` (Sage PEP model settings fixed by the validate run), `export_pep_source_data.py` (added here: writes the Source Data tables, pseudonymising the per-acquisition table). `cf_pep.py` imports `sage_q` and `sage_lda` from this folder (byte-identical to the iteration4/crossfit/design copies used in the run). Launcher `slurm/pep_crossfit.sbatch` (MODE=validate, run as array 0-19, or analyze). Only paths were changed (inputs addressed as LAB_ROOT/..., run directory from FASTALAKE_OUT_DIR); docstrings were added to functions.

Options: baseline = pooled LDA rescoring with picked-peptide q <= 0.01 (must reproduce the frozen pooled accepted sets; recorded per arm as baseline_equal_all); opt_a = baseline and the per-fold Sage posterior_error (log10 PEP) <= -2 on the same PSM, unset value 1.0 excluded; opt_b = baseline and PEP <= 0.01 re-estimated from the pooled discriminant scores with the Sage PSM PEP model (KDE, 1,000 bins, bandwidth 1, monotone), zero yield where pooled rescoring fell back.

Outputs (Source Data, `source_data/supplementary_methods_6_7/`): crossfit_pep_summary.tsv, crossfit_pep_summary_by_region.tsv, crossfit_pep_per_acquisition.tsv.

Resources (sacct): validate 6057983 20 s (8G request); pilot 6057984 23 s; negative control 6057988 1 s (refuses to run without the PEP model); array 6057990, 20 tasks of 1 CPU and 1.5G, 4.0 to 7.9 min each; analysis 6058010 18 s (2G).

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| PEP gate | log10 PEP <= -2 (PEP <= 0.01) on the same rank-one target PSM as the q gate | Sage stores log10(PEP) | PEP <= 0.001 in the non-cross-fitted confidence analysis (pipelines/sage_input_completeness) |
| unset PEP | posterior_error = 1.0 excluded | 1.0 is Sage's default for an unavailable estimate | kept in the q-only baseline |
| PEP model (opt_b) | KDE, 1,000 bins, bandwidth 1, monotone | validated against Sage's own fold outputs: p99 abs diff log10 PEP 3e-5, 0 gate flips over 30 fold searches (job 6057983) | opt_a uses Sage's per-fold PEP directly |
| bootstrap | 10,000 participant-cluster draws, seed 20261004 | as analyze_final.py | none |
