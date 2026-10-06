#!/bin/bash
# Resource flags and array ranges of the original submissions of histidine_boundary_localisation (measured values in README.md).
# Run from the repository root. Required: FASTALAKE_PIPELINES_DIR, FASTALAKE_OUT_DIR, FASTALAKE_PYTHON.
set -euo pipefail
# Original resources: -p p.hpcl8,p.hpcl92 --nice=10000 -c 1 --mem=32G -t 02:00:00 --nodes=1
# Original array ranges: available pilot `--array=0-3%1`, available full `--array=4-636%1`, non-specific pilot `--array=637-640%1`, CapScan `--array=641-<last>%1` (328 tasks completed), PREDICT repair `--array=332-636%1` with `FASTALAKE_SCRIPT=predict_index_wrapper.py`
sbatch -p p.hpcl8,p.hpcl92 --nice=10000 -c 1 --mem=32G -t 02:00:00 --nodes=1 --job-name=histidine_boundary_localisation --array="${ARRAY:?e.g. 0-1%2}" "$FASTALAKE_PIPELINES_DIR/histidine_boundary_localisation/slurm/run_task.sbatch"
