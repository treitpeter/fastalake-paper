#!/bin/bash
# Resource flags and array ranges of the original submissions of denovo_threshold_routes (measured values in README.md).
# Run from the repository root. Required: FASTALAKE_PIPELINES_DIR, FASTALAKE_OUT_DIR, FASTALAKE_PYTHON.
set -euo pipefail
# Original resources: -p p.hpcl8,p.hpcl92 --nice=10000 -c 1 --mem=16G -t 01:00:00 --nodes=1
# Original array ranges: pilot `--array=0-2%1`, full `--array=3-636%1`
sbatch -p p.hpcl8,p.hpcl92 --nice=10000 -c 1 --mem=16G -t 01:00:00 --nodes=1 --job-name=denovo_threshold_routes --array="${ARRAY:?e.g. 0-1%2}" "$FASTALAKE_PIPELINES_DIR/denovo_threshold_routes/slurm/run_task.sbatch"
