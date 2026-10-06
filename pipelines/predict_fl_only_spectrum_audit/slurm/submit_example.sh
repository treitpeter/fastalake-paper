#!/bin/bash
# Resource flags and array ranges of the original submissions of predict_fl_only_spectrum_audit (measured values in README.md).
# Run from the repository root. Required: FASTALAKE_PIPELINES_DIR, FASTALAKE_OUT_DIR, FASTALAKE_PYTHON.
set -euo pipefail
# Original resources: -p p.hpcl8 --nice=10000 -c 16 --mem=52500M -t 04:00:00 --nodes=1
# Original array ranges: `--array=0-295%30` (task index selects a line of `inputs/samples_296.txt`)
sbatch -p p.hpcl8 --nice=10000 -c 16 --mem=52500M -t 04:00:00 --nodes=1 --job-name=predict_fl_only_spectrum_audit --array="${ARRAY:?e.g. 0-1%2}" "$FASTALAKE_PIPELINES_DIR/predict_fl_only_spectrum_audit/slurm/run_task.sbatch"
