#!/bin/bash
# Resource flags and array ranges of the original submissions of predict_truncated_mzml_repair (measured values in README.md).
# Run from the repository root. Required: FASTALAKE_PIPELINES_DIR, FASTALAKE_OUT_DIR, FASTALAKE_PYTHON.
set -euo pipefail
# Original resources: -p p.hpcl8,p.hpcl92 --nice=10000 -c 8 --mem=480G -t 04:00:00 --nodes=1
# Original array ranges: `--array=0-1` for the search; the reconversion is a separate array over `raw_recovery/manifest.json` tasks
sbatch -p p.hpcl8,p.hpcl92 --nice=10000 -c 8 --mem=480G -t 04:00:00 --nodes=1 --job-name=predict_truncated_mzml_repair --array="${ARRAY:?e.g. 0-1%2}" "$FASTALAKE_PIPELINES_DIR/predict_truncated_mzml_repair/slurm/run_task.sbatch"
