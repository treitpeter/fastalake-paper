#!/bin/bash
# Resource flags and array ranges of the original submissions of capscan_nonspecific_search (measured values in README.md).
# Run from the repository root. Required: FASTALAKE_PIPELINES_DIR, FASTALAKE_OUT_DIR, FASTALAKE_PYTHON.
set -euo pipefail
# Original resources: -p p.hpcl93 --nice=10000 -c 8 --mem=640G -t 08:00:00 --nodes=1
# Original array ranges: pilots `--array=0-1%2` (task 0 capsule, task 1 stool), cohort `--array=2-330%10`
sbatch -p p.hpcl93 --nice=10000 -c 8 --mem=640G -t 08:00:00 --nodes=1 --job-name=capscan_nonspecific_search --array="${ARRAY:?e.g. 0-1%2}" "$FASTALAKE_PIPELINES_DIR/capscan_nonspecific_search/slurm/run_task.sbatch"
