#!/bin/bash
# Resource flags and array ranges of the original submissions of terminal_participant_inference (measured values in README.md).
# Run from the repository root. Required: FASTALAKE_PIPELINES_DIR, FASTALAKE_OUT_DIR, FASTALAKE_PYTHON.
set -euo pipefail
# Original resources: -p p.hpcl8,p.hpcl92 --nice=10000 -c 1 --mem=32G -t 04:00:00 --nodes=1
# Original array ranges: pilot `--array=0%1`, de novo `--array=1%1`, CapScan `--array=2%1`, PREDICT `--array=3%1`
sbatch -p p.hpcl8,p.hpcl92 --nice=10000 -c 1 --mem=32G -t 04:00:00 --nodes=1 --job-name=terminal_participant_inference --array="${ARRAY:?e.g. 0-1%2}" "$FASTALAKE_PIPELINES_DIR/terminal_participant_inference/slurm/run_task.sbatch"
