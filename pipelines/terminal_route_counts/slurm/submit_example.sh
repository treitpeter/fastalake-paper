#!/bin/bash
# Resource flags and array ranges of the original submissions of terminal_route_counts (measured values in README.md).
# Run from the repository root. Required: FASTALAKE_PIPELINES_DIR, FASTALAKE_OUT_DIR, FASTALAKE_PYTHON.
set -euo pipefail
# Original resources: -p p.hpcl8,p.hpcl92 --nice=10000 -c 1 --mem=16G -t 02:00:00 --nodes=1
# Original array ranges: pilot `--array=0-1%2`, CapScan `--array=2-331%1`, PREDICT `--array=332-636%1`; collectors via `FASTALAKE_SCRIPT=collector.py` with `--mem=16G -t 03:00:00`
sbatch -p p.hpcl8,p.hpcl92 --nice=10000 -c 1 --mem=16G -t 02:00:00 --nodes=1 --job-name=terminal_route_counts --array="${ARRAY:?e.g. 0-1%2}" "$FASTALAKE_PIPELINES_DIR/terminal_route_counts/slurm/run_task.sbatch"
