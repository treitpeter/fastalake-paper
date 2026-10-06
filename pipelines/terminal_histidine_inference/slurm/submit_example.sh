#!/bin/bash
# Resource flags and array ranges of the original submissions of terminal_histidine_inference (measured values in README.md).
# Run from the repository root. Required: FASTALAKE_PIPELINES_DIR, FASTALAKE_OUT_DIR, FASTALAKE_PYTHON.
set -euo pipefail
# Original resources: -p p.hpcl8,p.hpcl92 --nice=10000 -c 1 --mem=32G -t 08:00:00 --nodes=1
# Original array ranges: pilot `--array=0%1`, statistics `--array=5%1`, available CapScan `--array=1%1`, available PREDICT `--array=2%1`, full CapScan `--array=3%1`; audit with `FASTALAKE_SCRIPT=audit_inference_revision4.py`
sbatch -p p.hpcl8,p.hpcl92 --nice=10000 -c 1 --mem=32G -t 08:00:00 --nodes=1 --job-name=terminal_histidine_inference --array="${ARRAY:?e.g. 0-1%2}" "$FASTALAKE_PIPELINES_DIR/terminal_histidine_inference/slurm/run_task.sbatch"
