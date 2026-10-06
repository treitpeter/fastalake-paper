#!/bin/bash
# Resource flags and array ranges of the original submissions of predict_histidine_case_audit (measured values in README.md).
# Run from the repository root. Required: FASTALAKE_PIPELINES_DIR, FASTALAKE_OUT_DIR, FASTALAKE_PYTHON.
set -euo pipefail
# Original resources: -p p.hpcl8,p.hpcl92 --nice=10000 -c 1 --mem=16G -t 01:00:00 --nodes=1
# Original array ranges: no array; `--phase pilot` then `--phase full`
sbatch -p p.hpcl8,p.hpcl92 --nice=10000 -c 1 --mem=16G -t 01:00:00 --nodes=1 --job-name=predict_histidine_case_audit "$FASTALAKE_PIPELINES_DIR/predict_histidine_case_audit/slurm/run_task.sbatch" pilot   # then: full
