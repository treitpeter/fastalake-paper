#!/bin/bash
set -euo pipefail
exec bash ${PIPELINE_DIR}/slurm/job.sh sample --index "$SLURM_ARRAY_TASK_ID" --seed "${1:?seed}"
