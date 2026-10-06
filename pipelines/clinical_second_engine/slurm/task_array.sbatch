#!/bin/bash
# manifest task index = SLURM_ARRAY_TASK_ID + SE_OFFSET (MaxArraySize 1001)
set -euo pipefail
exec ${PIPELINE_DIR}/task.sh ${FASTALAKE_OUT_DIR} "$(( SLURM_ARRAY_TASK_ID + ${SE_OFFSET:-0} ))"
