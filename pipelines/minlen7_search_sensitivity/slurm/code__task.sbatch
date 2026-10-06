#!/bin/bash
#SBATCH --cpus-per-task=16
#SBATCH --nodes=1
set -euo pipefail
export OMP_NUM_THREADS=16
IDX=$(( ${TASK_OFFSET:-0} + ${SLURM_ARRAY_TASK_ID} ))
if [ -n "${IDX_LIST:-}" ]; then IDX=$(sed -n "$(( SLURM_ARRAY_TASK_ID + 1 ))p" "$IDX_LIST"); fi
echo "task $IDX on $(hostname) $(date -Is)"
exec /mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python "$RUN_DIR/code/run_task.py" "$IDX" ${PILOT_FLAG:-}
