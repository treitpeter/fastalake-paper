#!/bin/bash
set -euo pipefail
# sbatch copies this file into its spool; --chdir supplies the experiment path.
: "${SLURM_JOB_ID:?This launcher must run under Slurm}"
EXPERIMENT_DIR=$PWD
export PYTHONDONTWRITEBYTECODE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec ${FASTALAKE_PYTHON:-python3} "$EXPERIMENT_DIR/driver.py" "$@"
