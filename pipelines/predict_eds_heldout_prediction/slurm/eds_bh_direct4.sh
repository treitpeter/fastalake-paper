#!/bin/bash
# Four-test BH recompute (original run: Slurm job 6074462, 2 cores, 2G).
# Required: FASTALAKE_NESTED_REFIT_DIR (directory with model_metrics.tsv), FASTALAKE_OUT_DIR.
# Launch: srun -p p.hpcl92,p.hpcl8,p.hpcl93 --nice=10000 -c2 --mem=2G bash slurm/eds_bh_direct4.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "${FASTALAKE_OUT_DIR:?set FASTALAKE_OUT_DIR}"
cd "${FASTALAKE_OUT_DIR}"
"${FASTALAKE_PYTHON:-python}" -B "${HERE}/eds_bh_direct4.py" "${FASTALAKE_NESTED_REFIT_DIR:?set FASTALAKE_NESTED_REFIT_DIR}"
