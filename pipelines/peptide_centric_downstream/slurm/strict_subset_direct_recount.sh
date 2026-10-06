#!/bin/bash
# Fig. 5f direct-assignment recount (original run: Slurm job 6074394, 2 cores, 8G).
# Required: FASTALAKE_FIG5EF_RUN_DIR (fig5ef_run.py run directory with results/),
#           FASTALAKE_OUT_DIR (where recount.json and the per-participant tables are written).
# Launch: srun -p p.hpcl92,p.hpcl8,p.hpcl93 --nice=10000 -c2 --mem=8G bash slurm/strict_subset_direct_recount.sh
set -euo pipefail
HERE="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "${FASTALAKE_OUT_DIR:?set FASTALAKE_OUT_DIR}"
cd "${FASTALAKE_OUT_DIR}"
"${FASTALAKE_PYTHON:-python}" -B "${HERE}/strict_subset_direct_recount.py" "${FASTALAKE_FIG5EF_RUN_DIR:?set FASTALAKE_FIG5EF_RUN_DIR}"
