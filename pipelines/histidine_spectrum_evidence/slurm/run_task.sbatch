#!/bin/bash
# histidine_spectrum_evidence array task: one task of the manifest in $FASTALAKE_OUT_DIR/inputs (index = SLURM_ARRAY_TASK_ID).
# Required: FASTALAKE_PIPELINES_DIR, FASTALAKE_OUT_DIR, FASTALAKE_PYTHON. Optional: FASTALAKE_SCRIPT (default search_revision2.py).
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
: "${FASTALAKE_PIPELINES_DIR:?}" "${FASTALAKE_OUT_DIR:?}"
bash "$FASTALAKE_PIPELINES_DIR/capscan_nonspecific_search/slurm/assemble_run_dir.sh" "$FASTALAKE_OUT_DIR" capscan_nonspecific_search denovo_threshold_routes histidine_spectrum_evidence
exec "${FASTALAKE_PYTHON:-python3}" -B "$FASTALAKE_OUT_DIR/frozen/${FASTALAKE_SCRIPT:-search_revision2.py}" "$FASTALAKE_OUT_DIR" "${SLURM_ARRAY_TASK_ID}"
