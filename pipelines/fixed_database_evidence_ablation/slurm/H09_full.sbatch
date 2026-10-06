#!/bin/bash
# Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_WORK_DIR (working directory).
# Python code addresses inputs as LAB_ROOT/<relative path>, so link it into the working directory.
cd "${FASTALAKE_WORK_DIR:?set FASTALAKE_WORK_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
${FASTALAKE_LAB_ROOT}/reviews/quantification_mechanisms_2026-09-13/apt_env/bin/python -B ${FASTALAKE_LAB_ROOT}/Hartout_edits/H09_evidence_ablation/frozen/h09_evidence_ablation.py ${FASTALAKE_LAB_ROOT}/Hartout_edits/H09_evidence_ablation full "${SLURM_ARRAY_TASK_ID}"
