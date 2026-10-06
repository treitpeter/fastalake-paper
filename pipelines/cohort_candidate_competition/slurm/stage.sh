#!/bin/bash
# Stage the code as the original run had it: <run>/candidates/*.py, with the accumulation
# manifest.json in <run>/ (candidates/worker.py reads BASE/manifest.json, BASE = parent of candidates/).
# The pilot output <run>/candidates/pilot/COMPLETE.json must exist before full_gate.py runs.
set -euo pipefail
: "${PIPELINE_DIR:?set PIPELINE_DIR to this pipeline directory}" "${FASTALAKE_OUT_DIR:?set FASTALAKE_OUT_DIR to the run directory}"
mkdir -p "${FASTALAKE_OUT_DIR}/candidates"
cp "${PIPELINE_DIR}"/candidates/*.py "${PIPELINE_DIR}"/peptide.py "${FASTALAKE_OUT_DIR}/candidates/"
test -f "${FASTALAKE_OUT_DIR}/manifest.json" || { echo "manifest.json missing in FASTALAKE_OUT_DIR" >&2; exit 1; }
