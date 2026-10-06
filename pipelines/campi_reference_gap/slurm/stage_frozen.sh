#!/bin/bash
# Stage the code into the run directory the way the original run had it: <run>/frozen/.
# The analysis scripts import lib and peptide from <run>/frozen. Shared modules live in A03.
set -euo pipefail
: "${PIPELINE_DIR:?set PIPELINE_DIR to this pipeline directory}" "${FASTALAKE_OUT_DIR:?set FASTALAKE_OUT_DIR to the run directory}"
mkdir -p "${FASTALAKE_OUT_DIR}/frozen"
cp "${PIPELINE_DIR}"/*.py "${PIPELINE_DIR}"/*.params "${FASTALAKE_OUT_DIR}/frozen/" 2>/dev/null || true
SHARED="${PIPELINES_DIR:-${PIPELINE_DIR}/..}/campi_near_miss_audit"
cp "${SHARED}/lib.py" "${SHARED}/peptide.py" "${FASTALAKE_OUT_DIR}/frozen/"
