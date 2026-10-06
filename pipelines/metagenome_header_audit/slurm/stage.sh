#!/bin/bash
# Stage the code into the run directory the way the original run had it.
# worker.py locates LAB_ROOT as the 4th parent of its own directory (H.parents[3]) and reads
# manifest.json beside itself, so FASTALAKE_OUT_DIR must be
# ${FASTALAKE_DATA_DIR}/review/v15_integration/analyses/<run name> and hold the run's manifest.json
# (LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/manifest.json; sha256 in RECEIPTS.md).
set -euo pipefail
: "${PIPELINE_DIR:?set PIPELINE_DIR to this pipeline directory}" "${FASTALAKE_OUT_DIR:?set FASTALAKE_OUT_DIR to the run directory}"
mkdir -p "${FASTALAKE_OUT_DIR}"
cp "${PIPELINE_DIR}"/*.py "${FASTALAKE_OUT_DIR}/"
test -f "${FASTALAKE_OUT_DIR}/manifest.json" || { echo "manifest.json missing in FASTALAKE_OUT_DIR" >&2; exit 1; }
