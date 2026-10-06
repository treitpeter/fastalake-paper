#!/bin/bash
# Stage the code as the original run had it: <final>/discussion_counts/*.py and
# <final>/discussion_counts/spectra_db_ids/*.py, where <final> is the review "final" directory
# (common.py reads FINAL = parent of its own directory, e.g. FINAL/public347/A08_347 and
# FINAL/clustering_check). FASTALAKE_OUT_DIR is that discussion_counts directory;
# FASTALAKE_DATA_DIR is LAB_ROOT.
set -euo pipefail
: "${PIPELINE_DIR:?}" "${FASTALAKE_OUT_DIR:?}" "${FASTALAKE_DATA_DIR:?}"
mkdir -p "${FASTALAKE_OUT_DIR}/spectra_db_ids" "${FASTALAKE_OUT_DIR}/logs"
cp "${PIPELINE_DIR}"/*.py "${FASTALAKE_OUT_DIR}/"
cp "${PIPELINE_DIR}"/spectra_db_ids/*.py "${FASTALAKE_OUT_DIR}/spectra_db_ids/"
