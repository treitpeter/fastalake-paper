#!/bin/bash
# Stage the frozen search code into the run directory the way the original run had it: <run>/frozen/.
# The frozen code is the CAMPI second-engine code, unchanged (same sha256 apart from the
# interpreter and java paths, which the repository copies take from the environment):
# search.py and closed_fragger.params from campi_second_engine, lib.py and peptide.py from
# campi_near_miss_audit.
set -euo pipefail
: "${PIPELINE_DIR:?set PIPELINE_DIR to this pipeline directory}" "${FASTALAKE_OUT_DIR:?set FASTALAKE_OUT_DIR to the run directory}"
PDIR="${PIPELINES_DIR:-${PIPELINE_DIR}/..}"
mkdir -p "${FASTALAKE_OUT_DIR}/frozen"
cp "${PDIR}/campi_second_engine/search.py" "${PDIR}/campi_second_engine/closed_fragger.params" "${FASTALAKE_OUT_DIR}/frozen/"
cp "${PDIR}/campi_near_miss_audit/lib.py" "${PDIR}/campi_near_miss_audit/peptide.py" "${FASTALAKE_OUT_DIR}/frozen/"
