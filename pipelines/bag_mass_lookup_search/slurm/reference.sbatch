#!/bin/bash
set -euo pipefail
export PYTHONDONTWRITEBYTECODE=1
exec ${FASTALAKE_DATA_DIR}/python -B -u campi_reference.py
