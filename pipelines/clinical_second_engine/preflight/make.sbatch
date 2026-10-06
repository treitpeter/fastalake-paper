#!/bin/bash
set -euo pipefail
${FASTALAKE_PYTHON:-python3} -B ${FASTALAKE_OUT_DIR}/preflight/make_preflight.py ${FASTALAKE_OUT_DIR}
