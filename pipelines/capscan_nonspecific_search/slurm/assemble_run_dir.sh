#!/bin/bash
# Assemble RUN_DIR/frozen from the code of the named pipelines (copies all *.py; later pipelines win on name clashes).
# usage: assemble_run_dir.sh RUN_DIR PIPELINE_ID [PIPELINE_ID ...]
# Requires FASTALAKE_PIPELINES_DIR.
set -euo pipefail
RUN="${1:?run dir}"; shift
: "${FASTALAKE_PIPELINES_DIR:?}"
mkdir -p "$RUN/frozen"
for id in "$@"; do
  cp -f "$FASTALAKE_PIPELINES_DIR/$id"/*.py "$RUN/frozen/"
done
