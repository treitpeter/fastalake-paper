#!/bin/bash
# One task = one acquisition x arm. Usage: task.sh RUN_DIR INDEX
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
PY=${FASTALAKE_PYTHON:-python3}
RUN=$1; IDX=$2
CODE=$(dirname "$(readlink -f "$0")")
ID=$($PY -c "import json,sys;print(json.load(open(sys.argv[1]))['tasks'][int(sys.argv[2])]['id'])" "$RUN/inputs/manifest.json" "$IDX")
if [ -f "$RUN/results/$ID/RECEIPT.json" ]; then echo "SKIP receipt exists $ID"; exit 0; fi
$PY -B "$RUN/frozen/search.py" "$RUN" "$IDX"
$PY -B "$CODE/receipt.py" "$RUN" "$IDX"
