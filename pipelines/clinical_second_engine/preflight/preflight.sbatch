#!/bin/bash
set -uo pipefail
PY=${FASTALAKE_PYTHON:-python3}
P=${FASTALAKE_OUT_DIR}/preflight
case ${SLURM_ARRAY_TASK_ID} in
 0) ${PIPELINE_DIR}/task.sh $P/campi_repro 0 ;;
 1) ${PIPELINE_DIR}/task.sh $P/campi_repro 1 ;;
 2) if ${PIPELINE_DIR}/task.sh $P/neg_missing_mzml 0; then echo "NEGATIVE CONTROL DID NOT FAIL"; exit 1; else echo "NEGATIVE CONTROL FAILED AS REQUIRED"; fi ;;
 3) if ${PIPELINE_DIR}/task.sh $P/neg_wrong_fasta_hash 0; then echo "NEGATIVE CONTROL DID NOT FAIL"; exit 1; else echo "NEGATIVE CONTROL FAILED AS REQUIRED"; fi ;;
esac
