#!/bin/bash
# Healthy-stool common-FASTA identities and completeness at 5-50 and 9-50.
# Original run: srun job 6072769 (4 cores, 1.66 GB peak RSS, 35 s).
set -euo pipefail
: "${FASTALAKE_DATA_DIR:?}" "${FASTALAKE_OUT_DIR:?}"
srun -p p.hpcl92,p.hpcl8,p.hpcl93 --nice=10000 -c4 --mem=4G -t 30 \
    /usr/bin/time -v ${FASTALAKE_PYTHON:-python3} -B "${FASTALAKE_OUT_DIR}/SFig07/work/hstool_recount.py" "${FASTALAKE_OUT_DIR}/SFig07/work/hstool"
