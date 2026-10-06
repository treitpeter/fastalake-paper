#!/bin/bash
# PREDICT curation-arm accepted peptides at 5-50 and 9-50 (Supplementary Fig. 8).
# Writes FASTALAKE_OUT_DIR/SFig08/curation_counts_5_9.tsv and checks the 5-50 counts against
# source_data/SFig08/PREDICT_curation_counts.tsv. The argument is the number of worker processes.
set -euo pipefail
: "${FASTALAKE_DATA_DIR:?}" "${FASTALAKE_OUT_DIR:?}"
srun -p p.hpcl92,p.hpcl8,p.hpcl93 --nice=10000 -c8 --mem=16G -t 60 \
    ${FASTALAKE_PYTHON:-python3} -B "${FASTALAKE_OUT_DIR}/SFig08/build/recount_curation.py" 8
