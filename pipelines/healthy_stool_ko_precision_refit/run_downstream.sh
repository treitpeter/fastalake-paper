#!/bin/bash
# Run the frozen biology_refit drivers (byte-identical copies) in the shadow for one variant.
set -uo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg
V=$1
W=${FASTALAKE_DATA_DIR}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration1/fix_eggnog_native
S=$W/downstream/$V/B; PY=${FASTALAKE_DATA_DIR}/python
for a in PREDICT_noclust CAPSCAN_noclust HSTOOL_noclust; do test -f $W/profiles/$V/$a/COMPLETE.json || { echo "missing profiles $a"; exit 1; }; done
cd $S
for step in "crossomics.py refit" "capscan_functions.py" "eds_biochemistry.py refit" "nested_eds.py refit" "individual_functions.py"; do
  echo "=== $step $(date)"; $PY -u $S/$step && echo "STEP_PASS $step" || echo "STEP_FAIL $step"
done
