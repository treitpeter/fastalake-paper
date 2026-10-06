#!/bin/bash
#SBATCH -p p.hpcl92,p.hpcl8,p.hpcl93
#SBATCH --nice=10000
#SBATCH -c 4
#SBATCH --mem=3G
#SBATCH -t 00:15:00
#SBATCH -J crossomics_distance_sensitivity
# Parametrised by environment variables:
#   FASTALAKE_REPO      root of this repository
#   FASTALAKE_DATA_DIR  frozen inputs, laid out as in RECEIPTS.md (original relative paths)
#   FASTALAKE_OUT_DIR   writable work and output directory
#   FASTALAKE_PYTHON    interpreter with numpy, pandas, scipy, statsmodels (default: python)
set -euo pipefail
: "${FASTALAKE_DATA_DIR:?}" "${FASTALAKE_OUT_DIR:?}"
PY="${FASTALAKE_PYTHON:-python}"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg
PIPES="${FASTALAKE_REPO:?}/pipelines"
HERE="$PIPES/crossomics_distance_sensitivity"
# Production PREDICT rows (needs A84 cross_refit and A86 native sources under the data layout)
W="$FASTALAKE_OUT_DIR/B"; mkdir -p "$W"
cp "$HERE"/distance_sensitivity.py "$HERE"/common.py "$HERE"/experiment.py "$W"/
cd "$W"
"$PY" -u distance_sensitivity.py baseline
"$PY" -u distance_sensitivity.py refit
# 347 rerun (CapScan rows)
R="$FASTALAKE_OUT_DIR/A4"; mkdir -p "$R/code"
cp "$HERE"/*.py "$PIPES/native_crossomics_profiles/native_pyramids.py" "$R/code"/
cd "$R"
"$PY" -u code/make_configs.py
"$PY" -u code/distance_capscan347.py control331.json
"$PY" -u code/distance_capscan347.py fixed347.json
"$PY" -u code/collect.py
