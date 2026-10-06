#!/bin/bash
#SBATCH -p p.hpcl92,p.hpcl8,p.hpcl93
#SBATCH --nice=10000
#SBATCH -c 4
#SBATCH --mem=17G
#SBATCH -t 01:00:00
#SBATCH -J predict_verified_functional_profiles
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
HERE="$PIPES/predict_verified_functional_profiles"
W="$FASTALAKE_OUT_DIR/fix_eggnog_native"; mkdir -p "$W/scripts"
cp "$HERE"/*.py "$W/scripts"/
cd "$W/scripts"
"$PY" -u corrected_reproject.py PREDICT_noclust replay
"$PY" -u corrected_reproject.py PREDICT_noclust verified50
