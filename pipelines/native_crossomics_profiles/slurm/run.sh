#!/bin/bash
#SBATCH -p p.hpcl92,p.hpcl8,p.hpcl93
#SBATCH --nice=10000
#SBATCH -c 1
#SBATCH --mem=7G
#SBATCH -t 00:30:00
#SBATCH -J native_crossomics_profiles
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
HERE="$PIPES/native_crossomics_profiles"
W="$FASTALAKE_OUT_DIR/native"; mkdir -p "$W"
cp "$HERE"/native_pyramids.py "$W"/
cd "$W"
"$PY" -u native_pyramids.py prepare
"$PY" -u native_pyramids.py pilot
"$PY" -u native_pyramids.py run --cohort PREDICT
"$PY" -u native_pyramids.py run --cohort CAPSCAN
