#!/bin/bash
set -euo pipefail
export PYTHONNOUSERSITE=1
exec /mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python "$1/submit.py" gate "$2"
