#!/bin/bash
#SBATCH --cpus-per-task=2
#SBATCH --nodes=1
# usage: sbatch ... grouping.sbatch <driver.py>   (run331/grouping/group_fresh.py or run/grouping/group_fresh_347.py)
set -euo pipefail
PY=/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python
echo "driver $1 sha256 $(sha256sum $1 | cut -c1-64) on $(hostname) $(date -Is)"
/usr/bin/time -v $PY "$1" CAPSCAN_noclust
