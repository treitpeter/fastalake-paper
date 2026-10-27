#!/bin/bash
# usage: run_units.sh <repo dir> ; executes the CapScan-dependent units' notebooks (written back in that copy)
set -uo pipefail
cd "$1"
export MPLBACKEND=Agg
for u in figures/Fig4 figures/SFig06 figures/SFig09 text_statistics/supplementary_methods_6_7 text_statistics/methods text_statistics/abstract_discussion; do
  ${FASTALAKE_PYTHON:-python3} tools/run_notebooks.py --timeout 3600 $u > "$1/RUN_$(basename $u).log" 2>&1; echo "$u rc=$?"
done
