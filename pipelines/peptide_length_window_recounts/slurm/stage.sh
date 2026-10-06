#!/bin/bash
# Stage the code into the run directory in the layout the original run used:
#   <FASTALAKE_OUT_DIR>/<Unit>/{*.py, work/*.py, build/*.py}
# Several scripts locate their inputs and outputs relative to their own file.
set -euo pipefail
: "${PIPELINE_DIR:?set PIPELINE_DIR to this pipeline directory}" "${FASTALAKE_OUT_DIR:?set FASTALAKE_OUT_DIR to the run directory}"
put() { mkdir -p "${FASTALAKE_OUT_DIR}/$2"; cp "${PIPELINE_DIR}/$1" "${FASTALAKE_OUT_DIR}/$2/"; }
for f in panels_len950.py build_fig2_len950.py new_values_fig2.py; do put Fig2/$f Fig2; done
for f in panels_len950.py build_sfig03_len950.py new_values_sfig03.py; do put SFig03/$f SFig03; done
for f in panels_len950.py build_fig3_len950.py new_values.py; do put Fig3/$f Fig3; done
put Fig3/q3.py Fig3/work
for f in panels_len950.py build_sfig07_len950.py; do put SFig07/$f SFig07; done
for f in hstool_recount.py clin_recount.py; do put SFig07/$f SFig07/work; done
for f in build.py panels_len950.py; do put SFig06/$f SFig06/build; done
for f in build.py build_ylim.py; do put SFig09/$f SFig09/build; done
for f in capscan_regions.py same_vs_crossfit.py; do put Fig4/$f Fig4/build; done
for f in build.py recount_curation.py; do put SFig08/$f SFig08/build; done
# The paper repository (source_data/, figures/) is found through FASTALAKE_PAPER_REPO when the
# scripts run outside the repository.
echo "staged into ${FASTALAKE_OUT_DIR}; export FASTALAKE_PAPER_REPO=<repository root>"
