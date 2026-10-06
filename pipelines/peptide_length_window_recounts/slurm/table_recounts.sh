#!/bin/bash
# Recounts that read only repository Source Data (seconds each), then the figure-unit
# value tables (NEW_VALUES.tsv) and length-window renders. Run after stage.sh, and after
# hstool_recount.sh and clin_recount.sbatch for SFig07 and Fig3.
set -euo pipefail
: "${FASTALAKE_OUT_DIR:?}" "${FASTALAKE_PAPER_REPO:?set to the repository root}"
PY=${FASTALAKE_PYTHON:-python3}
export PYTHONDONTWRITEBYTECODE=1 MPLBACKEND=Agg
O=${FASTALAKE_OUT_DIR}
RUN="srun -p p.hpcl92,p.hpcl8,p.hpcl93 --nice=10000 -c2 --mem=8G -t 30"
$RUN $PY -B "$O/Fig2/new_values_fig2.py"
$RUN $PY -B "$O/SFig03/new_values_sfig03.py"
$RUN $PY -B "$O/Fig3/new_values.py"
$RUN $PY -B "$O/Fig4/build/capscan_regions.py" "$O/Fig4/build/capscan_regions.json"
$RUN $PY -B "$O/Fig4/build/same_vs_crossfit.py" "$O/Fig4/build/same_vs_crossfit.json"
(cd "$O/SFig06/build" && LEN950_MIN=5 LEN950_ROSTER=710 $RUN $PY -B build.py && LEN950_MIN=9 LEN950_ROSTER=726 $RUN $PY -B build.py)
(cd "$O/SFig09/build" && $RUN $PY -B build.py && $RUN $PY -B build_ylim.py)
(cd "$O/SFig08/build" && $RUN $PY -B build.py)
for m in control len950; do $RUN $PY -B "$O/Fig2/build_fig2_len950.py" $m; done
for m in control len950; do $RUN $PY -B "$O/SFig03/build_sfig03_len950.py" $m; done
# Fig. 3 composite: build_fig3_len950.py <panel_pdf> <out_stem> (see its docstring); not needed for the values.
$RUN $PY -B "$O/SFig07/build_sfig07_len950.py" "$O/SFig07/build"
