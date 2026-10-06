#!/bin/bash
# A2 preflight on a compute node: control, build, positive checks, negative controls (must exit 1).
A=$(cd "$(dirname "$0")/.." && pwd); PY=/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python
cd $A
echo "== prepare check331"; $PY run/code/prepare_p15.py check331; echo "exit=$?"
[ -f run/tasks.tsv ] || { echo "== prepare build"; $PY run/code/prepare_p15.py build; echo "exit=$?"; }
for m in searches grouping331 grouping347; do
  echo "== $m"; $PY scripts/preflight.py $m; echo "exit=$?"
  echo "== $m --negative (must be 1)"; $PY scripts/preflight.py $m --negative; echo "exit=$?"
done
$PY -c "import openpyxl, pandas; print('openpyxl', openpyxl.__version__, 'pandas', pandas.__version__)"
