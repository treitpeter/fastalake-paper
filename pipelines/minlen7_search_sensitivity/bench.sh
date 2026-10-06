#!/bin/bash
# Memory benchmark = real pilot searches (frozen run_task.py --pilot, writes run/pilot/) under /usr/bin/time -v.
A=$(cd "$(dirname "$0")/.." && pwd); PY=/mpcdf/soft/SLE_15/packages/x86_64/anaconda/3/2023.03/bin/python
export OMP_NUM_THREADS=16
for i in "$@"; do
  /usr/bin/time -v $PY $A/run/code/run_task.py $i --pilot > $A/logs/bench_$i.out 2> $A/logs/bench_$i.time; echo "idx $i exit=$?"
  grep -E "Maximum resident|Elapsed" $A/logs/bench_$i.time
done
