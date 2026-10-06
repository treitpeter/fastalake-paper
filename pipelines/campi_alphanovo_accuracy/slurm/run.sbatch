#!/bin/bash
#SBATCH -J an_acc
#SBATCH -p p.hpcl8,p.hpcl92,p.hpcl93
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1
#SBATCH --mem=640M
#SBATCH --time=00:30:00
#SBATCH -o logs/run_%j.out
#SBATCH -e logs/run_%j.err
# --mem: pilot S01 (55.7 MB predictions) peak RSS 276 MB; largest input S06 121.7 MB -> ~603 MB scaled, x1.05 = 634 MB
set -euo pipefail
/usr/bin/time -v ${FASTALAKE_PYTHON:-python3} -B "${PIPELINE_DIR}/accuracy.py"
