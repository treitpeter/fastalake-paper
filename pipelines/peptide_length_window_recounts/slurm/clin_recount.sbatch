#!/bin/bash
#SBATCH --job-name=len950_clin
#SBATCH --partition=p.hpcl92,p.hpcl8,p.hpcl93
#SBATCH --nice=10000
#SBATCH --cpus-per-task=2
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --array=0-9
#SBATCH --output=len950_clin_%A_%a.log
# Clinical own-FASTA vs common-FASTA recount at 5-50 and 9-50 (Supplementary Fig. 7 text values).
# Original run: job 6072755 (10 tasks). Writes FASTALAKE_OUT_DIR/SFig07/work/clin/part_<i>.tsv.
set -euo pipefail
: "${FASTALAKE_DATA_DIR:?}" "${FASTALAKE_OUT_DIR:?}"
mkdir -p "${FASTALAKE_OUT_DIR}/SFig07/work/clin"
/usr/bin/time -v ${FASTALAKE_PYTHON:-python3} -B "${FASTALAKE_OUT_DIR}/SFig07/work/clin_recount.py" "${SLURM_ARRAY_TASK_ID}" 10 "${FASTALAKE_OUT_DIR}/SFig07/work/clin"
