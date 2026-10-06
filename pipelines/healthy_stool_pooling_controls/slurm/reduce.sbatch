#!/bin/bash
#SBATCH --job-name=fig3_control_reduce
#SBATCH --partition=p.hpcl93
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time=01:00:00
#SBATCH --output=${FASTALAKE_DATA_DIR}/reduce_%j.out
#SBATCH --error=${FASTALAKE_DATA_DIR}/reduce_%j.err
set -euo pipefail
${FASTALAKE_DATA_DIR}/python ${FASTALAKE_DATA_DIR}/reduce.py
