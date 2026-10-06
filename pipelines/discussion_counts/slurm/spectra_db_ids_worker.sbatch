#!/bin/bash
#SBATCH --job-name=sdb_worker
#SBATCH --partition=p.hpcl92,p.hpcl8,p.hpcl93
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00
#SBATCH --array=0-99
#SBATCH --output=logs/%x_%A_%a.log
# pilot peak RSS 0.66 GB (PREDICT, 4 arms); 2G covers the largest MG FASTA (193 MB vs 81 MB in the pilot)
# Run from ${FASTALAKE_OUT_DIR}/spectra_db_ids (as slurm/spectra_db_ids.sh does).
/usr/bin/time -v ${FASTALAKE_PYTHON:-python3} -B "${FASTALAKE_OUT_DIR}/spectra_db_ids/worker.py" $SLURM_ARRAY_TASK_ID 100 "${FASTALAKE_OUT_DIR}/spectra_db_ids/parts"
