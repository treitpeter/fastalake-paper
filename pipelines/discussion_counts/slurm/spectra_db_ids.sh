#!/bin/bash
# Spectra, database size, candidate competition and identifications per acquisition (spectra_db_ids/).
# Original jobs: 6073318 tasks.py, 6073323-6073326 pilots, 6073330 worker array (100 tasks), 6073437 collect.py.
set -euo pipefail
bash "${PIPELINE_DIR}/slurm/stage.sh"
cd "${FASTALAKE_OUT_DIR}/spectra_db_ids"
PY=${FASTALAKE_PYTHON:-python3}
SR="srun -p p.hpcl92,p.hpcl8,p.hpcl93 --nice=10000 -c2"
$SR -J sdb_tasks --mem=2G -t 00:20:00 $PY -B tasks.py
mkdir -p logs
jid=$(sbatch --parsable "${PIPELINE_DIR}/slurm/spectra_db_ids_worker.sbatch")
srun -p p.hpcl92,p.hpcl8,p.hpcl93 --nice=10000 -c2 --dependency=afterok:${jid} -J sdb_collect --mem=2G -t 00:15:00 $PY -B collect.py
