#!/bin/bash
#SBATCH --job-name=FL_cf_pep
#SBATCH --partition=p.hpcl92
#SBATCH --nice=10000
#SBATCH --cpus-per-task=1
#SBATCH --mem=1500M
#SBATCH --time=00:30:00
#SBATCH --array=0-19
# Cross-fitted PEP sensitivity (as run: validate job 6057983, pilot 6057984, negative control 6057988,
# array 6057990, analysis 6058010). Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_OUT_DIR (run
# directory, holds PEP_MODEL.json and per_acq/), PIPELINE_DIR (this pipeline folder).
# MODE=validate runs "cf_pep.py validate 0" (8G); MODE=analyze runs cf_analyze.py (2G) after the array.
set -euo pipefail
cd "${FASTALAKE_OUT_DIR:?set FASTALAKE_OUT_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
[ -f PEP_MODEL.json ] || cp "${PIPELINE_DIR:?set PIPELINE_DIR}/PEP_MODEL.json" .
export PYTHONNOUSERSITE=1 PYTHONDONTWRITEBYTECODE=1
case "${MODE:-run}" in
  validate) ${FASTALAKE_PYTHON:-python} "${PIPELINE_DIR}/cf_pep.py" validate 0 ;;
  analyze)  ${FASTALAKE_PYTHON:-python} "${PIPELINE_DIR}/cf_analyze.py" ;;
  *)        ${FASTALAKE_PYTHON:-python} "${PIPELINE_DIR}/cf_pep.py" run "${SLURM_ARRAY_TASK_ID}" ;;
esac
