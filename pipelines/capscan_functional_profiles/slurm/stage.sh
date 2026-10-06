#!/bin/bash
# C2 Fig 5 annotation-repair lineage, one stage per call. usage: stage.sh <stage>
# Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_WORK_DIR (working directory).
# Python code addresses inputs as LAB_ROOT/<relative path>, so link it into the working directory.
cd "${FASTALAKE_WORK_DIR:?set FASTALAKE_WORK_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg
OUT=${FASTALAKE_LAB_ROOT}/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/coauthor_comment_check_v25_2026-10-05/compute/C2_participant15; PY=${FASTALAKE_PYTHON:-python}; S=$OUT/scripts
ME=$OUT/fig5/ra/evidence_all_cohorts; MA=$OUT/fig5/ra/annotation_repair; G=$OUT/fig5/grouping/cohorts/CAPSCAN_noclust
echo "STAGE $1 HOST $(hostname) JOB ${SLURM_JOB_ID:-none} TASK ${SLURM_ARRAY_TASK_ID:-none} $(date -Is)"
case "$1" in
 evidence_pilot) test -f $G/COMPLETE.json; $PY $S/prepare_stage.py evidence $G; cd $ME; $PY -u run_projection.py 2 pilot; $PY $S/gate_capscan.py $ME 347 ;;
 evidence_full) cd $ME; $PY -u run_projection.py 2 full ;;
 inventory) $PY $S/prepare_stage.py inventory; $PY $S/prepare_stage.py extraction; $PY $S/prepare_stage.py supplement_prep; $PY $S/prepare_stage.py supplement_run; cd $MA; $PY -u inventory.py 2 ;;
 extract) cd $MA; $PY -u extract_native.py "$SLURM_ARRAY_TASK_ID" ;;
 supplement_prep) cd $MA; $PY -u prepare_supplement.py ;;
 supplement_pilot) cd $MA; $PY -u run_supplement.py pilot ;;
 supplement_full) cd $MA; L=$(printf 'chunk_%03d' "$SLURM_ARRAY_TASK_ID"); if [ ! -f supplemental_functions/$L.faa ]; then echo "no $L, nothing to do"; exit 0; fi; $PY -u run_supplement.py full "$SLURM_ARRAY_TASK_ID" ;;
 supplement_collect) cd $MA; $PY -u run_supplement.py collect ;;
 adapter) $PY $S/prepare_stage.py adapter; cd $MA; $PY -u build_adapter.py ;;
 reproject_pilot) $PY $S/prepare_stage.py reprojection; cd $MA/evidence_repaired; $PY -u repair_projection.py 2 pilot; $PY $S/gate_capscan.py $MA/evidence_repaired 347 ;;
 reproject_full) cd $MA/evidence_repaired; $PY -u repair_projection.py 2 full ;;
 consensus_pilot) $PY $S/prepare_stage.py consensus; cd $MA/annotated_origin_consensus_2026-10-01; $PY -u project.py 2 pilot; $PY $S/gate_capscan.py $MA/annotated_origin_consensus_2026-10-01 347 ;;
 consensus_full) cd $MA/annotated_origin_consensus_2026-10-01; $PY -u project.py 2 full ;;
 profiles) $PY -u $OUT/fig5/mc/run_profiles_capscan.py ;;
 fig5bcd) $PY -u $S/capscan_fig5bcd.py $OUT/fig5/mc/profiles/CAPSCAN_noclust $OUT/roster/design_347.tsv $OUT/fig5/fig5bcd_347 ;;
 *) echo "unknown stage $1"; exit 2 ;;
esac
echo "STAGE $1 DONE $(date -Is)"
