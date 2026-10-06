#!/bin/bash
# Tasks A-E were run as interactive srun steps (logs/<task>.out, .err). The commands below
# reproduce them with the job names, cores and memory requests recorded by sacct
# (job ids in RECEIPTS.md). Each task reads only existing outputs.
set -euo pipefail
bash "${PIPELINE_DIR}/slurm/stage.sh"
cd "${FASTALAKE_OUT_DIR}"
PY=${FASTALAKE_PYTHON:-python3}
SR="srun -p p.hpcl92,p.hpcl8,p.hpcl93 --nice=10000"
$SR -J dcA  -c2  --mem=8G  -t 01:00:00 $PY -B a_accumulation.py      > logs/A.out  2> logs/A.err
$SR -J dcBp -c6  --mem=12G -t 00:30:00 $PY -B b_gmsc.py pilot        > logs/B_pilot.out 2> logs/B_pilot.err
$SR -J dcB  -c32 --mem=24G -t 03:00:00 $PY -B b_gmsc.py              > logs/B.out  2> logs/B.err
$SR -J dcBx -c2  --mem=1G  -t 00:05:00 $PY -B b_extra.py             > B_extra.txt
$SR -J dcCp -c32 --mem=16G -t 00:20:00 $PY -B c_homologs.py pilot    > logs/C_pilot.out 2> logs/C_pilot.err
$SR -J dcC  -c64 --mem=32G -t 04:00:00 $PY -B c_homologs.py          > logs/C.out  2> logs/C.err
$SR -J dcD  -c2  --mem=2G  -t 00:30:00 $PY -B d_clustering.py        > logs/D.out  2> logs/D.err
$SR -J dcE1 -c32 --mem=32G -t 02:00:00 $PY -B e_candidates.py e1     > logs/E1.out 2> logs/E1.err
$SR -J dcE2 -c16 --mem=48G -t 02:00:00 $PY -B e_candidates.py e2     > logs/E2.out 2> logs/E2.err
