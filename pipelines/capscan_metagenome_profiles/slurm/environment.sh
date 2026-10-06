#!/bin/bash
# Source in the batch shell. Keep software versions explicit.
# Required: FASTALAKE_LAB_ROOT (data root), FASTALAKE_WORK_DIR (working directory).
# Python code addresses inputs as LAB_ROOT/<relative path>, so link it into the working directory.
cd "${FASTALAKE_WORK_DIR:?set FASTALAKE_WORK_DIR}"
ln -sfn "${FASTALAKE_LAB_ROOT:?set FASTALAKE_LAB_ROOT}" LAB_ROOT
module load bwa/0.7.17 samtools/1.20 fastp/0.23.2 salmon/1.10.0
export PATH=${FASTALAKE_EGGNOG_ENV}/bin:${FASTALAKE_LAB_ROOT}/../megahit_install/MEGAHIT-1.2.9-Linux-x86_64-static/bin:$PATH
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export PYTHONPATH=${FASTALAKE_LAB_ROOT}/analyses/capscan_metag_analysis_2026-09-10/vendor/mOTUs_3.0.3${PYTHONPATH:+:$PYTHONPATH}
