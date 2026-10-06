# Open items

- PIPELINE: histidine_spectrum_evidence: the case selection receipts it audits come from the corrected PREDICT spectrum selection wave (`predict_spectrum_review` source, worker sha256 `9b14e329ba3d68535596cb41ca367a0977a79dfa0036421d0357fd2f89b491a5` for the CapScan twin `audit_capscan_spectrum_cases_20260915.py`); the selection worker for PREDICT is not among the listed code paths.
- The script uses `MM_REVIEW_ROOT` internally; check that `slurm/run_task.sbatch` maps it to `FASTALAKE_DATA_DIR`.
