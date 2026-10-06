# RECEIPTS

## Original Slurm jobs

| role | job |
|---|---|
| CapScan native | 5990406 |
| PREDICT native | 5990407 |
| pilot (failed strict gate, preserved) | 5990373 |

## Code

| original (LAB_ROOT-relative) | file here | sha256 original | sha256 here |
|---|---|---|---|
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/native/native_pyramids.py` | `native_pyramids.py` | `64ba668b266b9be42bb0a3f273878a2b93daba53022cc3c5a70ab04d4d56c8da` | `d1c53e6a591849b18ed06955d6d29b9b235ecb6fb9acefae28660f8e47290804` |
| `review/v15_integration/integration_2026-09-29/cohort_reconciliation/crossomics/native.py` | `native.py` | `940a2ed12c51edc6f6ca166b7ab4b6cb8ab6c4e1426c7bd47c9987bc0565ac9c` | `fa7c7600cc8fa261eb7be97de47a442ee1b2ea392ea4c7d1da67abcbf2ee9281` |

## Inputs (sha256 recorded by the original run)

Input manifest hash given in README.

## Outputs (sha256 recorded by the original run)

| output (LAB_ROOT-relative) | sha256 |
|---|---|
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/PREDICT/RNA_native_cells.npz` | `4b3a32611d4c2c27394d90711a36f242964dfa4de9e6972e0d1eca0580fb6828` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/PREDICT/DNA_native_cells.npz` | `e505dfc5b18d025bd7b29aca76e3b5a96b4c93b58e3ab34bf20bf15372306170` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/PREDICT/gene_sequence_features.tsv.gz` | `5300769e0f8bc414f9b08fc12e99e8cbb70d19c45c7d598ec2085a48c1dfa643` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/PREDICT/feature_quality.json` | `b0023383d2c0b70c82b9514fe636d53b9fe858491cbec093dc3a92b27343d33a` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/PREDICT/protein_acquisition_roster.tsv` | `85575889daf0e3c1400b884f59cc04ff6f18210e95500a90b60761bfb2f33290` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/PREDICT/correlations.tsv` | `d26cf7892de78b62a82bd1cfa3f3c9b8d00a732e59f8c8dca565f4e6bb7121f2` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/PREDICT/distance_matrices.npz` | `d579ad4338c51486f3a5d0fb3147da003c35cb9d56c19d70465057fefe43cf3e` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/PREDICT/Host_protein_features.tsv.gz` | `d6638461ffb4ffdbf66cea9e0586c5899d8728125d3f48a1362eef19b7bcdba3` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/PREDICT/Non-host_protein_features.tsv.gz` | `8f62ea39062d05d42e8fcac93b3b4359670dade5c15cbe77b2b2531a565fda76` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/CAPSCAN/protein_metabolite_specimen_roster.tsv` | `9068f667d0595ad2190ae0c7738c17cb41e68133e3b1026f71a3b9851093ba53` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/CAPSCAN/DNA_native_cells.npz` | `5d2a3f515e6c1d6682b7254a330771d10729195d654d74f068f559af09640084` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/CAPSCAN/gene_sequence_features.tsv.gz` | `4090abd2b912f0ebb18c78c0fa0009eb421f9cf322602fb67585c831bf979d4c` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/CAPSCAN/feature_quality.json` | `cb84929c610d7b8e7915fe07c5cb51363666a03802e19bb3c9a1b48e498f64b5` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/CAPSCAN/correlations.tsv` | `933b3775ad6e5017d16665a7f1e049561f3c44e3f1ef83b9a16650be4a260fd9` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/CAPSCAN/distance_matrices.npz` | `8cbb306cccd34b78639975fb549749c1f909b74aff31f6b8bcc8e51671ad5b5c` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/CAPSCAN/Host_protein_features.tsv.gz` | `6e36b3781be7da81ed9e5735ae7493202b8c9664db9f25e1d8e838bd2ec7d1d8` |
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/results/CAPSCAN/Non-host_protein_features.tsv.gz` | `f58eb15c5db7ccedb2adc549ed71521cdb577f5f6f86fa6d5ce2ea8fbd77791f` |

| output (LAB_ROOT-relative) | sha256 |
|---|---|
| `review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics/native/SOURCE_DATA_retained_signal_by_profile.tsv` | `761a6b1242d432ae9f5b6fb1cf30221b5cd9b45a8cc492996aeaac57b8125524` |

## Path aliases

Lab path segments are shortened in code and tables (`jess_integration_2026-09-29` to `integration_2026-09-29`, `v10_mm_philip_comment_audit_2026-09-22` to `v10_comment_audit_2026-09-22`, `coauthor_comment_check_v25_2026-10-05` to `comment_check_v25_2026-10-05`). `FASTALAKE_DATA_DIR` must hold the frozen inputs under the shortened names; the original names are the lab directory names.
