# predict_crossomics_refit

Internal ID: A84 (previous folder name A84_predict_crossomics_refit)

Do the functional profiles of host protein, non-host protein, DNA (metaG), RNA (metaT) and metabolites of PREDICT (77 patients) and CapScan (12 participants, 60 participant-region cells) agree in sample geometry, after the seed-orthologue correction of the non-host profiles? Used by Fig4 i, text line 139, SFig11 and supplementary methods 1-5 and 8-11. `check_prevalence.py` supplies the transform and Mantel-type evaluate functions that `crossomics.py` imports by parsing.

Status: code copied as it ran (path edits only), not re-run in this repository. See RECEIPTS.md for original job ids and hashes.

## Inputs

| input (LAB_ROOT-relative) | sha256 |
|---|---|
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/CAPSCAN_noclust/host_representative_intensity.tsv` | `4f0ee2b452a1b957bcf515c47d9602d2917aa6a7a911cecd725949fcba8539a7` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/CAPSCAN_noclust/microbial_KO.tsv` | `07b03f2a00328458fae14893409274a77eb30a1e69b736acd4f33f58e9df1eaf` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/CAPSCAN_noclust/microbial_OG.tsv` | `14d0eb1f7362e7650c41501385235b71b9c0ed21bc9816e49e2e38fb678dee29` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/CAPSCAN_noclust/microbial_Pathway.tsv` | `972e638c9bba8f4a2e6ddacc625e1011c950ed4e0d5d08150f81dec45c0ebdd6` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/PREDICT_noclust/host_representative_intensity.tsv` | `bf78627a6b6657263a4d52210f39bd966a90de096616ddf9a8f2f94e543fd1cb` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/PREDICT_noclust/microbial_KO.tsv` | `474088cfdbb5688b45e10fefa6007fab55a58b0d90b83dbd79fdcd0db108d672` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/PREDICT_noclust/microbial_OG.tsv` | `61df97764f988071c3bc2c1ffd5797f51f1e1460b066bf60a19603731b303b22` |
| `MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/crossomics_repaired__simplified/source_data/PREDICT_noclust/microbial_Pathway.tsv` | `51bfe74644fcb47ae5029f8a544748adca5789588b69401b84b346dc416fb61b` |
| `analyses/capscan_metag_analysis_2026-09-10/aggregate/metaG_KO.tsv` | `17b81e780ba40fb7d90993650f2cd5fa3194a96f945426fc3e59296d2414feab` |
| `analyses/capscan_metag_analysis_2026-09-10/aggregate/metaG_OG.tsv` | `82bc95cf8f010a8aca6caa93a5b0963fcffa8e5ffd5cb404bea439f46748279c` |
| `analyses/capscan_metag_analysis_2026-09-10/aggregate/metaG_Pathway.tsv` | `8d13d06f5399078c848ca50cac1fec1a1ca31148666a6e6022e87ad0c742763e` |
| `analyses/capscan_metag_analysis_2026-09-10/manifests/runs.tsv` | `0edba6053720e2ac1ccb81f15283aed1bf927416ecf6b772cc90951f63a96306` |
| `analyses/enterotypes_2026-09-07/CAPSCAN_metabolite_intensity.tsv` | `1ed653a2d738d6cfb8757a73e472d54f87514577caad29ebcfe038be45a6de13` |
| `analyses/enterotypes_2026-09-07/CAPSCAN_specimen_design.tsv` | `75bfffdb1007257acdf4f2c6b632356f24cf06074590231f6967786876d7ed15` |
| `analyses/enterotypes_2026-09-07/PREDICT_fullstack_primary_samples.tsv` | `5e967b0aa182f2dfd6b2d902d2d1b393dd96caf43ef619575ee5c5b972aa09a7` |
| `analyses/enterotypes_2026-09-07/PREDICT_metaG_KO.tsv` | `077b3f48920f8041e579acaf80d0ff285de4dc4188e8161593ebaf00ebcef034` |
| `analyses/enterotypes_2026-09-07/PREDICT_metaG_Pathway.tsv` | `cc692926acf327b0ccc9f9e71489d4a13a2138432e4d4cd6c1b6d4cab89d3d52` |
| `analyses/enterotypes_2026-09-07/PREDICT_metaT_KO.tsv` | `4cb6c0000a46daefeb7073944377d5779efc085b6b5d617851247e9fe7015a51` |
| `analyses/enterotypes_2026-09-07/PREDICT_metaT_Pathway.tsv` | `d8d267bc41a0610247e4c8fbb46e5c2fdf4dc131a31072cdfc306ccfc5fedb85` |
| `analyses/enterotypes_2026-09-07/PREDICT_metabolite_intensity.tsv` | `0663b71545e31e677368c1fcbfe7d2c503e987f4e5a4eff677da6d8848418cc9` |
| `analyses/enterotypes_2026-09-07/metaG_OG_TPM.tsv` | `c39df737891e118b92c3251aca87be433fae886fe16ceca2ad4e8eb55ebe1f73` |
| `analyses/enterotypes_2026-09-07/metaT_OG_TPM.tsv` | `351a6cbeead4b21b7b5aa40e4c27918170f8a6093b2fa76aab8dc3083f792266` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/full/functional_correlations.tsv` | `ea6dc8d014c6a12c396c4b6cb6beac641351e42b98654ae4b4214fcd0d1cd5bc` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/full/native_SOURCE_DATA_correlations.tsv` | `a02fb8ed34d2f9f68771b75aca4db137edd2a1965a138b57b528dcf6e901a693` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/functional_distance_matrices.npz` | `dd9eb6e595e474ae05d57614e6fbfe1567ba7382211019ffc12acff4d615b632` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/native/INPUT_MANIFEST.json` | `437c76019adcf128d7f2bdf5d7caaf97222b76c15bc1d0c1df8bc925d2df3fb8` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/native/native_pyramids.py` | `64ba668b266b9be42bb0a3f273878a2b93daba53022cc3c5a70ab04d4d56c8da` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/native_CAPSCAN_distance_matrices.npz` | `8cbb306cccd34b78639975fb549749c1f909b74aff31f6b8bcc8e51671ad5b5c` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/native_PREDICT_distance_matrices.npz` | `d579ad4338c51486f3a5d0fb3147da003c35cb9d56c19d70465057fefe43cf3e` |
| `review/v15_integration/integration_2026-09-29/source_data/current/Crossomics/provenance/check_prevalence.py` | `b05f65dac2aaf53f0f4d238df674e90d17d7632363e3c8d74bc5fcbabc352c10` |

## Outputs

`cross_baseline/` and `cross_refit/`: `functional_correlations.tsv`, `native_correlations.tsv`, `feature_coverage.tsv`, `signal_retention.tsv`, `functional_distances.npz` (77 x 77 PREDICT and CapScan distance matrices), `native_distances.npz`, `*_method_comparison.tsv`, rosters (hashes in RECEIPTS.md).

## Slurm resources and run time (original logs, `sacct`)

| step | job | cores | requested mem | peak RSS | wall |
|---|---|---|---|---|---|
| baseline (`crossomics.py baseline`) | 6044302 | 2 | 12G | 0.84 GB | 43 s |
| refit (`crossomics.py refit`, one of five steps) | 6048399 | 4 | 16G | 1.22 GB whole job | 4 min 9 s whole job |

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| permutations | 1999 | protocol | none |
| seed | 20260907 | fixed in advance | none |
| prevalence filter | 0 (primary, all observed) and 0.10 (sensitivity) | the 10% filter changes eligibility, so both are reported | native pyramids use 0.10 as primary (A86) |
| BH families | PREDICT 28 functional tests, CapScan 16 per prevalence; native 10 PREDICT, 6 CapScan | complete families, no selection | none |
| PREDICT unit | patient (77 earliest complete-stack stool specimens) | one specimen per patient | none |
| CapScan unit | participant-region cell, 12 DNA-complete participants x 5 regions | DNA availability | mixed 12/14-participant roster as sensitivity |
| CapScan permutation | whole participant blocks, region residualized | repeated measures | none |
| acquisitions to specimens | relative profiles, median within specimen, mean within participant-region | as the original producer | none |
| baseline replay | tolerance 1e-10 on rho, p, q | proves the refit changes only the non-host profiles | none |

## Notes

Code edits relative to the original: lab absolute paths replaced by `FASTALAKE_DATA_DIR`; lab-tree roots located by a marker file replaced by `FASTALAKE_DATA_DIR`. No statistics changed. The scripts keep their original integrity gates (`verify()`, `gate()`, hash checks against `*_CODE.json` and `*_INPUTS.json` manifests of the frozen run). Those manifests contain absolute lab paths and are not in the repository; the hashes they held are listed in RECEIPTS.md. Because the code here carries the path edits described below, its own sha256 differs from the original (both are in RECEIPTS.md), so a re-run needs the manifests regenerated for this copy.`common.py` of A70 provides the shared helpers and must be on the work path. Native DNA and RNA distances are taken from the hash-verified native source (A86). CapScan region labels are provisional pending the authoritative sample key (stated in the original run receipt). Group matrices (`groups/`) come from `build_groups.py`, outside this pipeline group. `crossomics.py` also executes `check_prevalence.py` by `functions(X/'provenance/check_prevalence.py', ...)`, so the file must sit at `provenance/check_prevalence.py` under the data directory layout; it has to be staged there before the run.
