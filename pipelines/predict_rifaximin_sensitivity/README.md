# predict_rifaximin_sensitivity

Internal ID: A92 (previous folder name A92_predict_rifaximin_sensitivity)

Does the antibiotic on-off result of A91 change if rifaximin given for hepatic encephalopathy is handled as exposure, patients with it are excluded, or only its off rows are set unknown? Four exposure definitions (A current rule, B rifaximin for HE as exposure, C HE patients excluded, D HE off rows excluded). Used by STab3 and supplementary methods 8-11; outputs hold counts and statistics only.

Status: code copied as it ran (path edits only), not re-run in this repository. See RECEIPTS.md for original job ids and hashes.

## Inputs

Controlled: the PREDICT metadata table (sha256 `fc5dbf07df2fd25c7def3e9d583b2745c9f237b12a6a9c6e78ecefb249e90238`, asserted by the script) and `acquisition_accounting.tsv` of A91; public: `clinical_paired_statistics.tsv` (A91) and the EDS patient order (A70).

## Outputs

`out/{paired_endpoints_by_scenario,paired_composition_by_scenario}.tsv`, `counts.json`, `COMPLETE.json` (hashes in RECEIPTS.md). The workbook assembly of STab3 lies in the unit's own code, not here.

## Slurm resources and run time (original logs, `sacct`)

Job 6048159 (main): 2 cores, 4G requested, peak RSS 2.5 MB (sampled), wall 7 s; log of job 6048158 (8 s) also exists in the original log folder.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| antibiotic state rule | re-derived from the metadata flags, asserted equal to the state used by A91 | same rule as the published analysis | scenarios B to D change it |
| bootstrap and sign-flip | same sizes, seeds and BH family as A91 | scenarios comparable to the main result | none |
| output | counts and statistics, no patient or sample identifiers | controlled data | none |
| negative control | `--negative-control` mode exists in the script | gate check | none |

## Notes

Code edits relative to the original: lab absolute paths replaced by `FASTALAKE_DATA_DIR`; lab-tree roots located by a marker file replaced by `FASTALAKE_DATA_DIR`. No statistics changed. The scripts keep their original integrity gates (`verify()`, `gate()`, hash checks against `*_CODE.json` and `*_INPUTS.json` manifests of the frozen run). Those manifests contain absolute lab paths and are not in the repository; the hashes they held are listed in RECEIPTS.md. Because the code here carries the path edits described below, its own sha256 differs from the original (both are in RECEIPTS.md), so a re-run needs the manifests regenerated for this copy.Controlled inputs are read from `FASTALAKE_CONTROLLED_DIR`. The original `run.sbatch` (partition list, 1 core, 4G, 30 min, nice 10000; `main` mode runs the script then `independent_check.py`) is replaced by the launcher. `independent_check.py` of the original folder (a second implementation used as a check) is not copied because it did not produce the used outputs.
