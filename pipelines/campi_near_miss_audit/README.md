# campi_near_miss_audit

Internal ID: A03 (previous folder name A03_campi_near_misses)

Question. For discrepant de novo predictions on the 15 CAMPI acquisitions, what happens to the same spectrum in the lookup, family assignment, parsimony and final search? Supports the near-miss audit (84.8%, 34,019 of 40,097; 20.3%), Fig2 line 81.

Code in this folder: analysis.py, lib.py, peptide.py, refine_families.py, report_v2.py.
Launchers (slurm/): family_refinement.sbatch, pilot.sbatch, prepare.sbatch, report.sbatch, samples.sbatch, stage_frozen.sh.
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

inputs/manifest.json of the run (sha256 below) lists the 15 prediction CSVs, the archived FL/REF Sage results, the production extractor and the REF FASTA; each entry carries path, size, mtime and sha256, re-verified by lib.verify at run time.

## Outputs

results/acquisition_summary.tsv, acquisition_validation.tsv, pooled_summary.tsv, reference_annotation_aliases.tsv, REPORT.md, COMPLETE.json (hashes below). Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

prepare 2 CPU, 4 GB request, 16 s. pilot 4 CPU, peak RSS 0.58 GB, 2.3 min. samples 14 tasks, peak RSS 1.60 GB, longest 3.5 min. report 14 s. family refinement peak RSS 0.29 GB, 2.4 min. Requests: 4 to 12 GB.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| query population | score > 0, 9-50 residues, not low complexity, top 30% by score | Production lookup filter, lib.load_queries | none tested here |
| sequence identity | I and L collapsed | Sage reports I/L-collapsed peptides | none |
| comparator | archived REF search of the same spectra | REF contains surrogate strains, so agreement is operational | stated in README of run |
| fixed fraction | 0.30 | same fraction as production | not varied |

## Notes

Shared modules: lib.py and peptide.py live here and are used by A07, A14 and A44 (slurm/stage_frozen.sh copies them into the run directory as frozen/).

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
