# campi_second_engine

Internal ID: A07 (previous folder name A07_campi_second_engine)

Question. Does the FASTA interface carry the result into a second search engine? 45 MSFragger 4.1 searches with Percolator 3.06.5 (15 CAMPI acquisitions x FL, REF, MG). SFig03 panel d.

Code in this folder: closed_fragger.params, collect.py, search.py.
Launchers (slurm/): aggregate.sbatch, mg_pilot_v5.sbatch, pilot_v4.sbatch, searches_v5.sbatch, stage_frozen.sh.
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

inputs/manifest.json of the run (target FASTAs per acquisition and arm, MSFragger jar, Percolator, mzML) with sha256 (below).

## Outputs

results/between_engine_identities.tsv and AGGREGATE_VALIDATION.json (hashes below). Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

searches v5: 42 tasks, 8 cores (16 hardware threads allocated), 32 GB request, peak RSS 14.05 GB, longest task 2.9 min, 0.9 task-hours in total. MG pilot 0.5 min. Earlier attempts (5920933, 5920934, 5920954, 5920955, 5920993, 5920994, 5921017) failed or were cancelled; they are listed in RECEIPTS.md for completeness.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| engine | MSFragger 4.1 closed search | independent second engine | not applicable |
| validation | Percolator 3.06.5 target-decoy, --seed 913, accept peptide q <= 0.01 | Percolator accepts seeds 1-20000 | none |
| decoys | one fully reversed sequence per target, rev_ prefix | shared target/decoy sequences are counted, not assumed absent | none |
| digestion/ranges | copied from the Sage search (5-50 residues, 500-5000 Da, fixed C, variable M) | matched chemistry | closed_fragger.params |
| N-terminal acetylation and adaptive tolerance | disabled | match Sage | none |

## Notes

Java and Percolator binaries are taken from PATH (java) and from the run manifest. MSFragger jar is not redistributed (open item).

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
