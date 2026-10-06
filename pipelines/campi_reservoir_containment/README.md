# campi_reservoir_containment

Internal ID: A42 (previous folder name A42_campi_reservoir_containment)

Question. What share of each CAMPI member reference proteome is present, exactly, in the 375,333,654-entry reservoir? 91.3-100% (line 75); member proteome sizes (line 611).

Code in this folder: overlap.py, preflight.sh.
Launchers (slurm/): full.sbatch.
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), PIPELINE_DIR (this folder), PIPELINES_DIR (parent of all pipeline folders), FASTALAKE_PYTHON. Requires the original input tree under FASTALAKE_DATA_DIR (relative layout unchanged); the public bundle mapping is an open item.

## Inputs

inputs/manifest.json lists the reservoir FASTA (136 GB), joint headers TSV (81 GB), reference FASTA, members TSV and the extractor binary with sha256 recorded by the pilot job.

## Outputs

results/full/protein_presence_by_proteome.tsv, protein_presence_detail.tsv, peptide_tryptic_presence_by_proteome.tsv and others (hashes below). Hashes are in RECEIPTS.md.

## Resources and run time (from sacct)

32 CPU, 155 GB request, peak RSS 130.0 GB (production extractor memory-maps the reservoir), wall 165.8 min. Request chosen as measured 145.6 GB (job 5895250) x 1.05.

## Choices

| parameter | value | why | alternatives |
|---|---|---|---|
| match rule | exact protein sequence presence, plus tryptic-context peptide presence | containment is exact, not similarity | none |
| controls | shuffled proteomes must give zero presence | selftest in preflight | none |
| expected entries | 375,333,654 | reservoir size asserted at run time | none |

## Notes

The manuscript tables reservoir_protein_presence_by_proteome.tsv is a later copy of results/full/protein_presence_by_proteome.tsv (renamed); this is stated, not verified byte for byte.

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
