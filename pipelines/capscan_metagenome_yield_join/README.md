# capscan_metagenome_yield_join

Internal ID: none; previous folder name UNASSIGNED_capsule_metagenome_join

How are the CapScan DNA yield, assembly size and same-spectrum peptide counts joined per acquisition (282 rows)?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `build_figures_v15.py` | call chain of the used outputs |
| `sanitize_assembly.py` | call chain of the used outputs |

Launchers: none, the step ran interactively or inside another launcher.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: CapScan per-acquisition DNA amount, MG database size (sanitizer n_kept, CAPSCAN_sanitize_report.tsv) and same-spectrum counts (Source Data Figure4 capscan_assem | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/Fig4, figures/SFig09. What is used: CapScan per-acquisition DNA amount, MG database size (sanitizer n_kept, CAPSCAN_sanitize_report.tsv) and same-spectrum counts (Source Data Figure4 capscan_assembly_dna_yield); 282-row join of peptide counts, BioSample, region, participant

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

No Slurm job ids are recorded for this step (see OPEN_ITEMS.md). It ran without a recorded accounting entry.

## Choices

Approach: sanitizer n_kept and 282-row join of peptide counts, BioSample, region and participant.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| sanitizer | sanitize_assembly.py report (n_kept) | assembly sizes must count kept sequences | none |
| join | by BioSample and acquisition | no neighbour substitution | none |

## Known limitations

See `OPEN_ITEMS.md`.
