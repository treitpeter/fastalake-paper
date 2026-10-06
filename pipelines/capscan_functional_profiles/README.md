# capscan_functional_profiles

Internal ID: A81 (previous folder name A81_capscan_functional_profiles_347)

What are the taxonomic and functional profiles, regional human share and own-stool similarity of the 347-acquisition CapScan set?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `capscan_fig5bcd.py` | call chain of the used outputs |
| `profiles.py` | call chain of the used outputs |

Launchers: `slurm/stage.sh`.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: taxonomy_signal.tsv and OG_other_unresolved.npz (grouping, evidence, annotation repair, consensus, profiles) for Fig5 b-d | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/Fig5, tables/SData1, text_statistics/supplementary_methods_8_11. What is used: taxonomy_signal.tsv and OG_other_unresolved.npz (grouping, evidence, annotation repair, consensus, profiles) for Fig5 b-d; own_stool_participants.tsv (CAPSCAN_REGIONAL_PAIRED15, 15 participants, 75 cells); regional human-share and own-stool tests, own_stool_distance_matrices.tsv

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

No Slurm job ids are recorded for this step (see OPEN_ITEMS.md). It ran without a recorded accounting entry.

## Choices

Approach: grouping, evidence, annotation repair, consensus and profiles (taxonomy_signal, OG_other_unresolved) and the Fig 5 b-d tests.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| cells | 15 participants, 75 own-stool paired cells (CAPSCAN_REGIONAL_PAIRED15) | participants with all five regions | none |
| profiles | KO, OG, pathway, host | as in the production arm | none |

## Known limitations

See `OPEN_ITEMS.md`.
