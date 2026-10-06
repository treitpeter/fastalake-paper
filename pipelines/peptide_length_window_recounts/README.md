# peptide_length_window_recounts

Internal ID: len950 (v30 9-50-residue recounts, 6 Oct 2026)

Question. What do the headline peptide counts, ratios and recoveries become when they are restricted to peptides of 9-50 residues? Every 9-50 value is a filter of the accepted 5-50-residue list (rank 1, target, peptide q <= 0.01, modifications stripped, I/L merged) by canonical peptide length. q values are not re-estimated and no search is repeated (Methods: "We obtained 9-50-residue values by filtering the accepted 5-50-residue lists without re-estimating q values").

Each script also recomputes the 5-50 value as a control; a recount is used only when the 5-50 control reproduces the published table or render exactly.

## Code

One subfolder per figure unit whose values were recounted. Scripts are as they ran; only absolute paths, interpreter and font paths were replaced by environment variables.

| folder | script | what it computes | reads |
|---|---|---|---|
| Fig2 | new_values_fig2.py | CAMPI peptide counts, reference identity recovery, member-specific peptides, UpSet intersections at 9-50 | source_data/Fig2 through panels_len950.py |
| Fig2 | panels_len950.py, build_fig2_len950.py | Fig. 2 panels with MIN_LENGTH switch; control (5-50) and len950 renders | source_data/Fig2 |
| SFig03 | new_values_sfig03.py, panels_len950.py, build_sfig03_len950.py | CAMPI member support and second-engine recovery at 9-50 | source_data/SFig03 |
| Fig3 | new_values.py, panels_len950.py, build_fig3_len950.py, q3.py | healthy-stool transfer and pooling medians at 9-50 | source_data/Fig3, SFig07 recount tables |
| Fig4 | capscan_regions.py, same_vs_crossfit.py | CapScan region gains and same-spectrum vs cross-fitted counts at 9-50 | source_data/SFig06, source_data/Fig4 |
| SFig06 | build.py, panels_len950.py | same-spectrum benchmark panels at 9-50 (roster 726) and 5-50 control (roster 710) | source_data/SFig06 |
| SFig07 | hstool_recount.py | healthy-stool common-FASTA identities and complete features at 9-50 | A21 and joint-search accepted identities and feature matrices |
| SFig07 | clin_recount.py | clinical own-FASTA vs common-FASTA accepted peptides at 5-50 and 9-50 | common-FASTA comparison manifest and Sage PSM tables |
| SFig07 | panels_len950.py, build_sfig07_len950.py | Supplementary Fig. 7 at 9-50 and its value tables | recount tables above |
| SFig08 | recount_curation.py, build.py | PREDICT curation-arm accepted peptides at 9-50 | A16 manifest, Sage PSM tables, own-metagenome peptide lists |
| SFig09 | build.py, build_ylim.py | CapScan acquisition counts, per-specimen values and tertiles at 9-50 | source_data/SFig09 |

Launchers (slurm/): stage.sh (restores the original run layout under FASTALAKE_OUT_DIR, which the scripts rely on), clin_recount.sbatch, hstool_recount.sh, curation_recount.sh, table_recounts.sh.
Environment variables: FASTALAKE_DATA_DIR (data root laid out as the lab tree), FASTALAKE_OUT_DIR (run directory), FASTALAKE_PAPER_REPO (this repository, for source_data/ and figures/), PIPELINE_DIR (this folder), FASTALAKE_PYTHON, FASTALAKE_FONT_DIR (Söhne fonts for renders), FASTALAKE_PY_MPL370 (Fig. 2 legacy-matplotlib render step only).

## Inputs

Repository Source Data (source_data/Fig2, Fig3, Fig4, SFig03, SFig06, SFig07, SFig08, SFig09) and, for the three recounts that start from search outputs, the manifests and tables below (sha256 in RECEIPTS.md):

- clin_recount.py: MM_review_analysis/end_of_day_2026-09-15/common_cache/comparison/inputs/manifest.json and the Sage results it names.
- hstool_recount.py: MM_review_analysis/A21_7_hstool_quant_comparison/{inputs/manifest_revision2.json, results/{original,union_individual,joint_corrected}/{accepted_identities.tsv.gz, same_file_id.features.npz}} and the joint_fileq outputs of the corrected joint search (job 6048637).
- recount_curation.py: MM_review_analysis/recovery_2026-09-15/public_comparison_recovery/A16/inputs/manifest.json, the PREDICT search-input inventory (SEARCH_INPUTS.json), own-metagenome peptide lists (R4K_ownmetag) and the five replayed searches (supplement6_repair).

## Outputs

Per unit: NEW_VALUES.tsv (each recounted value with its 5-50 control and source), plus recount tables (Fig2/upset_intersections_FL_9_50.tsv, SFig07/work/hstool/*_len.tsv, SFig07/work/clin/part_*.tsv, SFig08/curation_counts_5_9.tsv, SFig09/capscan_acquisitions_9_50.tsv, SFig06 plotted values, stats JSON). Hashes in RECEIPTS.md. The rebuilt figures themselves belong to the figure units.

## Resources and run time

- clin_recount: 10-task array (job 6072755), 4 cores allocated, 2 GB request, peak RSS 0.13-0.18 GB per task, 1.2-1.3 min per task.
- hstool_recount: srun job 6072769, 4 cores, peak RSS 1.66 GB, 35 s.
- All other scripts read only Source Data and finish in seconds to about a minute on 2 cores; they were run with srun and their job ids were not recorded.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| length window | 9-50 residues of the canonical (modification-stripped, I/L-merged) peptide | only 9-50-residue predictions were eligible for lookup | 5-50 values kept in Source Data and in every NEW_VALUES.tsv as control |
| q values | not re-estimated; 9-50 is a subset of the accepted 5-50 list | the searches estimated confidence over 5-50; refitting would change the comparison | none |
| acceptance | rank 1, label 1 (target), peptide_q <= 0.01 | same as the producers of the 5-50 tables | none |
| control | each 5-50 value recomputed and required to equal the published value or render | detects any change other than the length filter | none |
| SFig06 roster | 726 acquisitions at 9-50, 710 for the 5-50 control | roster of the v28f figure vs the earlier table | LEN950_ROSTER |

## Notes

Three layout probes used only to place relabelled text in the Fig. 3 composite (probe.py, probe2.py, probe3.py) are not included; they compute no value.

The code is the original code as it ran. Only absolute paths, interpreter paths and log paths were replaced by environment variables. Nothing was rerun for this folder.
