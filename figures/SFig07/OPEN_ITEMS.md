# Open items: SFig07

- AUTHOR: NUMBERS.tsv quotes and printed values are taken from the v30 text (`visible_v30f.txt`: Results line 57, Methods line 142, legend line 234, Supplementary Methods 7 line 395); the `v28d_line` column keeps the corresponding v28d line (115, 269, 437, 769). The v28d-only values (5-50 residues: 7.1%, 88.8%, 7.3%, 88.6%, 12,389, 236,835, 239,352, 239,013, 32.9%, 92.5%, 16,618, 211,589) and the untraced 5.9% PSM PEP sentence are no longer printed in v30 and were removed from NUMBERS.tsv (NUMBERS.tsv).
- AUTHOR: "Both retained about 94%" rounds 93.70% (common-FASTA individual) and 93.65% (joint, acquisition q) (NUMBERS.tsv rows 3-4).
- ZENODO: the Source Data package sheets for panels a and b are still the 5-50 tables; the 9-50 tables here have no package counterpart yet (PROVENANCE.md, Source data).
- TEXT-GAP: the manuscript does not state the software used for the medians (CPython 3.12.4, numpy 2.2.6, pandas 2.3.3 upstream and here; leading-species agreement CPython 3.10.9, numpy 1.26.4, pandas 2.3.3); Codex register rows S027-S032 (STATS.tsv).
- TEXT-GAP: the manuscript states no uncertainty interval for any number of this unit and no test; the code computes none (STATS.tsv, uncertainty column "none reported").
- TEXT-GAP: the legend does not say that panel b averages acquisitions within specimen and specimens within participant, as the upstream synthesis did (panel b).
- OPEN: the number of matched resolved functions behind the leading-species agreement (line 769) is not in the table supplied here; it is in A21 results/leaders/technical_function_agreement.tsv, not copied (STATS.tsv, n column).
- PIPELINE: healthy_stool_common_fasta, healthy_stool_joint_search_q, clinical_common_fasta and healthy_stool_clinical_length_recount_9_50 are not in pipelines/; their code is listed in UPSTREAM.json and was not ported, because it reads search outputs outside this repository (UPSTREAM.json).
- AUTHOR: panel b's CapScan participant values come from the 331-acquisition clinical synthesis, not the 347-acquisition CapScan roster; the participant count (15) is unchanged (panel b, PROVENANCE.md).
