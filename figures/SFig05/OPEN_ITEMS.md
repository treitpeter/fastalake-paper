# Open items: SFig05

- UNTRACED: The reservoir size 375,333,654 in the continuation legend (visible_v30f line 229) is not recomputed by this unit; the CAMPI fusion run records it as `reservoir_sequences` in its submission record (campi_fusion_375m, `status/SUBMITTED.json`), and the count belongs to the reservoir-size table of the database-size unit (NUMBERS.tsv row with printed value 375,333,654).
- MISMATCH: The CapScan fusion control in panels a-c pools 331 acquisitions (`fusion_acquisition_components.tsv`), not the 347-acquisition CapScan roster; either rerun the CapScan fusion control on the 347 roster or state in Methods (visible_v30f line 368) that this control used 331 CapScan acquisitions (panels a-c).
- AUTHOR: Line 30 of the v30 text reports only the donor-selected stratum (median 1.0% for 5-50 and 0.49% for 9-50 residues, 58 acquisitions) and cites Fig. 1e; if the self-selected stratum is reported alongside, the recomputed medians are 1.027% (5-50) and 0.510% (9-50), 5 acquisitions (STATS.tsv rows 3-4).
- TEXT-GAP: NUMBERS.tsv v30 columns give the v30 paragraph (`visible_v30f.txt`: 30, 45, 46, 136, 227, 229, 367, 373, 374, 376) and wording; one number is new in v30 (1.91%, paragraph 46; v28d columns NA).
- TEXT-GAP: The v30 legend (line 227) does not yet carry the v30 LEGEND_CHANGES.md edit for d ("Peptide-record donor assay: explicit-peptide self-selected ..."); the subpanel title it replaces was removed from the figure.
- AUTHOR: `fusion_acquisition_components.tsv` keys acquisitions by raw-file name (PREDICT `..._MicrobPred_<n>`, CapScan `..._T760_<n>`); confirm these run names are not tube identifiers before release.
- AUTHOR: The CAMPI cohort is labelled `CAMPI_MOCK` in `fusion_acquisition_components.tsv` and `CAMPI` in `fusion_pooled_current.tsv` and in panels a-c; the tables are left as produced.
- TEXT-GAP: Methods line 376 does not state the bootstrap construction details recovered from the code: percentile interval (2.5% and 97.5% quantiles), numpy default_rng seed 2026091051 shared across the eight estimates in fixed order, numpy 2.2.6 and scipy 1.16.3 (STATS.tsv rows 13-20).
- TEXT-GAP: The fusion control seeds (20260907, 20260908) and the FDRBench entrapment seeds (20260910, 20260911, 20260912) are not stated in the text (lines 45, 367, 374).
- TEXT-GAP: Line 376 says "approximately 30,000-protein CAMPI reference"; the reference FASTA has 29,673 records (`simulation_design.tsv`); the printed value matches at the nearest thousand.
- PIPELINE: healthy_stool_donor_entrapment: not in `pipelines/`; the donor and human-cohort fusion recount code and receipts are needed there (UPSTREAM.json).
- PIPELINE: campi_fusion_375m: not in `pipelines/`; the CAMPI fusion rerun code, Slurm launchers and receipts (jobs 5984954-5984958) are needed there.
- PIPELINE: restricted_reference_simulation: not in `pipelines/`; the simulation generate, evaluate and collect code and receipts are needed there.
- UNTRACED: The Slurm job id of the v30 restyle run (run_sf456.sh) is not recorded in its folder; output hashes and the pixel-identical regeneration are recorded instead (PROVENANCE.md).
