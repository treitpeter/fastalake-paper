# Open items: abstract_discussion

- UNTRACED: "At a matched 1% false discovery rate" (v28d line 35) is the pooled peptide q ≤ 0.01 acceptance applied upstream in the cross-fit rescoring (Methods line 775); the shipped table carries only accepted counts, so the 1% cannot be recomputed here (NUMBERS.tsv row 2).
- TEXT-GAP: the manuscript names the seed (20261004) and draw count but not the generator (NumPy default_rng, PCG64), the interval type (percentile, 2.5th and 97.5th percentiles), or the stratum order within the first stream (PREDICT, CapScan all, CapScan stool, intestinal fluid); the exact stool upper bound -2.26 depends on that order (Codex N0206 drew -2.47 with another order) (STATS.tsv rows 1-2).
- TEXT-GAP: four PREDICT acquisitions without a patient identifier are treated as their own bootstrap clusters (97 clusters = 93 patients + 4); the manuscript does not state this (STATS.tsv row 1).
- TEXT-GAP: the software of the bootstrap (Python 3.12.4, NumPy 2.2.6, pandas 2.3.3) is not stated; Methods line 777 names Python 3.10.9 and NumPy 1.26.4 only for the pooled rescoring (STATS.tsv all rows).
- TEXT-GAP: the Abstract "21–36% more in intestinal fluid" is the range of the four regional medians (colon 21.1%, jejunum 35.8%); a reader can link it only through Results line 123. The pooled intestinal-fluid median is 28.3% (Results line 133). Suggested v28f wording if the authors want it explicit: "21–36% more across four intestinal sampling regions" (AUTHOR decision; no value changes) (NUMBERS.tsv rows 6-7).
- TEXT-GAP: "median database-selection times of 1–3 min per acquisition" (line 193) is computed here as per-acquisition sum of the lookup-extraction (05_persample) and parsimony (06_parsimony) task times, median per cohort (1.23-3.21 min); the sum of the two stage medians gives 1.12-3.14 min. Both round to 1–3; the definition is not stated (NUMBERS.tsv rows 16-17).
- ZENODO: the Source Data package has no sheet for the Abstract cross-fitted gains (map rows ABS-02, ABS-02b, ABS-03, ABS-03b still anchor the superseded 9-14% and 47-64%; D7-l, D7-m, D7-n, D7-o, DIS-01, DIS-02, DIS-48 anchor deleted sentences) and no map row for the Discussion timing sentence; Supplementary_Figure_15/compute_stage_summary.tsv in the package still has 331 CapScan acquisitions, and Review_v16/Database_sizes/database_sizes_710.tsv still has 710 FASTAs. Request: add map rows pointing to source_data/Fig4/crossfit_per_acquisition.tsv (FL_pooled_9, MG_pooled_9), source_data/Fig1/acquisition_sizes_726.tsv and source_data/SFig15/compute_*.tsv.
- PIPELINE: capscan_crossfit_10fold: the 10-fold cross-fit (design/worker.py, sage_lda.py, sage_q.py, run/analyze_final.py) is not in pipelines/; shared request to the orchestrator and to the Fig4 unit to agree one id for it (this unit uses capscan_crossfit_10fold).
- AUTHOR: none of the Abstract or Discussion numbers here carries a reference to Source Data; the Abstract gains rely on Results lines 123 and 133 and Methods line 777.

## shared requests

- Fig4 unit: please keep `source_data/Fig4/crossfit_per_acquisition.tsv` with columns `cohort`, `sampling_site`, `participant_sort_rank`, `FL_pooled_9`, `MG_pooled_9` unchanged (sha256 5bcdbe88...), or tell this unit if it changes; this unit reads it.
- Orchestrator: assign a common upstream id for the 10-fold cross-fit (proposed `capscan_crossfit_10fold`).

## v30 reference (added 2026-10-06)

- UNTRACED: the v28d quote of NUMBERS.tsv rows 3 (1%), 4 (3.5), 5 (5.5%), 17 (1), 18 (3), 19 (1) was not located in the v30 text (`visible_v30f.txt`); the wording or the value changed or the sentence was removed. v30_status NOT_IN_V30; the v28d columns are unchanged.
