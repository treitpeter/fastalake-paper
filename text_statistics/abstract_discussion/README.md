# Abstract and Discussion statistics

The Abstract (v28d line 35) and the Discussion (v28d lines 181-193) print a few numbers that no figure panel carries: the cross-fitted peptide gains of FastaLake over matched metagenome databases, the count of selected databases smaller than the human reference proteome, the reservoir sizes and the computing times. This unit recomputes each of them from tables that other units ship in `source_data/`. It has no figure panel and no pixel check.

## Statistics

| statistic | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| cross-fitted gain | median over acquisitions of 100 x (FL / MG - 1), 9-50-residue peptides, databases selected from other spectra | 35 ("3.5–5.5% fewer ... in two stool cohorts and 21–36% more in intestinal fluid"); 183 | PREDICT 296 acquisitions (97 clusters); CapScan stool 56 (15 participants); four intestinal regions 51-65 each (13-14 participants) | none; 95% percentile bootstrap over participants | Statistic: cross-fitted gain over matched metagenome databases | source_data/Fig4/crossfit_per_acquisition.tsv | capscan_crossfit_10fold (not in pipelines/) |
| database sizes | selected FASTAs with fewer entries than the 20,431-entry human reference | 35, 181 ("All 726 databases ... smaller than the human reference proteome") | 726 acquisitions (FASTAs) | none (census) | Statistic: selected-database sizes against the human reference proteome | source_data/Fig1/acquisition_sizes_726.tsv; source_data/Fig1/human_reference.json | cohort_database_sizes |
| reservoir sizes | entries of the four study reservoirs | 193 ("233–375 million sequences") | 4 reservoirs | none | Statistic: reservoir sizes and computing time | source_data/Fig1/acquisition_sizes_726.tsv | cohort_database_sizes |
| computing time | lookup, selection and Sage wall times | 193 ("8–12 min per cohort", "1–3 min per acquisition", "less than 1 min") | Slurm tasks: lookup 1 per cohort; selection and Sage 15, 58, 306, 347 acquisitions | none | Statistic: reservoir sizes and computing time | source_data/SFig15/compute_stage_summary.tsv; source_data/SFig15/compute_task_records.tsv | common_resource_accounting |
| CAMPI members | members of the defined community | 35 ("eight-member community") | 8 species | none | Statistic: CAMPI community members | source_data/Fig2/member_summary.tsv | campi_reference_gap |

The tables live with the units that first needed them and are not duplicated (`source_data/README.md`); `source_data/abstract_discussion/README.md` lists them.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `MIN_LENGTH`, `MAX_LENGTH` | 9, 50 | the Abstract and Discussion report 9-50-residue peptides (Methods line 777) | 5-50 residues: Supplementary Fig. 6d,e and Methods line 777 |
| estimand | median over acquisitions of 100 x (FL / MG - 1) | as in the original analysis (`analyze_final.py`) | pooled-over-acquisition counts are not reported |
| stool strata | PREDICT (all matched acquisitions are stool) and CapScan `stool` | "two stool cohorts" | CapScan all 282: 18.6% (printed in no sentence of this unit) |
| intestinal fluid | range of the four regional medians (duodenum, jejunum, ileum, colon) | the Abstract prints a range; Results line 123 says "21–36% across the four intestinal sampling regions" | pooled intestinal fluid median 28.3% (Results line 133) |
| `B` | 10,000 | Methods line 777 | none |
| `SEED` | 20261004 | Methods line 777 | none |
| `CI_LEVEL` | 95, percentile interval | as in the original analysis | none |
| `MAIN_STRATA`, `REGION_STRATA` | two numpy PCG64 streams with this stratum order | reproduces the original intervals exactly; Methods line 777 states that regions use a separate stream | the stool interval differs between streams (-2.26 against -2.53 upper bound) |
| `SELECTION_STAGES` | `05_persample` + `06_parsimony`, summed per acquisition | "database selection" in the Discussion | sum of the two stage medians gives 1.12-3.14 min, also "1–3" |
| `LOOKUP_STAGE`, `SAGE_STAGE` | `02_evidence`, `07_sage` | "study-level lookup", "Sage searches" | none |

## Known limitations

- The cross-fitted counts come from the 578-acquisition matched roster (296 PREDICT, 282 CapScan). The 16 CapScan acquisitions of participant 15 have no matched metagenome and enter only as selection evidence.
- Three acquisitions used Sage's heuristic score after a singular pooled fit (`lda_fallback`); they are kept, as in the manuscript (Methods line 777).
- The 1% false discovery rate is a search parameter applied upstream; the shipped tables carry only accepted counts, so it is not recomputed here.
