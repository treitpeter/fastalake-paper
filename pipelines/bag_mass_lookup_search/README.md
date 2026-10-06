# bag_mass_lookup_search

Internal ID: A04 (previous folder name A04_bag_mass_lookup)

## Question

This pipeline produces: 48-acquisition bounded bag and mass-substitution searches (SEARCHES_COMPLETE.json, BOUNDS.json) for SFig02 panel e and lines 67, 561, 579; 336 arm-acquisition rows (bags_bounded48_*.tsv); original combined-route shuffled searches and candidate scans reused by STab6 paired nulls

Used by: figures/SFig02, tables/STab6, text_statistics/supplementary_methods_1_5.

## Code

The files are copies of the code that produced the used outputs, as it ran. The only change is mechanical: absolute lab paths were replaced by the token `LAB_ROOT` (Python and JSON) or `${FASTALAKE_DATA_DIR}` (shell and Slurm launchers), and personal identifiers were removed. Set `FASTALAKE_DATA_DIR` (inputs) and `FASTALAKE_OUT_DIR` (outputs) and substitute `LAB_ROOT` accordingly to run them. The original locations and sha256 of the unmodified files are in RECEIPTS.md.

| file | what it does (first docstring or comment line of the original) |
|---|---|
| `run_bounded.py` | Post hoc bounded-bag sweep on one finished acquisition of the native study. |
| `prepare.py` | Freeze the native AlphaNovo/FastaLake experiment; no submission in this step. |
| `submit.py` | Submit the frozen native AlphaNovo experiment with dependent full-cohort jobs. |
| `slurm/analysis.sbatch` | launcher |
| `analyze.py` | Recount completed searches and locate changed peptides in selected databases. |
| `campi_reference.py` | Corroborate S06 changes against its archived reference-database Sage search. |
| `check_methods.py` | Meaningful checks for sequence boundaries, context and classification. |
| `gain_loss_2026-09-24__prepare_submit.py` | Freeze every currently completed acquisition, then submit one attribution job. |
| `slurm/reference.sbatch` | launcher |
| `check_report.py` | Check report arithmetic, references, embedded PDF fonts and page boundaries. |
| `finalize.py` | no module docstring in original |
| `audit_results.py` | Independent table reconciliation; preserves unfinished/failed acquisitions. |
| `ions.py` | Explicit fragment definitions; no universal confidence threshold. |
| `fragment_validation__prepare_submit.py` | Freeze and run a CAMPI pilot followed by a dependency-gated five-case array. |
| `run.py` | Annotation replays and explicit raw-spectrum coverage for frozen comparisons. |
| `slurm/validate.sbatch` | launcher |
| `audit_and_plot.py` | no module docstring in original |

## Inputs

Inputs are named in the original code by `LAB_ROOT` paths (public sequence files, prediction files, search results of upstream pipelines). Their sha256 were not recorded by the original runs; see OPEN_ITEMS.md.

## Outputs

| output | sha256 |
|---|---|
| `SEARCHES_COMPLETE.json` | `14c3fe008397060e5574de25e9a2c406e15b7a6eaec226311b483a148888a846` |
| `SEARCHES_COMPLETE.json` | `26b400d9b47a7ab02ac55da5f27825c0404571936ca95cf84366351b1cf87ccc` |
| `SEARCHES_COMPLETE.json` | `ef37ca7874614965a4a91039e243bb822ac6a3f8ca0930ae3bf7380b612030d3` |
| `SEARCHES_COMPLETE.json` | `51ebfb0fa5b2b0c75fd6fb5ba8d2751a05014f5101251cadd1cdb3e802534800` |
| `SEARCHES_COMPLETE.json` | `931a3310fdd9cc79895290c34333a7a2b712abb3bd091dcb6c5567a54ea60dc6` |
| `SEARCHES_COMPLETE.json` | `412d9eeceae842a919a5ced5d229650270a9901edeecca07e8a466c1a1770d21` |
| `SEARCHES_COMPLETE.json` | `f7691330730a14ad5f8e6ad2ad8904aec931d077b6de6e9eef300a6bb00287fa` |
| `SEARCHES_COMPLETE.json` | `af5d975412324a771302bc1383c7dbfd4c269b588cccc5cc259eda69f020a8b1` |
| `SEARCHES_COMPLETE.json` | `423362c7a4d7571a9b9fe59d349cdafd46af8b29505039a43484125a5f313be9` |
| `SEARCHES_COMPLETE.json` | `751f5a168ffa723765391b08ff06de7f1d459da3f0e6784e097757eba1877f59` |
| `SEARCHES_COMPLETE.json` | `f4ed9b0a2d0572030200a858746156c1eb0e49dca8f6ccd9dc6ecc693949be0e` |
| `SEARCHES_COMPLETE.json` | `2d340c681e156fb92344ec73474e2910b4b3e8e57d1aca325c2ddadd87252e35` |
| `BOUNDS.json` | `bad6ca82de1dcb67bd995437726b40e27f9d266c29d3fb305d10b6f87d2688ca` |
| `BOUNDS.json` | `8930cef9174c226bd602d86a8e38c7eae7fe74554448bcd452a4d6ae93bc32fa` |
| `BOUNDS.json` | `413fe4e948c413f3d49a65e27c32beb68ce5687a88e986f432ac5622af73639f` |
| `BOUNDS.json` | `fddebfa999bbcba0d5d4fb18b4b8acd845d63450f312626ed8f9ac4ced09b56a` |
| `BOUNDS.json` | `0676ea151c6d2f97688f7faa323fde3e26e89a20e0bb9b017462be8941325e3b` |
| `BOUNDS.json` | `c582b8431ef5a0a9f8a919c7988ab73158e68478aa5782f442fca75605370e82` |
| `BOUNDS.json` | `8d94ae16193e0e579b264b74f70e2899c4d86ed85c615686d92fd680c1f607c8` |
| `BOUNDS.json` | `d8ef166b70c20a153ad0461bb8e90ccdafc20e4656479f7dc581d5d0234c5ff4` |
| `BOUNDS.json` | `715f356ec406d2f2da67399117574bf87204124d9d971512b176f402b7cff093` |
| `BOUNDS.json` | `c13e26c947a48a95fee2c6d8b23c5267b07d6bc074b69cd891c39ead37a39ced` |
| `BOUNDS.json` | `22af6b2b8a76632dc4ec5355d543ec701c46c0767f4f55fa39c66ec243395d2d` |
| `BOUNDS.json` | `911f17d166992e73e289fe97e814f282496d1bf60bb84a7a43635cb47e7a2905` |

## Slurm resources and run time

Read from the Slurm accounting database (sacct) for the job ids recorded in RECEIPTS.md.

| job | name | cores | requested memory | peak RSS (GB) | elapsed | state |
|---|---|---|---|---|---|---|
| 6000204 | FL_bag_gainloss | 4 | 16G | 0.4 | 00:01:16 | COMPLETED |
| 6000263 | FL_bag_ref_check | 4 | 4G | 0.0 | 00:00:11 | COMPLETED |

Requests in the launchers (`#SBATCH` lines):

- `slurm/analysis.sbatch`: no resource lines
- `slurm/reference.sbatch`: no resource lines
- `slurm/validate.sbatch`: no resource lines

## Choices

Values below are the top-level named constants of the code, with the comment written next to them in the original. A blank reason means the original gave none.

| parameter | value | why | alternatives |
|---|---|---|---|
| `STUDY` | `HERE.parent` | no comment in original code | not stated in code; see run_bounded.py |
| `REPO` | `STUDY / "fdr_revision/repo"` | no comment in original code | not stated in code; see run_bounded.py |
| `EXACT` | `1` | no comment in original code | not stated in code; see run_bounded.py |
| `MALL` | `1 << 10` | mass route, unbounded (reproduction only) | not stated in code; see run_bounded.py |
| `ARMS_B` | `{` | no comment in original code | not stated in code; see run_bounded.py |
| `ROOT` | `HERE.parents[2]` | no comment in original code | not stated in code; see prepare.py |
| `REPO` | `ROOT / 'worktrees/bags'` | no comment in original code | not stated in code; see prepare.py |
| `OLD` | `HERE.parent / 'full_reservoir_2026-09-23'` | no comment in original code | not stated in code; see prepare.py |
| `ROOT` | `HERE.parents[2]` | no comment in original code | not stated in code; see submit.py |
| `REPO` | `ROOT / 'worktrees/bags'` | no comment in original code | not stated in code; see submit.py |
| `ARMS` | `("exact", "bags", "mass", "combined")` | no comment in original code | not stated in code; see analyze.py |
| `CONTRASTS` | `(("exact", "bags"), ("exact", "mass"), ("exact", "combined"),` | no comment in original code | not stated in code; see analyze.py |
| `MOD` | `re.compile(r"\[[^\]]*\]")` | no comment in original code | not stated in code; see analyze.py |
| `SOURCES` | `{}` | no comment in original code | not stated in code; see analyze.py |
| `BASE` | `HERE.parent` | no comment in original code | not stated in code; see gain_loss_2026-09-24__prepare_submit.py |
| `PROTON` | `1.007276466621` | no comment in original code | not stated in code; see ions.py |
| `WATER` | `mass.calculate_mass(formula="H2O")` | no comment in original code | not stated in code; see ions.py |
| `ARMS` | `("exact", "bags", "mass", "combined")` | no comment in original code | not stated in code; see run.py |

## Adjacent-swap lookup on the full reservoir (added 2026-10-06)

Question. Does adding every single adjacent-residue swap of the retained queries of the representative healthy-stool acquisition (`20181231_QX7_PhGe_LC95_SA_WorkflowReplicates_2`) retrieve more authentic candidates from the 233,497,702-entry reservoir than it adds for a composition-preserving scrambled control (v30 Supplementary Methods paragraph 344)? Originally lab task R8 C2; the earlier A04 run of this pipeline applied the same swap definition to the 29,673-entry CAMPI reference only.

Code in `adjacent_swap/`: `build_manifest.py` (freezes inputs with sha256), `preflight.py` (re-verifies them at job start), `build_queries.py` (originals and originals plus swap-only variants; swap definition of A04: for each i, q[:i]+q[i+1]+q[i]+q[i+2:], skipping equal neighbours, after I/L collapse), `summarize.py` (proteins by originals, additional proteins by swaps only, rescued queries, swap-only proteins in the parsimony cover). Launchers: `slurm/adjacent_swap_pilot.sbatch`, `slurm/adjacent_swap_full.sbatch` (production fasta_extractor and parsimony_engine, `--min-length 9 --max-length 50 --normalize-il`, razor with hash-acc). Only paths were changed (FASTALAKE_OUT_DIR, FASTALAKE_DATA_DIR, PIPELINE_DIR, FASTALAKE_PYTHON); docstrings were added. Not rerun.

Inputs: 24,560 authentic retained queries (A02 REAL_queries.tsv) and 36,589 scrambled queries (seed 20260913, A02 SIM_20260913_queries.tsv); reservoir HSTOOL_lake.fasta (233,497,702 entries; full sha256 in inputs/reservoir.sha256.json of the run).

Result (Source Data `source_data/supplementary_methods_6_7/adjacent_swap_lookup.tsv`): authentic queries matched 6,671,447 proteins, and swap-only variants added 102,918 (b/a 0.0154); the scrambled control matched 11,401, and swap-only variants added 47,389 (b/a 4.16).

Resources (sacct, logs): pilot 5954786, 48 CPUs allocated, 10.9 min; full 5954787, 48 CPUs allocated (24 requested), 12.6 min, extractor peak RSS 66.2 GiB (mapped reservoir), parsimony 8.7 GiB.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| swap definition | every single adjacent transposition, equal neighbours skipped, after I/L collapse | A04 definition (frozen/benchmark.py line 42) | bag and mass lookup in this pipeline |
| acquisition | 20181231_QX7_PhGe_LC95_SA_WorkflowReplicates_2 | representative healthy-stool acquisition of Fig. 1 | none |
| control | composition-preserving scramble, seed 20260913 | one of the five A02 seeds, as A04 | scramble_end_to_end_workflow |
| lookup | fasta_extractor 9-50 residues, I/L normalised; parsimony razor, hash-acc | production settings | none |
