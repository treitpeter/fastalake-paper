# capscan_public_only_search

Internal ID: A08 (previous folder name A08_capscan_public_only)

Do the existing CapScan spectra identify as many peptides when searched against a database built only from public sequences (UHGG, GMGC, GMSC, HUMAN_SP; 233,497,702 sequences; no cohort metagenomes) as against the specimen metagenome?

Unit of analysis and outputs are listed below. Python code addresses inputs as `LAB_ROOT/<relative path>`: create a symlink `LAB_ROOT` in the working directory that points to the data root (the Slurm launchers do this from `FASTALAKE_LAB_ROOT` and `FASTALAKE_WORK_DIR`).

## Code

| file | role |
|---|---|
| `analysis.py` | call chain of the used outputs |
| `collect_a08.py` | call chain of the used outputs |
| `lib.py` | call chain of the used outputs |
| `peptide.py` | call chain of the used outputs |
| `search_helpers.py` | call chain of the used outputs |

Launchers: `slurm/evidence.sbatch`, `slurm/pilot.sbatch`, `slurm/samples.sbatch`.

## Inputs

| input | sha256 |
|---|---|
| upstream outputs (searches, rosters, annotations) used by the code above: CapScan public-sequence FastaLake searches (331) read by A16 for CapScan rows of SFig06 d,e | recorded in the RECEIPTS.md of the producing pipeline |

## Outputs

Used by: figures/SFig06, figures/SFig09. What is used: CapScan public-sequence FastaLake searches (331) read by A16 for CapScan rows of SFig06 d,e; per-acquisition FL_public and specimen-matched MG peptide counts (5-50 residues, q<=0.01) for 282 matched CapScan acquisitions, SFig09

Hashes of the outputs are in `RECEIPTS.md`.

## Slurm resources and run time (from sacct of the original jobs)

| jobs | cores | requested memory | peak RSS | longest task | summed task time | states |
|---|---|---|---|---|---|---|
| 681 tasks | 16, 2, 8 | 17G, 2G, 64G, 80G, 96G | 80.2 GB | 39 min 00 s | 37 h 20 min | COMPLETED 681 |

Partitions: p.hpcl8, p.hpcl92. Job ids: 5920984, 5920985, 5920986, 6065269, 6065270, 6065271, 6065329.

## Choices

Approach: public-only FastaLake database per acquisition: top-30% lookup of the pooled prediction evidence, 9-50 residues, I/L-equivalent, razor/hash-accession parsimony, then Sage with the frozen production settings including LFQ.

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| lookup depth | top 30% per acquisition | frozen production rule, so that only the sequence source changes | none run here |
| peptide length window for lookup | 9-50, I/L-equivalent | frozen production rule | 5-50 counted in the comparison (peptide_q01), 9-50 reported alongside |
| source classes | UHGG, GMGC, GMSC, HUMAN_SP; CAPSCAN_MG omitted | tests whether spectra can substitute for specimen sequencing | full 264.9M reservoir (production arm) |
| search engine settings | frozen production Sage settings, LFQ on | same settings as the production arm | none |
| matched roster | 282 specimen-metagenome matched acquisitions of 331 (347 in the rerun) | metagenome comparison needs a matched specimen | all 331 retained for other endpoints |
| evidence pool | all 331 production prediction files pooled | production rule | none |

## Known limitations

See `OPEN_ITEMS.md`.
