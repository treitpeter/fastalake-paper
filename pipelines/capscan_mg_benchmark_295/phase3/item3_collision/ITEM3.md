# Item 3: assembly name-collision audit (para 288) at 295 CapScan acquisitions / 263 specimens. DONE, 7 Oct 2026

Producer: `LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/worker.py` (sha256 331ad8e4…, = release `pipelines/metagenome_header_audit/worker.py`, function `audit`), output `full/CAPSCAN/header_audit.tsv` -> Source Data `supplementary_methods_1_5/header_audit.tsv`. Para 288 uses namespace `MG_assembly_local_headers` (alias count / distinct sequences; collisions / shared headers) and `MG_as_searched` (the "53").

Method: the worker is copied byte-identically into `lab/review/v15_integration/analyses/{run282,run295,neg295}/` (a mirror tree, so its `H.parent/confidence` and `H.parents[3]/MM_review_analysis` resolve through symlinks to the real lab paths). run282 uses the production manifest byte-identically (sha 02b02036…). run295 manifest (`build_manifest.py`, sha 69e1d874…) = the 282 production CAPSCAN rows unchanged + 13 participant-15 rows: FL_public from public347/A08_347, MG = mg295/db/_db FASTA, raw = assemblies_prefixed SAMN27256486–98, results/expected sets = mg295/mg_full/out. Specimen ID for the 13 = `15|<workbook_Main>` (basis aliquot_main_id, as production uses for participant 15); 13 distinct new specimens, 263 total. The worker's own gates passed for all 295 rows: FL/MG peptide sets reproduced at 5 and 9 residues, and MG FASTA = raw assembly ∪ human Swiss-Prot.

| claim (para 288, CapScan) | 282 value | 295 value | producer file | control result | job IDs |
|---|---|---|---|---|---|
| distinct sequences in >1 specimen under different names | 2,327,521 of 28,527,640 (8.2%) | 2,393,928 of 29,163,863 (8.2%) | worker.py `audit`, MG_assembly_local_headers | run282 header_audit.tsv byte-identical to production (d685580b…); also per_acquisition, accumulation, per_specimen, MG_SEQUENCE_PRESERVATION identical | control 6085730; 295 run 6086073 |
| shared assembly-local names denoting different sequences | 2,414,558 of 2,414,668 (>99.9%) | 2,494,126 of 2,494,236 (>99.9%, 99.996%) | same | same | same |
| sequences with >1 accession after renaming (as searched) | 53 | 53 (unchanged; does not depend on the roster here) | MG_as_searched | same | same |
| accessions denoting two sequences after renaming | 0 | 0 | MG_as_searched | same | same |

Negative control: neg295 (first new row's raw assembly swapped for another specimen's) failed at the worker's gate with `AssertionError: ('…_1800', 'MG raw+human set mismatch')`, no output (job 6085729).

Memory: control MaxRSS 17.85 GiB (6085730, 16:21); 295 run requested 20G = 17.85 × 295/282 × 1.05; measured 18.35 GiB (16:56). Partition p.hpcl92, --nice=10000. Manifest build 6085728. Collection: `srun … python compare.py && python make_edits.py` (logs/collect.log) -> RESULT.json, item3_edits.jsonl.

Edit: `item3_edits.jsonl`, 1 edit, anchored once in para 288. The "53 … in each cohort" wording stays correct. Not edited (outside scope): release README `pipelines/metagenome_header_audit` and Source Data `header_audit.tsv` still state 250/282; the 295 table is `lab/review/v15_integration/analyses/run295/full/CAPSCAN/header_audit.tsv` (PREDICT rows unchanged).
