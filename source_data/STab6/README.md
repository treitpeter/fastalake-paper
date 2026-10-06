# Source data: Supplementary Table 6

Shuffled-query control receipts. Each of the 60 acquisitions has one receipt
(`NULL_COMPARISON.json` in the original run). The three tables below are a lossless
tabular copy of the frozen 60-receipt snapshot that the published table was built from,
with the absolute file paths removed (the receipt and cohort columns are the two path
components the producer used). Extraction code: the `_work` script `tables/STab6/extract.py`
(not in this repository); it reads the snapshot, writes these files and confirmed that all
60 snapshot entries equal the receipts on disk on 6 October 2026 (Slurm job 6069108). No
clinical or personal data.

`null_receipt_rows.tsv` is the input of the table. The other two are used only by
`tables/STab6/panels.py:checks`.

## null_receipt_rows.tsv
sha256: 1469f96d71628e93bc71203785c6129ee34850cd4484288161d1d078d043c88c
rows: 4800 (one row per acquisition, seed, peptide q threshold, count level and lookup route; 80 per acquisition)
produced by: bag_mass_paired_nulls (`runtime/tools/run_bag_nulls.py`, function `report`)

| column | unit | meaning |
|---|---|---|
| receipt | id | acquisition (raw-file name, or CAMPI sample id S01..S12, S11_Fraction1..4) |
| cohort | id | CAMPI, HSTOOL (healthy stool), CAPSCAN, PREDICT |
| seed | id | shuffle seed, 1 or 2 |
| method | id | lookup route: exact, bags, mass, combined |
| q | dimensionless | peptide (or PSM) q threshold: 0.001, 0.005, 0.01, 0.02, 0.05 |
| level | id | `canonical_peptide` (distinct canonical peptides) or `PSM` (spectrum-peptide assignments) |
| accepted | count of observations | identifications at q in this route's search |
| exact_accepted | count of observations | identifications at q in the exact route's search of the same seed |
| gained | count of observations | in this route but not in the exact route |
| lost | count of observations | in the exact route but not in this route |
| net | count of observations | accepted - exact_accepted |
| ratio_to_exact | dimensionless | accepted / exact_accepted; empty when exact_accepted is 0 |
| selected_proteins | count of protein entries | protein entries selected into this route's search database |

## null_receipts.tsv
sha256: 90e99de942008bea6c3b5a42de51cb2d34a94eaa8f1830c6bfa2f15d97e36d4f
rows: 60 (one row per acquisition)
produced by: bag_mass_paired_nulls

| column | unit | meaning |
|---|---|---|
| receipt | id | acquisition |
| cohort | id | cohort |
| status | id | receipt status; all `NULL_DIAGNOSTIC_COMPLETE_REVIEW_REQUIRED` |
| workflow_fdr | id | all `UNVERIFIED` |
| same_fdr_sensitivity_claim_supported | boolean | all False |
| n_rows | count of rows | rows in the receipt (80) |

## null_receipt_results.tsv
sha256: 9ba9cc6b2a42ad48348466742d9f48eaf66d720a9b21a4e1379d6c0895600be9
rows: 720 (one row per acquisition and Sage search arm; 12 arms)
produced by: bag_mass_paired_nulls

| column | unit | meaning |
|---|---|---|
| receipt | id | acquisition |
| cohort | id | cohort |
| arm | id | exact, bags, mass, combined (original searches); shuffle_1, shuffle_2 (combined route, seeds 1 and 2); shuffle_<seed>_<route> (added null arms) |
| bytes | bytes | size of the arm's `results.sage.tsv` |
| sha256 | id | sha256 of the arm's `results.sage.tsv` |
