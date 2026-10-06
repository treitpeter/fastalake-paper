# Supplementary Table 6

Shuffled-query controls. For each acquisition, the de novo queries were shuffled with two
seeds and each shuffled query set went through database selection and the Sage search by
four lookup routes: exact, recorded order bags, mass-equivalent substitutions, and the
combined route. The table counts canonical peptide observations accepted at peptide
q ≤ 0.01 and the protein entries each route selected. It is a frozen snapshot of the first
60 completed acquisition receipts (15 CAMPI, 43 healthy stool, 1 CapScan, 1 PREDICT), not
the full planned set; see `OPEN_ITEMS.md`.

## Panels

| panel | what it shows | claim (v28d line) | n and unit | test | notebook cell | source_data file | upstream pipeline |
|---|---|---|---|---|---|---|---|
| Summary sheet | accepted observations and selected protein entries per route, summed over acquisition-seed runs | exact and bag 0, mass 13,512, combined 13,518; entries 55,229 / 78,239 / 561,126 / 561,331 (485); same mass count (561) | 120 acquisition-seed runs per route (60 acquisitions x 2 seeds) | none (counts) | Table | `source_data/STab6/null_receipt_rows.tsv` | bag_mass_paired_nulls (not in `pipelines/`) |
| Acquisition_seed_rows sheet | the 480 rows the summary sums | 60 acquisitions, two seeds (485, 561) | 480 rows: 60 acquisitions x 2 seeds x 4 routes | none | Table | `source_data/STab6/null_receipt_rows.tsv` | bag_mass_paired_nulls |
| Read me sheet | scope and interpretation text | 485 | NA | none | Table | constant `README_ROWS` in `panels.py` | NA |

The notebook rebuilds both TSVs byte for byte and the workbook cell for cell, repeats the
producer checks (section `Producer checks`) and prints every number in `NUMBERS.tsv`.

## Choices

| parameter | value | why | alternatives and where reported |
|---|---|---|---|
| `Q_CUTOFF` | 0.01 | peptide q threshold used throughout the paper | each receipt also holds q = 0.001, 0.005, 0.02, 0.05 (in `null_receipt_rows.tsv`); not reported in the manuscript |
| `LEVEL` | `canonical_peptide` | counts distinct canonical peptides per acquisition and seed | PSM level is in `null_receipt_rows.tsv`; not reported |
| `METHODS` | exact, bags, mass, combined | the four lookup routes compared in Supplementary Methods 2 (line 561) | none |
| `SEEDS` | 1, 2 | the two seeded shuffles reused from the original combined null searches | none |
| `N_RECEIPTS_FROZEN` | 60 | receipts complete when the snapshot was taken (25 September 2026) | 646 of 710 planned receipts are complete now; see `OPEN_ITEMS.md` |

## Known limitations

- The 60 acquisitions are the receipts that had finished when the snapshot was frozen, not
  a designed subset; the cohort mix is uneven (15 CAMPI, 43 healthy stool, 1 CapScan,
  1 PREDICT) and the published Read me sheet does not state it.
- Accepted counts are canonical peptide-acquisition-seed observations, not globally
  distinct sequences. Selected protein entries are summed over runs, not deduplicated.
- Shuffled queries can retrieve real proteins, so accepted hits are not known false
  identifications. Every receipt records workflow FDR as unverified. These are
  nominal-threshold diagnostics.
- Acquisition identifiers are raw-file names, as in the Source Data package.
