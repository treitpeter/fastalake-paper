"""Write NUMBERS.tsv and STATS.tsv of this unit from ``stats.collect`` (run once, on Slurm).

The notebook recomputes the same values and asserts that they equal these files.

Usage, from the repository root::

    python text_statistics/supplementary_methods_8_11/write_tables.py
"""
import csv
import sys
from pathlib import Path

UNIT = Path(__file__).resolve().parent
REPO = UNIT.parents[1]
sys.path[:0] = [str(REPO), str(UNIT)]
import stats  # noqa: E402

R, D = stats.collect(REPO / 'source_data', UNIT.name)
for name, rows, cols in [('NUMBERS.tsv', stats.number_rows(R), ['v28d_line', 'quote', 'printed_value', 'recomputed_value', 'code', 'source_data', 'status']),
                         ('STATS.tsv', stats.stats_rows(D), stats.STATS_COLUMNS)]:
    with open(UNIT / name, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=cols, delimiter='\t', lineterminator='\n')
        w.writeheader(); w.writerows(rows)
    print(name, len(rows), 'rows')
