"""Per-acquisition medians of the GMSC-only absent-from-MG share (reads B_per_acquisition.tsv)."""
import csv, statistics
from pathlib import Path
O = Path(__file__).resolve().parent
R = [r for r in csv.DictReader(open(O / 'B_per_acquisition.tsv'), delimiter='\t') if r.get('mg_fasta')]
print('rows read (matched)', len(R), 'min MG records read', min(int(r['mg_records_read']) for r in R))
for c in ['PREDICT', 'CAPSCAN']:
    S = [r for r in R if r['cohort'] == c]
    for k in [9, 5]:
        v = [100 * int(r[f'GMSC_only_{k}_50_absent_from_MG']) / int(r[f'GMSC_only_{k}_50']) for r in S if int(r[f'GMSC_only_{k}_50'])]
        a = [int(r[f'GMSC_only_{k}_50_absent_from_MG']) for r in S]
        print(c, f'{k}-50', 'acq', len(S), 'median absent count', statistics.median(a), 'median absent share %', round(statistics.median(v), 2),
              'median absent as % of all accepted', round(statistics.median(100 * int(r[f'GMSC_only_{k}_50_absent_from_MG']) / int(r[f'peptides_{k}_50']) for r in S), 3))
