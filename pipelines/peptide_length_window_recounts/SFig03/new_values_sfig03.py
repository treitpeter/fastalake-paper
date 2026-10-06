"""NEW_VALUES.tsv for Supplementary Figure 3 (len950): panels c, d medians and Results line 89, 5-50 vs 9-50."""
import json
from pathlib import Path
H = Path(__file__).resolve().parent
c, n = (json.loads((H / 'work' / m / 'BUILD_LOG.json').read_text()) for m in ('control', 'len950'))
src = {'c': 'source_data/SFig03/CAMPI_reference_recovery.tsv', 'd': 'source_data/SFig03/CAMPI_within_MSFragger.tsv'}
filt = ', min_length == 9 rows (accepted 5-50 list filtered to length >= 9), median reference_recovery x 100 over 15 acquisitions'
rows = []
for a, b in zip(c['panel_medians'], n['panel_medians']):
    assert (a['panel'], a['arm']) == (b['panel'], b['arm'])
    rows.append((a['panel'], f"median reference identities recovered %, {a['arm']}", f"{a['median_percent']:.1f}", f"{b['median_percent']:.1f}", src[a['panel']] + filt))
for k, q in [('ratio_REF', 'median FL/REF peptide-count ratio (MSFragger)'), ('ratio_MG', 'median FL/MG peptide-count ratio (MSFragger)'),
             ('recovery_pooled', 'median comparator identity recovery %, REF and MG pooled (n = 30)')]:
    f = (lambda v: f'{100*v:.1f}') if k == 'recovery_pooled' else (lambda v: f'{v:.2f}')
    rows.append(('d/text 89', q, f(c['line89'][k]), f(n['line89'][k]), src['d'] + ', min_length == 9 rows; panels_len950.stat_line89_msfragger'))
rows.append(('a,b', 'pooled FL member-specific peptides (window comparison panels)', 'unchanged', 'unchanged', 'source_data/SFig03/CAMPI_member_support_recount.tsv (both windows kept)'))
with open(H / 'NEW_VALUES.tsv', 'w') as fh:
    fh.write('panel\tquantity\told_5_50\tnew_9_50\tsource_file_and_filter\n')
    for r in rows: fh.write('\t'.join(r) + '\n')
print(open(H / 'NEW_VALUES.tsv').read())
