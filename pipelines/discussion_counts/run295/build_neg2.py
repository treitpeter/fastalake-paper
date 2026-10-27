"""Negative control 2: the 13 new rows, with only the MG Sage results file (path, bytes, mtime) of rows 0 and 1 swapped;
FASTA and expected peptide sets stay with their acquisition. The worker must fail its peptide-set assertion."""
import json, copy
from pathlib import Path
D = Path(__file__).resolve().parent; SUB = 'review/v15_integration/analyses/sequence_accumulation_2026-09-29'
m = json.loads((D / 'trees/r295' / SUB / 'manifest.json').read_text())
n = copy.deepcopy(m); n['rows'] = copy.deepcopy(m['rows'][578:]); assert len(n['rows']) == 13
a, b = n['rows'][0]['arms']['MG'], n['rows'][1]['arms']['MG']
for k in ['results', 'results_bytes', 'results_mtime_ns']:
    a[k], b[k] = b[k], a[k]
n['cohorts'] = {'CAPSCAN': 13}; n['negative_control'] = 'MG results of rows 0 and 1 swapped; expected sets kept'
(D / 'trees/neg2' / SUB / 'manifest.json').write_text(json.dumps(n, indent=1) + '\n'); print('neg2 written', n['rows'][0]['sample'], n['rows'][1]['sample'])
