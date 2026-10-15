#!/usr/bin/env python3
"""Item 3 collection: control (run282 header_audit.tsv byte-identical to production full/CAPSCAN/header_audit.tsv),
negative control (neg295 must have failed with the MG set-mismatch assertion), and the 295 values.
Writes RESULT.json."""
import csv, hashlib, json
from pathlib import Path
B = Path(__file__).resolve().parent
A = B / 'lab/review/v15_integration/analyses'
PROD = Path('LAB_ROOT/review/v15_integration/analyses/sequence_accumulation_2026-09-29/full/CAPSCAN')
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
rd = lambda p: {r['namespace']: r for r in csv.DictReader(open(p), delimiter='\t')}
res = {}
c = A / 'run282/full/CAPSCAN'
res['control'] = {f: {'production': sha(PROD / f), 'rerun': sha(c / f), 'identical': sha(PROD / f) == sha(c / f)}
                  for f in ['header_audit.tsv', 'MG_SEQUENCE_PRESERVATION.json', 'per_acquisition.tsv', 'accumulation.tsv', 'per_specimen.tsv']}
res['control']['header_audit_sha_expected_d685580b'] = sha(c / 'header_audit.tsv').startswith('d685580b')
neglogs = sorted(B.glob('logs/neg295_*.log'))
txt = neglogs[-1].read_text() if neglogs else ''
res['negative_control'] = {'log': str(neglogs[-1]) if neglogs else None, 'failed_on_set_mismatch': 'MG raw+human set mismatch' in txt,
                           'no_output': not (A / 'neg295/full/CAPSCAN/header_audit.tsv').exists()}
n = A / 'run295/full/CAPSCAN'
res['complete295'] = json.loads((n / 'COMPLETE.json').read_text())
old, new = rd(c / 'header_audit.tsv'), rd(n / 'header_audit.tsv')
res['values'] = {k: {'282': old[k], '295': new[k]} for k in new}
loc, srch = new['MG_assembly_local_headers'], new['MG_as_searched']
res['para288_295'] = {'alias_seqs': int(loc['shared_sequences_with_multiple_headers']), 'distinct_seqs': int(loc['distinct_sequences']),
                      'alias_pct_all': 100 * float(loc['alias_rate_among_all_sequences']),
                      'collision_headers': int(loc['shared_headers_with_multiple_sequences']), 'shared_headers': int(loc['headers_shared_across_specimens']),
                      'collision_pct_shared': 100 * float(loc['collision_rate_among_shared_headers']),
                      'as_searched_multi_accession_seqs': int(srch['shared_sequences_with_multiple_headers']),
                      'as_searched_header_collisions': int(srch['shared_headers_with_multiple_sequences'])}
(B / 'RESULT.json').write_text(json.dumps(res, indent=1) + '\n')
print(json.dumps({k: res[k] for k in ['control', 'negative_control', 'para288_295']}, indent=1))
