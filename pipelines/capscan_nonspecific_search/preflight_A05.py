"""Preflight for the A05 participant-15 searches: every check search.py makes before Sage starts,
plus sha256 of FASTA and mzML recorded, the Sage binary version, and a negative control that
must FAIL (wrong FASTA hash; missing mzML)."""
import json, pathlib, subprocess, sys, copy, hashlib
RUN = pathlib.Path(sys.argv[1])
sys.path.insert(0, str(RUN/'frozen')); from lib import verify, sha

def check(m):
    for f in m['frozen_code']: verify(f)
    verify(m['sage'])
    for t in m['tasks']:
        assert not (RUN/'results'/t['id']/'COMPLETE.json').exists()
        for k in ['fasta', 'raw', 'source_config', 'config']: verify(t[k])
        c = json.loads(pathlib.Path(t['config']['path']).read_text())
        assert c['database']['fasta'] == t['fasta']['path'] and c['mzml_paths'] == [t['raw']['path']]
        assert c['database']['generate_decoys'] is True and c['database']['enzyme']['cleave_at'] == '' and not c['database']['enzyme'].get('semi_enzymatic')
        assert pathlib.Path(t['raw']['path']).stat().st_size > 1e8
        for k in ['baseline_search', 'semi_search']: verify(t[k])
    return True

m = json.loads((RUN/'inputs/manifest.json').read_text())
assert check(m)
v = subprocess.run([m['sage']['path'], '--version'], capture_output=True, text=True).stdout.strip()
assert v == m['actual_version'], v
neg = []
bad = copy.deepcopy(m); bad['tasks'][0]['fasta']['sha256'] = '0' * 64
try: check(bad); neg.append('wrong_fasta_hash: NOT CAUGHT')
except AssertionError: neg.append('wrong_fasta_hash: failed as required')
bad = copy.deepcopy(m); bad['tasks'][3]['raw']['path'] = bad['tasks'][3]['raw']['path'] + '.missing'
try: check(bad); neg.append('missing_mzml: NOT CAUGHT')
except (AssertionError, FileNotFoundError): neg.append('missing_mzml: failed as required')
assert all('as required' in x for x in neg), neg
hashes = {t['sample']: dict(fasta=sha(t['fasta']['path']), mzml=sha(t['raw']['path'])) for t in m['tasks']}
out = dict(status='PASS', tasks=len(m['tasks']), sage_version=v, negative_controls=neg, input_sha256=hashes)
(RUN/'PREFLIGHT.json').write_text(json.dumps(out, indent=2) + '\n'); print(json.dumps({k: out[k] for k in ['status', 'tasks', 'sage_version', 'negative_controls']}))
