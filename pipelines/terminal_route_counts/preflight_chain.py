"""Preflight for the downstream Fig. 6 stages, with negative controls that must fail."""
import json, pathlib, sys, copy
F6 = pathlib.Path(sys.argv[1]).resolve()
sys.path.insert(0, str(F6/'A05/frozen')); from lib import verify
J = lambda p: json.loads(pathlib.Path(p).read_text())
P15 = {l.split('\t')[0] for l in (F6.parent/'roster/participant15_added.tsv').read_text().splitlines()[1:]}

def per_acq(m, codekey, extra_top=()):
    for x in m[codekey]: verify(x)
    for k in extra_top:
        for x in (m[k] if isinstance(m[k], list) else [m[k]]): verify(x)
    for t in m['tasks']:
        assert t['sample'] in P15
        verify(t['fasta']); verify(t['prediction'])
        for x in t.get('inputs', []): verify(x)
        for r in t.get('routes', {}).values(): verify(r['search']); verify(r['config'])
        if 'baseline_peptides' in t: verify(t['baseline_peptides'])
        if 'raw' in t: verify(t['raw'])
        assert t['design']['participant'] == '15' and t['design']['region'] in ('duodenum', 'jejunum', 'ileum', 'colon')
    return True

def collector(m, n, pending_ok):
    for x in m['code'] + m['inputs']: verify(x)
    assert len(m['tasks']) == n and len({t['sample'] for t in m['tasks']}) == n
    for t in m['tasks']:
        recs = [t['denovo_receipt'], t['routes_receipt']] if 'denovo_receipt' in t else list(t['receipts'].values())
        for r in recs:
            p = pathlib.Path(r)
            if t['sample'] in P15 and pending_ok: assert str(F6) in str(p); continue
            assert J(p)['status'] == 'PASS', r
    assert sum(t['design']['participant'] == '15' and t['design']['region'] != 'stool' for t in m['tasks']) == (16 if n == 347 else 0)
    return True

res = {}
res['A26'] = per_acq(J(F6/'A26/inputs/manifest_revision2.json'), 'code', ['source_inputs', 'protocol'])
res['A20'] = per_acq(J(F6/'A20/inputs/manifest.json'), 'frozen_code', ['protocol', 'inputs'])
res['A27'] = per_acq(J(F6/'A27/inputs/manifest.json'), 'code', ['inputs'])
for n, d, f in [(331, 'A24_baseline331', 'manifest.json'), (347, 'A24_fixed347', 'manifest.json'),
                (331, 'A28_baseline331', 'manifest_revision4.json'), (347, 'A28_fixed347', 'manifest_revision4.json')]:
    res[d] = collector(J(F6/d/'inputs'/f), n, True)
# negative controls
neg = {}
m = J(F6/'A20/inputs/manifest.json'); bad = copy.deepcopy(m); bad['tasks'][0]['fasta']['sha256'] = '0' * 64
try: per_acq(bad, 'frozen_code'); neg['A20_wrong_fasta_hash'] = 'NOT CAUGHT'
except AssertionError: neg['A20_wrong_fasta_hash'] = 'failed as required'
bad = copy.deepcopy(m); bad['tasks'][2]['prediction']['path'] += '.missing'
try: per_acq(bad, 'frozen_code'); neg['A20_missing_prediction'] = 'NOT CAUGHT'
except (AssertionError, FileNotFoundError): neg['A20_missing_prediction'] = 'failed as required'
m = J(F6/'A28_baseline331/inputs/manifest_revision4.json'); bad = copy.deepcopy(m); bad['tasks'][5]['receipts']['available'] += '.missing'
try: collector(bad, 331, True); neg['A28_missing_receipt'] = 'NOT CAUGHT'
except (AssertionError, FileNotFoundError): neg['A28_missing_receipt'] = 'failed as required'
m = J(F6/'A24_fixed347/inputs/manifest.json'); bad = copy.deepcopy(m); bad['tasks'] = bad['tasks'][:-1]
try: collector(bad, 347, True); neg['A24_dropped_task'] = 'NOT CAUGHT'
except AssertionError: neg['A24_dropped_task'] = 'failed as required'
assert all(v == 'failed as required' for v in neg.values()), neg
out = dict(status='PASS', stages=res, negative_controls=neg)
(F6/'PREFLIGHT_CHAIN.json').write_text(json.dumps(out, indent=2) + '\n'); print(json.dumps(out))
