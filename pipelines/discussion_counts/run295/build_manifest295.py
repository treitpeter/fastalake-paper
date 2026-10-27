"""Roster change only: manifest_295 = the production accumulation manifest (578 rows, every row unchanged)
+ the 13 participant-15 MG-matched acquisitions (final/mg295/inputs/ACQ13.tsv) as CAPSCAN rows.
FL_public for the 13 = public347/A08_347 (razor FASTA, Sage results, peptide_q01.txt = 5-50, length9_50.txt = 9-50).
MG for the 13 = mg295/db/_db FASTA (raw = prefixed assembly .faa), mg295/mg_full/out Sage results, ACCEPTED_{5,9}_50.txt.
Specimen = workbook_Main (design_347 rule 'main:<Main>'; no device_id recorded for participant 15), written as '15|main:<Main>'.
Also writes the negative-control manifest (13 rows, MG results of the first two swapped) and runs the roster gate
positively on manifest_295 and negatively on three wrong rosters.
usage: python build_manifest295.py"""
import csv, copy, hashlib, json
from pathlib import Path

D = Path(__file__).resolve().parent
F = D.parent
LAB = Path('LAB_ROOT')
ACC = LAB / 'review/v15_integration/analyses/sequence_accumulation_2026-09-29'
SUB = 'review/v15_integration/analyses/sequence_accumulation_2026-09-29'
A08 = F / 'public347/A08_347/results'
MGOUT = F / 'mg295/mg_full/out'
ACQ13 = F / 'mg295/inputs/ACQ13.tsv'
PROD_SHA = '02b020362f3d65708d5a187effda21e9e83973fd20e07d4420ad6df9152f1249'


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


def exp(p):
    st = Path(p).stat()
    return {'path': str(p), 'bytes': st.st_size, 'mtime_ns': st.st_mtime_ns, 'sha256': sha(p)}


def arm(fasta, raw, results, config, mzml, e5, e9):
    st = Path(results).stat()
    for p in [fasta, raw, results, config, mzml, e5, e9]:
        assert Path(p).is_file(), p
    return {'fasta': str(fasta), 'raw_fasta': str(raw), 'results': str(results), 'results_bytes': st.st_size, 'results_mtime_ns': st.st_mtime_ns,
            'config': str(config), 'config_sha256': sha(config), 'mzml_paths': [str(mzml)], 'expected_sets': {'5': exp(e5), '9': exp(e9)}}


def gate(m, prod, acq):
    """Roster gate: 296 PREDICT rows and 282 CAPSCAN rows identical to production, plus exactly the 13 ACQ13 samples."""
    rows = m['rows']; P = [r for r in rows if r['cohort'] == 'PREDICT']; C = [r for r in rows if r['cohort'] == 'CAPSCAN']
    assert len(rows) == 591 and len(P) == 296 and len(C) == 295, (len(rows), len(P), len(C))
    assert rows[:578] == prod['rows'], 'production rows changed'
    s = [r['sample'] for r in C]; assert len(set(s)) == 295
    assert set(s) == {r['sample'] for r in prod['rows'] if r['cohort'] == 'CAPSCAN'} | {a['sample'] for a in acq}
    spec = {r['design']['specimen_id'] for r in C}; assert len(spec) == 263, len(spec)
    assert len({r['design']['subject_id'] for r in C}) == 15
    for r in rows[578:]:
        a = {x['sample']: x for x in acq}[r['sample']]
        assert r['arms']['MG']['fasta'] == a['mg_db'] and r['arms']['MG']['raw_fasta'] == a['assembly_faa']
        assert r['arms']['MG']['results'] == str(MGOUT / r['sample'] / 'results.sage.tsv')
        assert r['arms']['FL_public']['results'] == str(A08 / r['sample'] / 'search/results.sage.tsv')
        assert r['design']['specimen_id'] == f"15|main:{a['workbook_Main']}" and r['design']['group'] == a['region']
    return True


prod = json.loads((ACC / 'manifest.json').read_text()); assert sha(ACC / 'manifest.json') == PROD_SHA
acq = list(csv.DictReader(open(ACQ13), delimiter='\t')); assert len(acq) == 13 and {a['participant'] for a in acq} == {'15'}
new = []
for a in acq:
    s = a['sample']; fl = A08 / s; mg = MGOUT / s
    new.append({'cohort': 'CAPSCAN', 'sample': s,
                'design': {'subject_id': '15', 'group': a['region'], 'aliquot': 'pellet', 'specimen_id': f"15|main:{a['workbook_Main']}", 'specimen_basis': 'workbook_Main'},
                'arms': {'FL_public': arm(fl / 'parsimony' / f'{s}_razor.fasta', fl / 'parsimony' / f'{s}_razor.fasta', fl / 'search/results.sage.tsv', fl / 'search/results.json',
                                          a['mzml_real'], fl / 'peptide_q01.txt', fl / 'length9_50.txt'),
                         'MG': arm(a['mg_db'], a['assembly_faa'], mg / 'results.sage.tsv', mg / 'results.json', a['mzml'], mg / 'ACCEPTED_5_50.txt', mg / 'ACCEPTED_9_50.txt')}})
m = copy.deepcopy(prod); m['rows'] = prod['rows'] + new
m['cohorts'] = {'PREDICT': 296, 'CAPSCAN': 295}
m['roster_change'] = {'base_manifest': str(ACC / 'manifest.json'), 'base_sha256': PROD_SHA, 'added_rows': 13, 'added_from': str(ACQ13), 'added_sha256': sha(ACQ13)}
assert gate(m, prod, acq)
out = D / 'trees/r295' / SUB / 'manifest.json'; out.write_text(json.dumps(m, indent=1) + '\n')

# negative controls on the roster gate: each wrong roster must be rejected
neg = {}
bad = copy.deepcopy(m); bad['rows'] = bad['rows'][:-1]; neg['drop_one_of_13'] = bad
bad = copy.deepcopy(m); bad['rows'][-1]['design']['specimen_id'] = bad['rows'][-2]['design']['specimen_id']; neg['merge_two_specimens'] = bad
bad = copy.deepcopy(m); bad['rows'][-1]['arms']['MG'], bad['rows'][-2]['arms']['MG'] = bad['rows'][-2]['arms']['MG'], bad['rows'][-1]['arms']['MG']; neg['swap_MG_of_two'] = bad
gate_neg = {}
for k, b in neg.items():
    try:
        gate(b, prod, acq); gate_neg[k] = 'ACCEPTED (control FAILED)'
    except AssertionError as e:
        gate_neg[k] = 'rejected: ' + str(e)[:120]
assert all(v.startswith('rejected') for v in gate_neg.values()), gate_neg

# negative control for the producer itself: 13-row manifest with the MG arms of the first two rows swapped
n = copy.deepcopy(prod); n['rows'] = copy.deepcopy(new); n['rows'][0]['arms']['MG'], n['rows'][1]['arms']['MG'] = n['rows'][1]['arms']['MG'], n['rows'][0]['arms']['MG']
n['cohorts'] = {'CAPSCAN': 13}; n['negative_control'] = 'MG arms of rows 0 and 1 swapped; worker must fail its set assertions'
(D / 'trees/neg' / SUB / 'manifest.json').write_text(json.dumps(n, indent=1) + '\n')
rec = {'manifest_295': str(out), 'manifest_295_sha256': sha(out), 'rows': len(m['rows']), 'capscan_rows': 295, 'capscan_specimens': 263,
       'gate_positive': 'PASS', 'gate_negative': gate_neg, 'neg_manifest_sha256': sha(D / 'trees/neg' / SUB / 'manifest.json')}
(D / 'MANIFEST295.json').write_text(json.dumps(rec, indent=1) + '\n'); print(json.dumps(rec, indent=1))
