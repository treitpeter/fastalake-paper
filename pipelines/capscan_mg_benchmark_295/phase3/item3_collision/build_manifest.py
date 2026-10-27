#!/usr/bin/env python3
"""Item 3 (para 288 assembly name-collision audit): build the 295-acquisition CapScan manifest for the unchanged
sequence_accumulation worker.py (sha256 331ad8e4...).

run295/manifest.json = the 282 production CAPSCAN rows, byte-identical as JSON objects, + 13 participant-15 rows in the
same schema:
  FL_public: public347/A08_347/results/<s>/ (razor FASTA, search/results.sage.tsv, peptide_q01.txt, length9_50.txt),
             the same layout as the production A08 rows.
  MG:        fasta = mg295/db/_db/<s>.metag_hashed_plus_human.fasta; raw_fasta = assemblies_prefixed/CAPSCAN/<SAMN>.faa;
             results = mg295/mg_full/out/<s>/results.sage.tsv; expected sets = ACCEPTED_5_50.txt / ACCEPTED_9_50.txt.
  design:    subject 15, group = region, specimen_id '15|<workbook_Main>' with basis aliquot_main_id (the basis production
             uses for participant 15; each of the 13 has its own BioSample, so 13 new specimens -> 263).
neg295/manifest.json = run295 with the first new row's MG raw_fasta swapped for another new specimen's assembly and
moved to the front; the worker's 'MG raw+human set mismatch' gate must fail.
"""
import csv, hashlib, json, os
from pathlib import Path
B = Path(__file__).resolve().parent
MG295 = B.parents[1]
FINAL = MG295.parent
LAB = Path('LAB_ROOT')
ORIG = LAB / 'review/v15_integration/analyses/sequence_accumulation_2026-09-29/manifest.json'
A08 = FINAL / 'public347/A08_347/results'
RUNS = B / 'lab/review/v15_integration/analyses'


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''): h.update(b)
    return h.hexdigest()


def eset(p):
    st = os.stat(p); return {'path': str(p), 'bytes': st.st_size, 'mtime_ns': st.st_mtime_ns, 'sha256': sha(p)}


def arm(fasta, raw, results, config, mzml, e5, e9):
    st = os.stat(results)
    return {'fasta': str(fasta), 'raw_fasta': str(raw), 'results': str(results), 'results_bytes': st.st_size,
            'results_mtime_ns': st.st_mtime_ns, 'config': str(config), 'config_sha256': sha(config), 'mzml_paths': [str(mzml)],
            'expected_sets': {'5': eset(e5), '9': eset(e9)}}


m = json.loads(ORIG.read_text())
cap = [r for r in m['rows'] if r['cohort'] == 'CAPSCAN']
assert len(cap) == 282 and len({r['design']['specimen_id'] for r in cap}) == 250
new = []
for a in csv.DictReader(open(MG295 / 'inputs/ACQ13.tsv'), delimiter='\t'):
    s = a['sample']; d = A08 / s; mg = MG295 / 'mg_full/out' / s
    for p in [d / 'COMPLETE.json', mg / 'COUNTS.json']: assert p.exists(), p
    new.append({'cohort': 'CAPSCAN', 'sample': s,
                'design': {'subject_id': a['participant'], 'group': a['region'], 'aliquot': 'main',
                           'specimen_id': f"{a['participant']}|{a['workbook_Main']}", 'specimen_basis': 'aliquot_main_id'},
                'arms': {'FL_public': arm(d / f'parsimony/{s}_razor.fasta', d / f'parsimony/{s}_razor.fasta', d / 'search/results.sage.tsv',
                                          d / 'search/results.json', a['mzml_real'], d / 'peptide_q01.txt', d / 'length9_50.txt'),
                         'MG': arm(a['mg_db'], a['assembly_faa'], mg / 'results.sage.tsv', mg / 'results.json', a['mzml'],
                                   mg / 'ACCEPTED_5_50.txt', mg / 'ACCEPTED_9_50.txt')}})
assert len(new) == 13
specs = {r['design']['specimen_id'] for r in cap}
assert not specs & {r['design']['specimen_id'] for r in new} and len({r['design']['specimen_id'] for r in new}) == 13
assert not {r['arms']['MG']['raw_fasta'] for r in cap} & {r['arms']['MG']['raw_fasta'] for r in new}
base = {k: v for k, v in m.items() if k != 'rows'}
base.update(source_manifest=str(ORIG), source_manifest_sha256=sha(ORIG), cohorts={'CAPSCAN': 295},
            note='282 production CAPSCAN rows unchanged + 13 participant-15 rows (mg295/phase3/item3_collision/build_manifest.py)')
(RUNS / 'run295/manifest.json').write_text(json.dumps(dict(base, rows=cap + new), indent=1) + '\n')
bad = json.loads(json.dumps(new[0])); bad['arms']['MG']['raw_fasta'] = new[1]['arms']['MG']['raw_fasta']
(RUNS / 'neg295/manifest.json').write_text(json.dumps(dict(base, rows=[bad] + cap + new[1:]), indent=1) + '\n')
print('run295 rows', len(cap + new), 'specimens', len(specs) + 13, 'sha', sha(RUNS / 'run295/manifest.json'))
