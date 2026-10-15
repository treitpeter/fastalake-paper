#!/usr/bin/env python3
"""roster_mg295.tsv = production iteration3/disjoint/roster.tsv (578 rows, byte-identical lines) + the 13
participant-15 MG-matched acquisitions as rows i = 578..590, same columns, same derivation as
iteration3/disjoint/build_roster.py: split task from MANIFEST_347 (index 710+), odd-scan sha from the
spectrum_split347 PARITY.json, FL HALF_DISJOINT config from spectrum_split347, MG FASTA = mg295/db/_db,
MG config = the mg295 same-spectrum MG search (mg_full/out/<s>/results.json, written by job (b)).
Archived MG expected sets / frozen MG counts do not exist for the 13 (left empty; phase 2 reads mg_full)."""
import csv, json, hashlib
from pathlib import Path
M = Path(__file__).resolve().parents[1]
PROD = M.parents[1] / 'iteration3/disjoint/roster.tsv'
PROD_SHA = '33eb88f8877f260bba17b1c6ed20866c4fbebad94c49220d15d55fb26e38d5ab'
OUT = M / 'disjoint/roster_mg295.tsv'
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
assert sha(PROD) == PROD_SHA
lines = PROD.read_text().splitlines(keepends=True)
hdr = lines[0].rstrip('\n').split('\t'); assert len(lines) == 579
acq = list(csv.DictReader(open(M / 'inputs/ACQ13.tsv'), delimiter='\t'))
b13 = {r['sample']: r for r in json.load(open(M / 'db/BUILD13.json'))}
new = []
for k, a in enumerate(acq):
    s = a['sample']; db = b13[s]; assert db['status'] == 'PASS' and db['mg_db'] == a['mg_db']
    assert sha(a['mg_db']) == db['mg_db_sha256']
    new.append(dict(i=578 + k, cohort='CAPSCAN', sample=s, split_index=a['split_index'], subject_id=a['participant'], group=a['region'],
                    csv=a['csv'], mzml_split=a['mzml'], odd_mzml_sha256=a['odd_mzml_sha256'], fl_disjoint_config=a['fl_disjoint_config'],
                    mg_frozen_config=a['mg_search_dir'] + '/results.json', mg_fasta=a['mg_db'], mg_fasta_bytes=db['mg_db_bytes'],
                    mg_mzml_frozen=a['mzml'], mg_expected_5_sha256='', mg_expected_9_sha256='',
                    frozen_n_fl_5=a['fl_a08_peptides_5_50'], frozen_n_mg_5='', frozen_n_fl_9=a['fl_a08_peptides_9_50'], frozen_n_mg_9=''))
assert list(new[0]) == hdr
with open(OUT, 'w', newline='') as f:
    f.writelines(lines)
    w = csv.DictWriter(f, fieldnames=hdr, delimiter='\t', lineterminator='\n'); w.writerows(new)
ro = list(csv.DictReader(open(OUT), delimiter='\t')); assert len(ro) == 591 and [int(r['i']) for r in ro] == list(range(591))
assert open(OUT).read().startswith(PROD.read_text())
json.dump(dict(n=591, production_rows=578, new_rows=13, roster_sha256=sha(OUT), production_roster=str(PROD), production_sha256=PROD_SHA),
          open(M / 'disjoint/ROSTER.json', 'w'), indent=1)
print(sha(OUT))
