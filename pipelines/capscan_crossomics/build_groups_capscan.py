"""C2 copy of biology_refit/build_groups.py (mode full), CAPSCAN only. Edits: paths are arguments
(consensus full/CAPSCAN_noclust, grouping cohort dir, profiles dir, output dir); PREDICT loop and code/inputs
manifest gating (verify/gate/receipt of the biology_refit folder) removed. Computation unchanged.
Fix 2026-10-05 (6063599 failed, csv field limit): production never runs build_groups.py directly; it runs it through
biology_refit/large_csv_groups.py, which sets csv.field_size_limit(100000000) first (function_evidence_json reaches
~407 kB, 149 rows in 331 and 151 in 347 exceed the 131072 default). The copy omitted that launcher; the same limit is now
set here. Parsing only; computation unchanged.
usage: build_groups_capscan.py <consensus_arm_dir> <grouping_cohort_dir> <profiles_arm_dir> <out_dir>"""
import csv, gzip, hashlib, json, sys
from array import array
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse
csv.field_size_limit(100000000)
def human(t): return t['rank'] == 'species' and t['supported_lineage'][-1:] == ['s__Homo sapiens']
src, g, prof, dest = map(Path, sys.argv[1:5]); dest.mkdir(parents=True, exist_ok=False)
for p in [src/'COMPLETE.json', g/'COMPLETE.json', prof/'COMPLETE.json']:
    assert json.loads(p.read_text())['status'] == 'PASS', p
z = np.load(src/'sample_totals.npz'); samples = list(z['samples']); sx = {s: i for i, s in enumerate(samples)}
host = {}
with gzip.open(src/'peptide_evidence.tsv.gz', 'rt') as f:
    for r in csv.DictReader(f, delimiter='\t'): host[r['canonical_peptide']] = human(json.loads(r['taxonomy_json']))
rr = array('i'); cc = array('i'); vv = array('d'); hh = array('d'); groups = {}; totals = np.zeros(len(samples)); counts = np.zeros(len(samples), int)
with gzip.open(g/'direct/feature_assignments.tsv.gz', 'rt') as f:
    for r in csv.DictReader(f, delimiter='\t'):
        s = r['sample']; pep = r['canonical_peptide'].upper().replace('I', 'L'); ishost = host[pep]; v = float(r['intensity']); j = sx[s]; key = r['study_group']
        if key not in groups: groups[key] = len(groups)
        rr.append(j); cc.append(groups[key]); vv.append(v); hh.append(v if ishost else 0.); totals[j] += v; counts[j] += 1
shape = (len(samples), len(groups)); coords = (np.frombuffer(rr, dtype=np.int32), np.frombuffer(cc, dtype=np.int32))
full = sparse.coo_matrix((np.frombuffer(vv, dtype=np.float64), coords), shape=shape).tocsr()
hm = sparse.coo_matrix((np.frombuffer(hh, dtype=np.float64), coords), shape=shape).tocsr(); hm.eliminate_zeros(); other = full - hm; other.eliminate_zeros()
assert np.min(other.data, initial=0) >= -1e-5
selected = np.arange(len(samples))
assert np.array_equal(counts[selected], z['features'][selected]); assert np.allclose(totals[selected], z['intensity'][selected], rtol=1e-9)
tx = pd.read_csv(prof/'taxonomy_signal.tsv', sep='\t').set_index('sample').loc[samples]
assert np.allclose(hm.sum(1).A.ravel(), tx['Human-assigned'].to_numpy(), rtol=1e-9)
prior = pd.read_csv(g/'reaggregated/simplified_study_group_matrix.tsv.gz', sep='\t', index_col=0).fillna(0).T
names = list(groups); assert np.allclose(full.toarray(), prior.reindex(index=samples, columns=names, fill_value=0).to_numpy(), rtol=1e-9, atol=1e-5)
assert np.allclose(prior.loc[samples].sum(1), totals, rtol=1e-9)
for name, m in [('supported_human_groups', hm), ('other_unresolved_groups', other)]:
    keep = m.sum(0).A.ravel() > 0; sparse.save_npz(dest/(name + '.npz'), m[:, keep])
    (dest/(name + '_axes.json')).write_text(json.dumps(dict(samples=samples, features=list(np.array(names)[keep])), indent=2) + '\n')
pd.DataFrame(dict(sample=samples, total=totals, human=hm.sum(1).A.ravel(), other=other.sum(1).A.ravel())).to_csv(dest/'coverage.tsv', sep='\t', index=False)
(dest/'COMPLETE.json').write_text(json.dumps(dict(status='PASS', cohort='CAPSCAN', acquisitions=len(samples), full_group_matrix_reproduced=True,
    human_total_reproduced=True, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=2) + '\n')
print('CAPSCAN groups PASS', len(samples), 'acquisitions', len(groups), 'groups', flush=True)
