"""Per-acquisition FASTA size rows with the A01 producer's own scan() (frozen/audit_v2.py, imported unchanged).
Usage: sizes.py repro | new | neg"""
import sys, json, csv, pathlib, hashlib, importlib.util
OUT = pathlib.Path(__file__).resolve().parents[1]
A01 = pathlib.Path('LAB_ROOT/MM_review_analysis/A01_4_1_database_sizes')
SRC = A01 / 'frozen/audit_v2.py'
print('audit_v2.py sha256', hashlib.sha256(SRC.read_bytes()).hexdigest())
spec = importlib.util.spec_from_file_location('audit_v2', SRC); a = importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
reservoirs = json.loads((A01 / 'inputs/audit_sources.json').read_text())['reservoirs']
human = json.loads((A01 / 'results/human_reference.json').read_text())
ARM = pathlib.Path('LAB_ROOT/canonical/publication_data/arms/noclust_2026-08-07/cohorts/CAPSCAN/4_parsimony/razor_hash-acc')

def row(sample, path=None):
    new = a.scan(path or ARM / sample / 'parsimony.fasta', True)   # same call as audit_v2.main
    return dict(cohort='CapScan', sample=sample, reservoir_entries=reservoirs['CAPSCAN']['entries'],
                below_human_reference=new['entries'] < human['entries'], **new)

def fmt(r):  # csv.DictWriter string rendering, as in audit_v2.table
    return {k: str(v) for k, v in r.items()}

mode = sys.argv[1]
orig = {r['sample']: r for r in csv.DictReader(open(A01 / 'results/acquisition_sizes.tsv'), delimiter='\t')}
man = json.loads((OUT / 'a09_run/inputs/manifest.json').read_text())['tasks']
if mode == 'repro':
    res = []
    for t in man[:3]:
        r = fmt(row(t['sample'])); same = r == orig[t['sample']]
        res.append(dict(sample=t['sample'], identical=same)); print('REPRO', t['sample'], 'identical' if same else 'DIFF')
        assert same, (r, orig[t['sample']])
    (OUT / 'checks/sizes_repro.json').write_text(json.dumps(res, indent=1) + '\n'); print('SIZES REPRO PASS')
elif mode == 'neg':
    try:
        row('NEGATIVE_CONTROL', ARM / 'does_not_exist' / 'parsimony.fasta'); print('NEG FAILED TO FAIL'); sys.exit(3)
    except FileNotFoundError as e:
        print('NEG correctly failed:', type(e).__name__, e)
elif mode == 'new':
    rows = [row(t['sample']) for t in man[3:]]
    assert len(rows) == 16
    with open(OUT / 'work/acquisition_sizes_new16_fullpath.tsv', 'w') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n'); w.writeheader(); w.writerows(rows)
    print('SIZES NEW16 WRITTEN')
