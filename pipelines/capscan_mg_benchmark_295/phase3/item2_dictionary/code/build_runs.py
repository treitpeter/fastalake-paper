#!/usr/bin/env python3
"""Build the two grouping run folders for item 2 (metadata only: symlinks + manifests; no data read).

ctl282: the original CAPSCAN_two_arm_282 manifest with only output_root / samples_file / search_root /
        runtime_module / output relocated into phase3. search_inputs records are the ORIGINAL ones, so the
        unchanged dispatch.py preflight also proves the 282 inputs have not changed (size + mtime).
r295:   the same, plus the 13 participant-15 acquisitions (mg295/inputs/ACQ13.tsv):
        FL    = canonical noclust_2026-08-07 CAPSCAN 5_sage/razor_hash-acc/<s> (the same arm/root as the 282)
        metaG = mg295/mg_full/out/<s> (same-spectrum MG search, production config copy)
        Because the 13 metaG directories are not under the production metaG root, source_search_root
        points to a symlink farm (sources_<arm>/) whose links name the true directories; dispatch.py
        resolves both sides, so its link check is unchanged.
"""
import csv, hashlib, json, os
from pathlib import Path

FL_ROOT = Path('LAB_ROOT')
ORIG = FL_ROOT / 'simplified_grouping_results/CAPSCAN_two_arm_282_2026-09-11_immediate'
HERE = Path(__file__).resolve().parents[1]
MG295 = HERE.parents[1]
ACQ13 = [r['sample'] for r in csv.DictReader((MG295 / 'inputs/ACQ13.tsv').open(), delimiter='\t')]
NEW_SRC = {'FL': FL_ROOT / 'canonical/publication_data/arms/noclust_2026-08-07/cohorts/CAPSCAN/5_sage/razor_hash-acc',
           'metaG': MG295 / 'mg_full/out'}


def sha(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def rec(p):
    st = Path(p).stat()
    return dict(path=str(p), bytes=st.st_size, mtime_ns=st.st_mtime_ns)


def build(name, extra):
    run = HERE / 'runs' / name
    run.mkdir(parents=True, exist_ok=False)
    (run / 'code').mkdir()
    for f in ['study.py']:
        (run / 'code' / f).write_bytes((HERE / 'code' / f).read_bytes())
    m = json.loads((ORIG / 'manifest.json').read_text())
    names = (ORIG / 'inputs/acquisitions.txt').read_text().splitlines()
    assert len(names) == 282 and not set(names) & set(extra)
    allnames = sorted(names + extra)
    (run / 'inputs').mkdir()
    (run / 'inputs/acquisitions.txt').write_text('\n'.join(allnames) + '\n' if not extra else '\n'.join(allnames) + '\n')
    if not extra:
        assert (run / 'inputs/acquisitions.txt').read_bytes() == (ORIG / 'inputs/acquisitions.txt').read_bytes()
    m.update(output_root=str(run), samples_file=str(run / 'inputs/acquisitions.txt'),
             samples_sha256=sha(run / 'inputs/acquisitions.txt'),
             runtime_module=str(run / 'code/study.py'), runtime_module_sha256=sha(run / 'code/study.py'),
             comparison_population=len(allnames), manuscript_benchmark_population=len(allnames),
             phase3_note=f'{name}: relocated copy of {ORIG}/manifest.json; 13 added acquisitions: {len(extra)}; '
                         'peptide_reference is hash-checked only by dispatch.py (not read by grouping) and is kept as the original 282 file')
    assert m['runtime_module_sha256'] == '5f69c457d3636f5a1a6c9064a11a5703faab146da6a4bf1096c63dbf0359cdcd'
    for t in m['tasks']:
        arm = t['id']
        src_orig = Path(t['source_search_root'])
        inputs = run / 'inputs' / arm
        inputs.mkdir()
        farm = run / f'sources_{arm}'
        farm.mkdir()
        for s in allnames:
            real = (src_orig / s) if s in names else (NEW_SRC[arm] / s)
            assert real.is_dir(), real
            (farm / s).symlink_to(real)
            (inputs / s).symlink_to(real)
        new_records = []
        for s in extra:
            d = NEW_SRC[arm] / s
            cfg = json.loads((d / 'results.json').read_text())
            for p in [d / 'results.json', d / 'results.sage.tsv', d / 'lfq.tsv', Path(cfg['database']['fasta'])]:
                new_records.append(rec(p))
        t.update(search_root=str(inputs), source_search_root=str(farm), output=str(run / 'cohorts' / arm / 'direct'),
                 acquisitions=len(allnames), search_inputs=t['search_inputs'] + new_records)
    (run / 'manifest.json').write_text(json.dumps(m, indent=1) + '\n')
    print(name, len(allnames), 'acquisitions; manifest written')


if __name__ == '__main__':
    assert len(ACQ13) == 13
    build('ctl282', [])
    build('r295', ACQ13)
