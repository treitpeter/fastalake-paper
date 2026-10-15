#!/usr/bin/env python3
"""(a) Matched MG FASTA = specimen assembly proteins, renamed CAPSCAN_MG_<sha256(seq)>, + human Swiss-Prot.

Exactly the production chain, per specimen:
  1. canonical/publication/arms/rehash_metag_headers.py PASS 2 for one sample: read the prefixed .faa with
     the production read_fasta (imported), keep the first record per sequence hash in file order, write
     '>CAPSCAN_MG_<sha256>' and the sequence wrapped at 60 columns. (Per-sample output depends only on
     that sample's own records; PASS 1's global map is the same pure sha256.)
  2. canonical/publication/arms/build_metag_human_fasta.py: byte-concatenate the hashed FASTA and
     databases/isoforms/human/uniprot/swissprot_human_reviewed_2026_03.fasta, adding a newline after a
     part that does not end in one; header-count check.

usage: build_mg_fasta.py control   -> rebuild the control specimen (existing 282 pair) into db/control/,
                                      require sha256(step1) == production assemblies_hashed file and
                                      sha256(step2) == production metag_hashed_human/_db file; exit 1 otherwise.
       build_mg_fasta.py build13   -> build the 13 into db/_db/, require step 1 to equal the production
                                      assemblies_hashed file for that specimen; record sha256 of every input/output.
       build_mg_fasta.py negative  -> the control rebuild with one residue changed must FAIL the sha gate.
"""
import csv, hashlib, importlib.util, json, os, shutil, sys, time
from pathlib import Path

LAB = Path('LAB_ROOT')
ARMS = LAB / 'canonical/publication/arms'
HUMAN = Path('LAB_ROOT/databases/isoforms/human/uniprot/swissprot_human_reviewed_2026_03.fasta')
HUMAN_SHA = None
W = Path(__file__).resolve().parent
INP = W.parent / 'inputs'
COHORT = 'CAPSCAN'


def load(name):
    spec = importlib.util.spec_from_file_location(name, ARMS / f'{name}.py')
    m = importlib.util.module_from_spec(spec); spec.loader.exec_module(m); return m


REHASH = load('rehash_metag_headers')
HUMBUILD = load('build_metag_human_fasta')


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


def rehash_one(src, dst, mutate=False):
    """rehash_metag_headers.py PASS 2, one sample (lines 172-182 of the production script)."""
    recs = list(REHASH.read_fasta(src))
    if mutate:   # negative control: change one residue of the first record
        h, s = recs[0]; recs[0] = (h, ('A' if s[1] != 'A' else 'G').join([s[:1], s[2:]]))
    seqs_in = {s for _, s in recs}
    by_hash = {}
    for header, s in recs:
        hx = hashlib.sha256(s.encode()).hexdigest()
        by_hash.setdefault(hx, s)
    with open(dst, 'w') as fh:
        for hx, s in by_hash.items():
            fh.write(f'>{COHORT}_MG_{hx}\n')
            for j in range(0, len(s), 60):
                fh.write(s[j:j + 60] + '\n')
    assert set(by_hash.values()) == seqs_in, 'V1 round-trip failed'
    return len(recs), len(by_hash)


def add_human(mg, dst):
    """build_metag_human_fasta.py main(), the per-sample write (lines 125-137 of the production script)."""
    tmp = str(dst) + '.part'
    with open(tmp, 'wb') as w:
        for part in (mg, HUMAN):
            with open(part, 'rb') as rd:
                shutil.copyfileobj(rd, w, length=1 << 22)
            with open(part, 'rb') as rd:
                rd.seek(-1, os.SEEK_END)
                if rd.read(1) != b'\n':
                    w.write(b'\n')
    os.replace(tmp, dst)
    mg_n = HUMBUILD.count_headers(mg); hu_n = HUMBUILD.count_headers(HUMAN); tot = HUMBUILD.count_headers(dst)
    assert tot == mg_n + hu_n, (tot, mg_n, hu_n)
    return mg_n, hu_n, tot


def build(sample, faa, prod_hashed, outdir, prod_db=None, mutate=False):
    outdir.mkdir(parents=True, exist_ok=True)
    t = time.time()
    hashed = outdir / (Path(faa).name.rsplit('.', 1)[0] + '.fasta')
    before, after = rehash_one(faa, hashed, mutate)
    db = outdir / f'{sample}.metag_hashed_plus_human.fasta'
    mg_n, hu_n, tot = add_human(hashed, db)
    rec = dict(sample=sample, assembly_faa=str(faa), assembly_faa_sha256=sha(faa), records_before=before, records_after=after,
               rebuilt_hashed=str(hashed), rebuilt_hashed_sha256=sha(hashed), production_hashed=str(prod_hashed),
               production_hashed_sha256=sha(prod_hashed), human=str(HUMAN), human_sha256=sha(HUMAN),
               mg_db=str(db), mg_db_sha256=sha(db), mg_db_bytes=db.stat().st_size, mg_entries=mg_n, human_entries=hu_n, total_entries=tot,
               seconds=time.time() - t)
    rec['hashed_identical_to_production'] = rec['rebuilt_hashed_sha256'] == rec['production_hashed_sha256']
    if prod_db:
        rec['production_mg_db'] = str(prod_db); rec['production_mg_db_sha256'] = sha(prod_db)
        rec['mg_db_identical_to_production'] = rec['mg_db_sha256'] == rec['production_mg_db_sha256']
    return rec


def main(mode):
    if mode in ('control', 'negative'):
        c = next(csv.DictReader(open(INP / 'CONTROL.tsv'), delimiter='\t'))
        out = W / ('control' if mode == 'control' else 'negative')
        rec = build(c['sample'], c['assembly_faa'], c['production_hashed_fasta'], out, c['production_mg_db'], mutate=mode == 'negative')
        ok = rec['hashed_identical_to_production'] and rec['mg_db_identical_to_production']
        rec['status'] = 'PASS' if ok else 'FAIL'
        json.dump(rec, open(out / f'{mode.upper()}.json', 'w'), indent=1)
        print(json.dumps(rec, indent=1))
        if mode == 'negative':
            # the gate must reject a database that differs by one residue
            print('NEGATIVE CONTROL', 'FAILED AS REQUIRED' if not ok else 'DID NOT FAIL')
            sys.exit(0 if not ok else 1)
        sys.exit(0 if ok else 1)
    rows = list(csv.DictReader(open(INP / 'ACQ13.tsv'), delimiter='\t'))
    only = sys.argv[2:]
    recs = []
    for r in rows:
        if only and r['sample'] not in only:
            continue
        rec = build(r['sample'], r['assembly_faa'], r['production_hashed_fasta'], W / '_db')
        assert Path(rec['mg_db']) == Path(r['mg_db']), (rec['mg_db'], r['mg_db'])
        # mzML and second-engine raw hashed here (Slurm), for the PLAN inputs table
        rec['mzml'] = r['mzml_real']; rec['mzml_sha256'] = sha(r['mzml_real'])
        rec['se_raw_sha256_recorded'] = r['se_raw_sha256']; rec['se_raw_identical'] = rec['mzml_sha256'] == r['se_raw_sha256'] if os.path.realpath(r['se_raw']) == r['mzml_real'] else None
        rec['status'] = 'PASS' if rec['hashed_identical_to_production'] else 'FAIL'
        recs.append(rec); print(json.dumps(rec), flush=True)
    name = 'BUILD13.json' if not only else f"BUILD_{'_'.join(s.rsplit('_', 1)[-1] for s in only)}.json"
    json.dump(recs, open(W / name, 'w'), indent=1)
    assert recs and all(r['status'] == 'PASS' for r in recs), 'a rebuilt hashed FASTA differs from production'


if __name__ == '__main__':
    main(sys.argv[1])
