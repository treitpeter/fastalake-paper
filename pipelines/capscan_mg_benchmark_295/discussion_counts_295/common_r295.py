"""Shared readers for the Discussion counts. Read-only on every input."""
import hashlib, json, re, csv
from pathlib import Path

LAB = Path('LAB_ROOT')
OUT = Path(__file__).resolve().parent
FINAL = OUT.parent
ACC = OUT / 'trees/r295/review/v15_integration/analyses/sequence_accumulation_2026-09-29'
A08 = FINAL / 'public347/A08_347'
HST = LAB / 'canonical/publication_data/arms/noclust_2026-08-07/cohorts/HSTOOL'
PRD = LAB / 'canonical/publication_data/arms/abl_nometag_noclust/cohorts/PREDICT'
HST_ROSTER = LAB / 'simplified_grouping/publication/analyses/simplified_grouping_2026-09-10/hstool_samples.txt'
_MOD = re.compile(r"\[[^\]]*\]")


def canonical(p):
    """Identical to review/v15_integration/analyses/confidence/peptide.py::canonical."""
    return _MOD.sub('', p).replace('I', 'L')


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


def accepted(results):
    """rank 1, label 1 (target), peptide_q<=0.01, mods stripped, I->L, distinct.
    Returns (set 5-50, set 9-50, rows_read)."""
    s5 = set(); rows = 0
    with open(results) as f:
        h = f.readline().rstrip('\r\n').split('\t'); ix = {k: i for i, k in enumerate(h)}
        ip, il, iq, ir = ix['peptide'], ix['label'], ix['peptide_q'], ix['rank']
        for line in f:
            v = line.rstrip('\r\n').split('\t'); rows += 1
            if v[ir] != '1' or v[il] != '1' or float(v[iq]) > 0.01:
                continue
            p = canonical(v[ip]); assert re.fullmatch('[A-Z]+', p), p
            if 5 <= len(p) <= 50:
                s5.add(p)
    return s5, {p for p in s5 if len(p) >= 9}, rows


def fasta(path):
    acc = None; chunks = []
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if acc is not None:
                    yield acc, ''.join(chunks)
                acc = line[1:].split()[0] if line[1:].strip() else ''; chunks = []
            else:
                chunks.append(line.strip())
    if acc is not None:
        yield acc, ''.join(chunks)


def source(acc):
    """Source rule (reservoir header prefixes; see COUNTS.md section B)."""
    if acc.startswith('MGYG'):
        return 'UHGG'
    if acc.startswith('GMGC10.'):
        return 'GMGC'
    if acc.startswith('GMSC10.'):
        return 'GMSC'
    if acc.startswith('sp|'):
        return 'SwissProt_human'
    return 'OTHER'


def table(path, rows):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter='\t'); w.writeheader(); w.writerows(rows)


def dump(path, x):
    Path(path).write_text(json.dumps(x, indent=1, default=str) + '\n')
