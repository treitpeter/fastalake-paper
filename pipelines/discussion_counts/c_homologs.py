"""Task C: for the worked HSTOOL acquisition, count proteins containing each accepted peptide (I/L-equivalent exact substring)
in the searched FASTA, the retained per-acquisition candidate FASTA, the pooled evidence FASTA and the full 233,497,702-sequence reservoir.
usage: python c_homologs.py [pilot]   (pilot = first 64 MB of each file only, for timing)"""
import os, sys, json, statistics
from collections import Counter
from multiprocessing import Pool
import ahocorasick
from common import HST, LAB, OUT, accepted, table, dump

S = '20181231_QX7_PhGe_LC95_SA_WorkflowReplicates_2'
PILOT = len(sys.argv) > 1 and sys.argv[1] == 'pilot'
RES = HST / '5_sage/razor_hash-acc' / S / 'results.sage.tsv'
DBS = {'searched': HST / '4_parsimony/razor_hash-acc' / S / 'parsimony.fasta',
       'candidate': HST / '3_persample' / S / 'sample_lake.fasta',
       'evidence': HST / '1_evidence/evidence_lake.fasta',
       'reservoir': LAB / 'canonical/publication_data/cohorts/HSTOOL/0_lake/HSTOOL_lake.fasta'}
s5, s9, NROWS = accepted(RES)
A = ahocorasick.Automaton()
for p in s5:
    A.add_word(p, p)
A.make_automaton()


def chunk(args):
    path, start, end = args
    cnt = Counter(); nrec = 0
    with open(path, 'rb') as f:
        if start:
            f.seek(start - 1); f.readline()
        pos = f.tell(); line = f.readline()
        while line and not line.startswith(b'>'):
            pos = f.tell(); line = f.readline()
        # now line is a header at pos (or EOF)
        while line and pos < end:
            seq = []; pos2 = f.tell(); nxt = f.readline()
            while nxt and not nxt.startswith(b'>'):
                seq.append(nxt.strip()); pos2 = f.tell(); nxt = f.readline()
            nrec += 1
            s = b''.join(seq).decode().replace('I', 'L')
            for p in {p for _, p in A.iter(s)}:
                cnt[p] += 1
            line, pos = nxt, pos2
    return cnt, nrec


def count(path, nproc):
    size = os.path.getsize(path)
    if PILOT:
        size = min(size, 64 << 20)
    n = max(1, min(nproc * 8, size // (8 << 20) or 1))
    bounds = [(str(path), size * i // n, size * (i + 1) // n) for i in range(n)]
    tot = Counter(); recs = 0
    with Pool(nproc) as pool:
        for c, r in pool.imap_unordered(chunk, bounds):
            tot.update(c); recs += r
    return tot, recs


if __name__ == '__main__':
    nproc = int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))
    print('PSM rows read', NROWS, 'accepted 5-50', len(s5), '9-50', len(s9), flush=True)
    C = {}; recs = {}
    for k, path in DBS.items():
        C[k], recs[k] = count(path, nproc); print(k, 'records read', recs[k], flush=True)
    rows = [{'peptide': p, 'length': len(p), **{f'n_{k}': C[k].get(p, 0) for k in DBS}} for p in sorted(s5)]
    table(OUT / ('C_per_peptide' + ('_PILOT' if PILOT else '') + '.tsv'), rows)
    # was each accepted peptide itself a de novo evidence peptide of this acquisition?
    hits = set()
    if not PILOT:
        with open(HST / '3_persample' / S / 'sample_hits.tsv') as f:
            f.readline()
            for line in f:
                hits.add(line.split('\t', 1)[0].replace('I', 'L'))
    out = {'sample': S, 'psm_rows_read': NROWS, 'records_read': recs, 'pilot': PILOT}
    for L in [9, 5]:
        R = [r for r in rows if r['length'] >= L]; out[f'{L}-50'] = {'peptides': len(R)}
        for k in DBS:
            v = sorted(r[f'n_{k}'] for r in R); q = statistics.quantiles(v, n=4) if len(v) > 1 else v
            out[f'{L}-50'][k] = {'median': statistics.median(v), 'q1': q[0], 'q3': q[-1], 'mean': sum(v) / len(v), 'max': v[-1],
                                 'n_eq0': sum(x == 0 for x in v), 'n_eq1': sum(x == 1 for x in v), 'n_ge2': sum(x >= 2 for x in v),
                                 'n_ge10': sum(x >= 10 for x in v), 'n_ge100': sum(x >= 100 for x in v), 'n_ge1000': sum(x >= 1000 for x in v)}
        if not PILOT:
            out[f'{L}-50']['accepted_that_were_denovo_evidence_peptides'] = sum(1 for r in R if r['peptide'] in hits)
    dump(OUT / ('C_RESULT' + ('_PILOT' if PILOT else '') + '.json'), out); print(json.dumps(out, indent=1))
