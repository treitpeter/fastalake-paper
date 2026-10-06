"""Per acquisition, per arm, from EXISTING outputs only (no search):
  MS2 spectra in the mzML (cvParam MS:1000511 'ms level' value 2; checked: MS1+MS2+other == <spectrumList count>)
  searched FASTA: proteins (records) and residues (sum of sequence lengths)
  Sage scored_candidates: median over rank-one reported spectra (all, and paired with the partner arm by scannr)
  accepted peptides: rank 1, label 1 (target), peptide_q <= 0.01, mods stripped, I->L, distinct; 5-50 and 9-50
  proteins with >= 1 accepted 9-50 peptide (Sage 'proteins' column, rev_ removed): count, distinct 9-50 peptides per protein
      (mean, median), sequence coverage per protein (residues covered by the accepted 9-50 peptides assigned to it,
      located in the searched FASTA sequence) -> median per acquisition
usage: python worker.py TASK_INDEX NTASKS OUTDIR"""
import json, re, subprocess, sys
from pathlib import Path
import numpy as np, pandas as pd

MOD = re.compile(r'\[[^\]]*\]')
PAIRS = [('FLpub', 'MG'), ('FLcmg', 'FLunion'), ('FLown', 'FLunion')]


def strip(p):
    return MOD.sub('', p)


def ms2(path):
    lv = {}
    p = subprocess.run(['grep', '-F', 'MS:1000511', path], capture_output=True, text=True, check=True)
    for line in p.stdout.splitlines():
        v = re.search(r'value="(\d+)"', line).group(1); lv[v] = lv.get(v, 0) + 1
    c = subprocess.run(['grep', '-m1', '-o', '<spectrumList count="[0-9]*"', path], capture_output=True, text=True, check=True).stdout
    n = int(re.search(r'(\d+)', c).group(1))
    return {'ms2_spectra': lv.get('2', 0), 'ms1_spectra': lv.get('1', 0), 'spectrumList_count': n, 'ms_level_total_equals_spectrumList': sum(lv.values()) == n}


def read(path):
    d = pd.read_csv(path, sep='\t', usecols=['peptide', 'proteins', 'scannr', 'rank', 'label', 'peptide_q', 'scored_candidates'],
                    dtype={'peptide': str, 'proteins': str, 'scannr': str})
    n = len(d); d = d[d['rank'] == 1]
    assert d.scannr.is_unique, path
    return d, n


def fasta_scan(path, want):
    n = 0; res = 0; seqs = {}; acc = None; ch = []
    def flush():
        nonlocal n, res
        if acc is None:
            return
        s = ''.join(ch); n += 1; res += len(s)
        if acc in want:
            seqs[acc] = s
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                flush(); acc = line[1:].split()[0]; ch = []
            else:
                ch.append(line.strip())
    flush()
    return n, res, seqs


def arm(info, df):
    out = {}
    acc = df[(df['label'] == 1) & (df['peptide_q'] <= 0.01)].copy()
    acc['bare'] = acc['peptide'].map(strip); acc['canon'] = acc['bare'].str.replace('I', 'L')
    assert acc['canon'].str.fullmatch('[A-Z]+').all()
    L = acc['canon'].str.len()
    s5 = set(acc.loc[(L >= 5) & (L <= 50), 'canon']); s9 = {p for p in s5 if len(p) >= 9}
    out['peptides_5_50'] = len(s5); out['peptides_9_50'] = len(s9)
    if info.get('expected'):
        for k, s in (('5', s5), ('9', s9)):
            e = set(Path(info['expected'][k]).read_text().split())
            assert e == s, (info['results'], k, len(e), len(s))
        out['reproduces_existing_sets'] = True
    a9 = acc[(L >= 9) & (L <= 50)]
    # protein -> (bare peptides for location, canonical peptides for counting)
    pb, pc = {}, {}
    for bare, canon, prots in zip(a9['bare'], a9['canon'], a9['proteins']):
        for pr in prots.split(';'):
            if pr.startswith('rev_'):
                continue
            pb.setdefault(pr, set()).add(bare); pc.setdefault(pr, set()).add(canon)
    n, res, seqs = fasta_scan(info['fasta'], set(pb))
    out['db_proteins'] = n; out['db_residues'] = res
    missing = [p for p in pb if p not in seqs]; assert not missing, (info['fasta'], missing[:5])
    cov = []; unplaced = 0
    for pr, peps in pb.items():
        s = seqs[pr]; m = np.zeros(len(s), bool)
        for p in peps:
            i = s.find(p); found = i >= 0
            if not found:  # Sage may list a protein for the I/L-equivalent form only
                s2 = s.replace('I', 'L'); p2 = p.replace('I', 'L'); i = s2.find(p2)
                while i >= 0:
                    m[i:i + len(p)] = True; found = True; i = s2.find(p2, i + 1)
            else:
                while i >= 0:
                    m[i:i + len(p)] = True; i = s.find(p, i + 1)
            unplaced += not found
        cov.append(m.mean())
    k = np.array([len(v) for v in pc.values()])
    out['proteins_with_9_50_peptide'] = len(pb)
    out['peptides_per_protein_mean'] = float(k.mean()) if len(k) else None
    out['peptides_per_protein_median'] = float(np.median(k)) if len(k) else None
    out['coverage_median'] = float(np.median(cov)) if cov else None
    out['peptide_protein_pairs_unplaced'] = unplaced
    out['median_scored_candidates_all_rank1'] = float(df['scored_candidates'].median())
    out['rank1_spectra_reported'] = len(df)
    return out


def run(r):
    row = {'cohort': r['cohort'], 'acquisition': r['acquisition'], 'mzml': r['mzml'], 'mzml_same_file_all_arms': r['mzml_same_file_all_arms'],
           'flpub_source': r.get('flpub_source', '')}
    row.update(ms2(r['mzml']))
    D = {}
    for a, info in r['arms'].items():
        df, n = read(info['results']); D[a] = df
        row[f'{a}_psm_rows_read'] = n; row[f'{a}_results'] = info['results']; row[f'{a}_fasta'] = info['fasta']
        for k, v in arm(info, df).items():
            row[f'{a}_{k}'] = v
        print(r['acquisition'], a, 'rows read', n, flush=True)
    for a, b in PAIRS:
        if a not in D or b not in D:
            continue
        x = D[a].set_index('scannr'); y = D[b].set_index('scannr'); both = x.index.intersection(y.index)
        row[f'paired_spectra_{a}_{b}'] = len(both)
        row[f'{a}_median_scored_candidates_paired_vs_{b}'] = float(x.loc[both, 'scored_candidates'].median())
        row[f'{b}_median_scored_candidates_paired_vs_{a}'] = float(y.loc[both, 'scored_candidates'].median())
    return row


if __name__ == '__main__':
    i, n, out = int(sys.argv[1]), int(sys.argv[2]), Path(sys.argv[3]); out.mkdir(parents=True, exist_ok=True)
    tasks = json.loads((Path(__file__).resolve().parent / 'tasks.json').read_text())
    rows = [run(r) for r in tasks[i::n]]
    with open(out / f'part_{i:03d}.json', 'w') as f:
        json.dump(rows, f)
    print('task', i, 'acquisitions', len(rows))
