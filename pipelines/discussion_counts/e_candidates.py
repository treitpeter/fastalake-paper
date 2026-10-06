"""Task E: Sage 0.14.6 scored_candidates per rank-one spectrum.
E1: FL per-sample vs matched MG (manifest 578 searches = 289 pairs... see COUNTS.md), all paired spectra and paired spectra whose rank-one
    peptide (canonical, I/L, mods stripped) is 9-50 in both arms. Reproduces the existing candidates/full/per_acquisition.tsv first.
E2: HSTOOL 58: FL per-sample FASTA vs FL study union (A06 one-acquisition searches) and vs the 58-file joint union search.
usage: python e_candidates.py e1|e2"""
import csv, json, os, statistics, sys
from multiprocessing import Pool
from pathlib import Path
import numpy as np, pandas as pd
from common import ACC, LAB, OUT, canonical, table, dump

COLS = ['scannr', 'rank', 'label', 'scored_candidates', 'peptide']


def read(path, extra=()):
    df = pd.read_csv(path, sep='\t', usecols=COLS + list(extra), dtype={'scannr': str, 'peptide': str})
    n = len(df); df = df[df['rank'] == 1].copy()
    df['plen'] = df['peptide'].map(lambda p: len(canonical(p)))
    return df, n


def pair(f, g):
    f = f.set_index('scannr'); g = g.set_index('scannr'); assert f.index.is_unique and g.index.is_unique
    both = f.index.intersection(g.index); a = f.loc[both]; b = g.loc[both]
    m9 = ((a.plen >= 9) & (a.plen <= 50) & (b.plen >= 9) & (b.plen <= 50)).to_numpy()
    m5 = ((a.plen >= 5) & (a.plen <= 50) & (b.plen >= 5) & (b.plen <= 50)).to_numpy()
    x, y = a.scored_candidates.to_numpy(), b.scored_candidates.to_numpy()
    out = {'paired_spectra': len(both), 'median_a': float(np.median(x)), 'median_b': float(np.median(y)), 'median_ratio_b_over_a': float(np.median(y / x)),
           'paired_spectra_5_50': int(m5.sum()), 'paired_spectra_9_50': int(m9.sum())}
    if m9.any():
        out.update({'median_a_9_50': float(np.median(x[m9])), 'median_b_9_50': float(np.median(y[m9])), 'median_ratio_b_over_a_9_50': float(np.median(y[m9] / x[m9]))})
    return out


def e1(r):
    f, nf = read(r['arms']['FL_public']['results']); g, ng = read(r['arms']['MG']['results'])
    return {'cohort': r['cohort'], 'sample': r['sample'], 'subject_id': r['design']['subject_id'], 'rows_read_fl': nf, 'rows_read_mg': ng, **pair(f, g)}


A06 = LAB / 'MM_review_analysis/A06_7_hstool_second_pass/results'
JOINT = LAB / 'sample_specific_merge_all_raw/2026-09-11/cohorts/HSTOOL_noclust/search/results.sage.tsv'


def e2(s):
    f, nf = read(A06 / f'{s}__individual/search/results.sage.tsv'); g, ng = read(A06 / f'{s}__study_union/search/results.sage.tsv')
    return {'sample': s, 'rows_read_individual': nf, 'rows_read_union': ng, **pair(f, g)}


if __name__ == '__main__':
    nproc = int(os.environ.get('SLURM_CPUS_PER_TASK', '1')); mode = sys.argv[1]
    if mode == 'e1':
        man = json.loads((ACC / 'manifest.json').read_text())
        with Pool(nproc) as p:
            rows = list(p.imap(e1, man['rows']))
        prev = {x['sample']: x for x in csv.DictReader(open(ACC / 'candidates/full/per_acquisition.tsv'), delimiter='\t')}
        bad = [r['sample'] for r in rows if (float(prev[r['sample']]['median_candidates_fl_paired']), float(prev[r['sample']]['median_candidates_mg_paired']), int(prev[r['sample']]['paired_reported_spectra'])) != (r['median_a'], r['median_b'], r['paired_spectra'])]
        print('acquisitions', len(rows), 'mismatches vs existing per_acquisition.tsv', len(bad), bad[:5])
        table(OUT / 'E1_per_acquisition.tsv', rows); summ = []
        for c in ['PREDICT', 'CAPSCAN']:
            R = [r for r in rows if r['cohort'] == c]; d = {'cohort': c, 'acquisitions': len(R)}
            for k in ['median_a', 'median_b', 'median_ratio_b_over_a', 'median_a_9_50', 'median_b_9_50', 'median_ratio_b_over_a_9_50']:
                d['cohort_median_of_' + k] = statistics.median(r[k] for r in R)
            d['spectra_paired_total'] = sum(r['paired_spectra'] for r in R); d['spectra_paired_5_50_total'] = sum(r['paired_spectra_5_50'] for r in R)
            d['spectra_paired_9_50_total'] = sum(r['paired_spectra_9_50'] for r in R)
            d['psm_rows_read_total'] = sum(r['rows_read_fl'] + r['rows_read_mg'] for r in R)
            summ.append(d)
        table(OUT / 'E1_summary.tsv', summ); print(json.dumps(summ, indent=1)); dump(OUT / 'E1_RESULT.json', {'mismatches': bad, 'summary': summ})
    else:
        S = sorted({p.name.split('__')[0] for p in A06.iterdir()}); assert len(S) == 58
        with Pool(nproc) as p:
            rows = list(p.imap(e2, S))
        # joint 58-file union search: key (filename stem, scannr)
        J = []; nrows = 0
        for ch in pd.read_csv(JOINT, sep='\t', usecols=COLS + ['filename'], dtype={'scannr': str, 'peptide': str}, chunksize=2_000_000):
            nrows += len(ch); ch = ch[ch['rank'] == 1]; J.append(ch[['filename', 'scannr', 'scored_candidates', 'peptide']])
        J = pd.concat(J); J['plen'] = J['peptide'].map(lambda p: len(canonical(p))); print('joint rows read', nrows, 'rank1', len(J), flush=True)
        for r in rows:
            s = r['sample']; f, _ = read(A06 / f'{s}__individual/search/results.sage.tsv')
            g = J[J['filename'] == s + '.mzML'].drop(columns='filename')
            j = pair(f, g); r.update({'joint_' + k: v for k, v in j.items()})
        table(OUT / 'E2_hstool_per_acquisition.tsv', rows)
        d = {'acquisitions': len(rows), 'joint_psm_rows_read': nrows}
        for k in ['median_a', 'median_b', 'median_ratio_b_over_a', 'median_a_9_50', 'median_b_9_50', 'median_ratio_b_over_a_9_50']:
            d['individual_vs_union_' + k] = statistics.median(r[k] for r in rows)
            d['individual_vs_joint_' + k] = statistics.median(r['joint_' + k] for r in rows)
        d['paired_spectra_total'] = sum(r['paired_spectra'] for r in rows); d['joint_paired_spectra_total'] = sum(r['joint_paired_spectra'] for r in rows)
        dump(OUT / 'E2_RESULT.json', d); print(json.dumps(d, indent=1))
