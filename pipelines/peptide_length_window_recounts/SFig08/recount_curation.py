"""Recount PREDICT_curation_counts.tsv (A60) at 5-50 and 9-50 residues from the upstream search outputs.

Same sources and acceptance rule as A60 collect_parallel.py (phase 'full'):
  FL_cohortMG, MG, FL_public : results.sage.tsv named in the A16 manifest configs; rank 1, label 1 (target),
                               peptide_q <= 0.01, modifications stripped, I -> L.
  MG_own_curated             : R4K_ownmetag/<s>.peptides.txt.gz (target, peptide_q <= 0.01; the Sage
                               work directory was deleted by R4K_ownmetag.sbatch, this list is what remains),
                               except the 5 replayed samples: supplement6_repair/results/<idx>/results.sage.tsv.
9-50 = the accepted 5-50 set filtered by length of the stripped, I/L-merged sequence. No q re-estimation.
Output: curation_counts_5_9.tsv (sample, arm, peptides_5_50, peptides_9_50) and a check against the repo table.
"""
import os  # release: paths come from environment variables
import csv, gzip, json, re, sys, collections
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import pandas as pd

csv.field_size_limit(64 * 1024 * 1024)
LAB = Path(os.environ['FASTALAKE_DATA_DIR'])
S6 = LAB / 'MM_review_analysis/end_of_day_2026-09-15/supplement6_repair'
A16P = LAB / 'MM_review_analysis/recovery_2026-09-15/public_comparison_recovery/A16/inputs/manifest.json'
INV = LAB / 'reviews/supplementary_completion_2026-09-13/assignment_exports/PREDICT_296_MLAB/SEARCH_INPUTS.json'
R4K = LAB / 'canonical/final_revision/thread_selection_math/out/R4K_ownmetag'
REPO_TABLE = LAB / 'release_prep/fastalake-paper_v28/source_data/SFig08/PREDICT_curation_counts.tsv'
OUT = Path(__file__).resolve().parents[1]


def seq(p):
    return re.sub(r'\[[^\]]*\]', '', p)


def accepted_psm(p, sample):
    out = set()
    with Path(p).open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            assert Path(r['filename']).name == sample + '.mzML'
            if r['rank'] == '1' and r['label'] == '1' and float(r['peptide_q']) <= .01:
                s = seq(r['peptide']).replace('I', 'L'); assert 5 <= len(s) <= 50
                out.add(s)
    return out


def work(args):
    s, cfgs, own_task_idx = args
    rows = []
    for arm in ['FL_cohortMG', 'MG', 'FL_public']:
        a = accepted_psm(Path(cfgs[arm]).parent / 'results.sage.tsv', s)
        rows.append((s, arm, len(a), sum(len(x) >= 9 for x in a), str(Path(cfgs[arm]).parent / 'results.sage.tsv')))
    if own_task_idx is None:
        p = R4K / (s + '.peptides.txt.gz')
        with gzip.open(p, 'rt') as f:
            own = {seq(x.strip()).replace('I', 'L') for x in f if x.strip()}
        assert all(5 <= len(x) <= 50 for x in own)
    else:
        p = S6 / 'results' / str(own_task_idx) / 'results.sage.tsv'
        own = accepted_psm(p, s)
        ref = set((S6 / 'results' / str(own_task_idx) / 'accepted_canonical.txt').read_text().splitlines())
        assert own == ref, (s, len(own), len(ref))
    rows.append((s, 'MG_own_curated', len(own), sum(len(x) >= 9 for x in own), str(p)))
    return rows


def main():
    m = json.loads((S6 / 'manifest.json').read_text())
    bad = {t['sample'] for t in m['tasks'] if not t['pilot']}
    own_idx = {t['sample']: t['index'] for t in m['tasks'] if t['arm'] == 'own_curated'}
    a16 = {t['sample']: t for t in json.loads(A16P.read_text())['tasks'] if t['cohort'] == 'PREDICT'}
    matched = sorted(t['sample'] for t in json.loads(INV.read_text())); assert len(matched) == 296
    jobs = [(s, {arm: a16[s]['configs'][arm]['path'] for arm in ['FL_cohortMG', 'MG', 'FL_public']},
             own_idx[s] if s in bad else None) for s in matched]
    rows = []
    with ProcessPoolExecutor(max_workers=int(sys.argv[1]) if len(sys.argv) > 1 else 8) as ex:
        for r in ex.map(work, jobs):
            rows.extend(r)
    df = pd.DataFrame(rows, columns=['sample', 'arm', 'peptides_5_50', 'peptides_9_50', 'source'])
    df.to_csv(OUT / 'curation_counts_5_9.tsv', sep='\t', index=False)
    old = pd.read_csv(REPO_TABLE, sep='\t')
    mrg = old.merge(df, on=['sample', 'arm'], how='outer', indicator=True)
    assert (mrg._merge == 'both').all() and len(mrg) == 1184, mrg._merge.value_counts()
    nd = int((mrg.peptides != mrg.peptides_5_50).sum())
    print('rows', len(mrg), 'differences at 5-50 vs repo table:', nd)
    assert nd == 0
    assert (df.peptides_9_50 <= df.peptides_5_50).all()
    print('RECOUNT PASS')


if __name__ == '__main__':
    main()
