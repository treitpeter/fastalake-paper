"""Task B: GMSC (smORF) share of selected proteins and of accepted peptides; GMSC-only peptides absent from the matched MG database.
usage: python b_gmsc.py [pilot]"""
import json, statistics, sys, os
from collections import Counter
from multiprocessing import Pool
from pathlib import Path
import ahocorasick
from common import ACC, A08, HST, PRD, HST_ROSTER, OUT, accepted, fasta, source, table, dump

PILOT = len(sys.argv) > 1 and sys.argv[1] == 'pilot'


def searched_fasta(results):
    return json.loads((Path(results).parent / 'results.json').read_text())['database']['fasta']


def tasks():
    man = json.loads((ACC / 'manifest.json').read_text()); byname = {}
    for r in man['rows']:
        byname[(r['cohort'], r['sample'])] = r
    T = []
    for s in Path(HST_ROSTER).read_text().split():
        T.append({'cohort': 'HSTOOL', 'sample': s, 'results': str(HST / '5_sage/razor_hash-acc' / s / 'results.sage.tsv'), 'mg': None, 'check9': None})
    for d in sorted((PRD / '4_parsimony/razor_hash-acc').iterdir()):
        s = d.name; r = byname.get(('PREDICT', s))
        res = r['arms']['FL_public']['results'] if r else str(PRD / '5_sage/razor_hash-acc' / s / 'results.sage.tsv')
        T.append({'cohort': 'PREDICT', 'sample': s, 'results': res, 'mg': r['arms']['MG']['fasta'] if r else None,
                  'check9': r['arms']['FL_public']['expected_sets']['9']['path'] if r else None})
    for x in json.loads((A08 / 'inputs/manifest.json').read_text())['samples']:
        s = x['sample']; r = byname.get(('CAPSCAN', s))
        T.append({'cohort': 'CAPSCAN', 'sample': s, 'results': str(A08 / 'results' / s / 'search/results.sage.tsv'), 'mg': r['arms']['MG']['fasta'] if r else None,
                  'check9': str(A08 / 'results' / s / 'length9_50.txt'), 'check9_production': r['arms']['FL_public']['expected_sets']['9']['path'] if r else None})
    return T


def work(t):
    s5, s9, nrows = accepted(t['results'])
    fa = searched_fasta(t['results'])
    A = ahocorasick.Automaton()
    for p in s5:
        A.add_word(p, p)
    A.make_automaton()
    src_n = Counter(); pep_src = {}; nprot = 0
    for acc, seq in fasta(fa):
        nprot += 1; so = source(acc); src_n[so] += 1
        for _, p in A.iter(seq.replace('I', 'L')):
            pep_src.setdefault(p, set()).add(so)
    row = {'cohort': t['cohort'], 'sample': t['sample'], 'searched_fasta': fa, 'results': t['results'], 'psm_rows_read': nrows,
           'selected_proteins': nprot, **{f'proteins_{k}': src_n.get(k, 0) for k in ['UHGG', 'GMGC', 'GMSC', 'SwissProt_human', 'OTHER']}}
    for k in [9, 5]:
        S = s9 if k == 9 else s5
        row[f'peptides_{k}_50'] = len(S)
        row[f'unmapped_{k}_50'] = sum(1 for p in S if p not in pep_src)
        row[f'any_GMSC_{k}_50'] = sum(1 for p in S if 'GMSC' in pep_src.get(p, ()))
        row[f'GMSC_only_{k}_50'] = sum(1 for p in S if pep_src.get(p) == {'GMSC'})
    if t.get('check9'):
        exp = set(Path(t['check9']).read_text().splitlines()); row['set9_equals_reference'] = exp == s9
    if t.get('check9_production'):
        row['set9_equals_production_A08'] = set(Path(t['check9_production']).read_text().splitlines()) == s9
    gonly = sorted(p for p in s5 if pep_src.get(p) == {'GMSC'})
    found = set(); mgrec = 0
    if t['mg'] and gonly:
        B = ahocorasick.Automaton()
        for p in gonly:
            B.add_word(p, p)
        B.make_automaton()
        for _, seq in fasta(t['mg']):
            mgrec += 1
            for _, p in B.iter(seq.replace('I', 'L')):
                found.add(p)
    if t['mg']:
        row['mg_fasta'] = t['mg']; row['mg_records_read'] = mgrec
        for k in [9, 5]:
            g = [p for p in gonly if len(p) >= k]
            row[f'GMSC_only_{k}_50_absent_from_MG'] = sum(1 for p in g if p not in found)
    d = OUT / 'B_gmsc_only' / t['cohort']; d.mkdir(parents=True, exist_ok=True)
    with open(d / (t['sample'] + '.tsv'), 'w') as f:
        f.write('peptide\tlength\tin_matched_MG\n')
        for p in gonly:
            f.write(f"{p}\t{len(p)}\t{'NA' if not t['mg'] else (p in found)}\n")
    return row


if __name__ == '__main__':
    T = tasks()
    print('tasks', Counter(t['cohort'] for t in T), flush=True)
    if PILOT:
        T = [next(t for t in T if t['cohort'] == c and (t['mg'] or c == 'HSTOOL')) for c in ['HSTOOL', 'PREDICT', 'CAPSCAN']]
    with Pool(int(os.environ.get('SLURM_CPUS_PER_TASK', '1'))) as pool:
        rows = []
        for i, r in enumerate(pool.imap_unordered(work, T)):
            rows.append(r); print(i + 1, len(T), r['cohort'], r['sample'], flush=True)
    rows.sort(key=lambda r: (r['cohort'], r['sample']))
    table(OUT / ('B_per_acquisition' + ('_PILOT' if PILOT else '') + '.tsv'), rows)
    if PILOT:
        print(json.dumps(rows, indent=1)); sys.exit()
    man = json.loads((ACC / 'manifest.json').read_text()); matched = {(r['cohort'], r['sample']) for r in man['rows']}
    summ = []
    for c, sub in [('HSTOOL', 'all58'), ('PREDICT', 'all306'), ('PREDICT', 'matched296'), ('CAPSCAN', 'all347'), ('CAPSCAN', 'matched282')]:
        R = [r for r in rows if r['cohort'] == c and (sub.startswith('all') or (c, r['sample']) in matched)]
        tot = lambda k: sum(r[k] for r in R)
        d = {'cohort': c, 'subset': sub, 'acquisitions': len(R), 'selected_proteins_total': tot('selected_proteins'), 'GMSC_proteins_total': tot('proteins_GMSC'),
             'GMSC_protein_share_pooled_pct': 100 * tot('proteins_GMSC') / tot('selected_proteins'),
             'GMSC_proteins_median_per_acq': statistics.median(r['proteins_GMSC'] for r in R),
             'GMSC_protein_share_median_per_acq_pct': statistics.median(100 * r['proteins_GMSC'] / r['selected_proteins'] for r in R),
             'OTHER_source_proteins_total': tot('proteins_OTHER')}
        for k in [9, 5]:
            d[f'peptides_{k}_50_total(acq-summed)'] = tot(f'peptides_{k}_50'); d[f'unmapped_{k}_50_total'] = tot(f'unmapped_{k}_50')
            d[f'GMSC_only_{k}_50_total(acq-summed)'] = tot(f'GMSC_only_{k}_50')
            d[f'GMSC_only_{k}_50_median_per_acq'] = statistics.median(r[f'GMSC_only_{k}_50'] for r in R)
            d[f'GMSC_only_{k}_50_share_median_per_acq_pct'] = statistics.median(100 * r[f'GMSC_only_{k}_50'] / r[f'peptides_{k}_50'] for r in R)
            d[f'any_GMSC_{k}_50_share_median_per_acq_pct'] = statistics.median(100 * r[f'any_GMSC_{k}_50'] / r[f'peptides_{k}_50'] for r in R)
            # distinct GMSC-only peptides across the cohort subset
            U = set(); UA = set(); UP = set()
            for r in R:
                for line in open(OUT / 'B_gmsc_only' / c / (r['sample'] + '.tsv')).read().splitlines()[1:]:
                    p, L, inmg = line.split('\t')
                    if int(L) >= k:
                        U.add(p)
                        if inmg == 'False': UA.add(p)
                        if inmg == 'True': UP.add(p)
            d[f'GMSC_only_{k}_50_distinct'] = len(U)
            if sub.startswith('matched'):
                d[f'GMSC_only_{k}_50_absent_from_MG_total(acq-summed)'] = tot(f'GMSC_only_{k}_50_absent_from_MG')
                d[f'GMSC_only_{k}_50_absent_from_MG_share_pooled_pct'] = 100 * tot(f'GMSC_only_{k}_50_absent_from_MG') / max(1, tot(f'GMSC_only_{k}_50'))
                d[f'GMSC_only_{k}_50_distinct_absent_in_every_matched_MG'] = len(UA - UP)
        summ.append(d)
    table(OUT / 'B_summary.tsv', summ); print(json.dumps(summ, indent=1))
    print('set9 reference mismatches:', [(r['cohort'], r['sample']) for r in rows if r.get('set9_equals_reference') is False])
    print('set9 production-A08 mismatches:', [(r['cohort'], r['sample']) for r in rows if r.get('set9_equals_production_A08') is False])
