import os
"""Independent recount of the CAMPI candidate-clustering sensitivity. Reads only; writes into OUT."""
import csv, gzip, hashlib, json, re, statistics, sys
from pathlib import Path
RUN = Path(os.environ['FASTALAKE_DATA_DIR'], 'review/v15_integration/analyses/clustering_90_95_2026-09-29')
H08CMP = Path(os.environ['FASTALAKE_DATA_DIR'], 'MM_review_analysis/H08_S02_recovery_2026-09-14/independent_full/full/identification_comparison.tsv')
OUT = Path(__file__).resolve().parent
ONLY = sys.argv[1] if len(sys.argv) > 1 else None
ARMS = ['unclustered', 'cluster_099', 'cluster_097', 'cluster_095', 'cluster_090', 'cluster_075', 'cluster_050']
OLD = {'unclustered', 'cluster_097', 'cluster_099'}
m = json.loads((RUN / 'manifest.json').read_text())
tasks = [t for t in m['tasks'] if ONLY is None or t['sample'] == ONLY]
tag = '' if ONLY is None else '_' + ONLY

def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(1 << 20), b''): h.update(b)
    return h.hexdigest()

def accepted(path):
    """rank-1 target PSMs, peptide_q<=0.01, mods stripped (letters only), I->L, 5-50 aa, distinct."""
    out = set(); n_rank_gt1 = 0
    with open(path) as f:
        hdr = f.readline().rstrip('\n').split('\t')
        ip, il, iq, ir = hdr.index('peptide'), hdr.index('label'), hdr.index('peptide_q'), hdr.index('rank')
        for line in f:
            c = line.rstrip('\n').split('\t')
            if c[il] != '1': continue
            if int(c[ir]) != 1: n_rank_gt1 += 1; continue
            if not (float(c[iq]) <= 0.01): continue
            s = ''.join(ch for ch in c[ip] if 'A' <= ch <= 'Z').replace('I', 'L')
            if 5 <= len(s) <= 50: out.add(s)
    return out, n_rank_gt1

def n_entries(p):
    n = 0
    with open(p, 'rb') as f:
        for line in f:
            if line[:1] == b'>': n += 1
    return n

def cfg_core(p):
    d = json.loads(Path(p).read_text()); d.pop('output_paths', None); fasta = d['database'].pop('fasta')
    return hashlib.sha256(json.dumps(d, sort_keys=True).encode()).hexdigest(), fasta

def linclust_params(stdout):
    first = [l for l in Path(stdout).read_text().splitlines() if l.startswith('linclust ')][0].split()
    g = lambda k: first[first.index(k) + 1]
    return {'min_seq_id': g('--min-seq-id'), 'c': g('-c'), 'cov_mode': g('--cov-mode'), 'cluster_mode': g('--cluster-mode')}

sets, rows, cfgrows = {}, [], []
for t in tasks:
    S = t['sample']; prior = Path(t['prior_search']); h04 = Path(t['prior_h04'])
    for arm in ARMS:
        d = prior / arm if arm in OLD else RUN / 'results' / S / 'searches' / arm
        tsv = d / 'results.sage.tsv'
        peps, nr = accepted(tsv); sets[S, arm] = peps
        chash, fasta = cfg_core(d / 'config.json')
        rhash, rfasta = cfg_core(d / 'results.json') if (d / 'results.json').exists() else (None, None)
        if arm == 'unclustered': lc = {}
        else:
            src = (h04 if arm in OLD else RUN / 'results' / S) / (arm + '.stdout')
            lc = linclust_params(src); lc['stdout'] = str(src)
            exp_fasta = (h04 if arm in OLD else RUN / 'results' / S) / (arm + '_parsimony') / (S + '_razor.fasta')
            lc['fasta_is_expected'] = str(Path(fasta)) == str(exp_fasta)
        cfgrows.append({'sample': S, 'arm': arm, 'tsv': str(tsv), 'tsv_sha256': sha(tsv), 'rank_gt1_target_rows': nr,
                        'config_core_sha256': chash, 'results_json_core_sha256': rhash, 'results_json_fasta_same': rfasta == fasta if rfasta else None,
                        'searched_fasta': fasta, 'entries': n_entries(fasta), **lc})
    base = sets[S, 'unclustered']
    # baseline equality checks
    saved = set((prior / 'unclustered' / 'accepted_canonical_5_50.txt').read_text().splitlines())
    for arm in ARMS:
        p = sets[S, arm]
        rows.append({'sample': S, 'arm': arm, 'peptides': len(p), 'baseline_peptides': len(base), 'ratio': len(p) / len(base),
                     'reduction_pct': 100 * (1 - len(p) / len(base)), 'lost': len(base - p), 'gained': len(p - base),
                     'baseline_equals_saved_H08_list': base == saved})

entries = {(r['sample'], r['arm']): r['entries'] for r in cfgrows}
for r in rows: r['entries'] = entries[r['sample'], r['arm']]; r['entry_ratio'] = r['entries'] / entries[r['sample'], 'unclustered']

# compare with collector table, row by row
coll = {(r['sample'], r['arm']): r for r in csv.DictReader(open(RUN / 'figures/per_acquisition_identifications.tsv'), delimiter='\t')}
h08 = {(r['sample'], r['arm']): r for r in csv.DictReader(open(H08CMP), delimiter='\t')}
diffs = []
for r in rows:
    c = coll[r['sample'], r['arm']]
    for mine, theirs in [('peptides', 'peptides'), ('baseline_peptides', 'baseline_peptides'), ('lost', 'lost'), ('gained', 'gained'), ('entries', 'protein_entries')]:
        if int(c[theirs]) != r[mine]: diffs.append(('collector', r['sample'], r['arm'], mine, r[mine], c[theirs]))
    for mine, theirs in [('ratio', 'peptide_count_ratio'), ('entry_ratio', 'protein_entry_ratio')]:
        if abs(float(c[theirs]) - r[mine]) > 1e-9: diffs.append(('collector', r['sample'], r['arm'], mine, r[mine], c[theirs]))
    if r['arm'] in OLD:
        h = h08[r['sample'], r['arm']]
        for mine, theirs in [('peptides', 'accepted'), ('lost', 'lost_vs_unclustered'), ('gained', 'gained_vs_unclustered')]:
            if int(h[theirs]) != r[mine]: diffs.append(('H08_identification_comparison', r['sample'], r['arm'], mine, r[mine], h[theirs]))

# KO recovery from the frozen fixed mapping and my sets
pko = {}
with gzip.open(RUN / 'figures/fixed_peptide_KO_mapping.tsv.gz', 'rt') as f:
    for r in csv.DictReader(f, delimiter='\t'): pko[r['peptide']] = set(filter(None, r['compatible_KOs'].split(';')))
kocoll = {(r['sample'], r['arm']): r for r in csv.DictReader(open(RUN / 'figures/per_acquisition_KO_coverage.tsv'), delimiter='\t')}
korows = []; missing_from_mapping = 0
for t in tasks:
    S = t['sample']; kb = set().union(*(pko.get(p, set()) for p in sets[S, 'unclustered']))
    for arm in ARMS:
        miss = [p for p in sets[S, arm] if p not in pko]; missing_from_mapping += len(miss)
        k = set().union(*(pko.get(p, set()) for p in sets[S, arm]))
        rec = len(k & kb) / len(kb)
        korows.append({'sample': S, 'arm': arm, 'KOs': len(k), 'baseline_KOs': len(kb), 'KO_recovery': rec, 'peptides_absent_from_mapping': len(miss)})
        c = kocoll[S, arm]
        if abs(float(c['KO_recovery']) - rec) > 1e-9 or int(c['compatible_KOs']) != len(k): diffs.append(('collector_KO', S, arm, 'KO_recovery/KOs', (rec, len(k)), (c['KO_recovery'], c['compatible_KOs'])))

def summ(rs, key, arm): return statistics.median([r[key] for r in rs if r['arm'] == arm])
summary = []
for arm in ARMS[1:]:
    sub = [r for r in rows if r['arm'] == arm]
    summary.append({'arm': arm, 'n': len(sub), 'n_ratio_lt1': sum(r['ratio'] < 1 for r in sub), 'median_reduction_pct': summ(rows, 'reduction_pct', arm),
                    'median_entry_ratio': summ(rows, 'entry_ratio', arm), 'min_entry_ratio': min(r['entry_ratio'] for r in sub), 'max_entry_ratio': max(r['entry_ratio'] for r in sub),
                    'median_KO_recovery_pct': 100 * summ(korows, 'KO_recovery', arm)})

# negative controls (only meaningful with all 15)
neg = {}
if ONLY is None:
    red = sorted((r['reduction_pct'], r['sample']) for r in rows if r['arm'] == 'cluster_095')
    vals = [v for v, _ in red]; s_low = red[0][1]  # acquisition with the smallest reduction
    planted = [r['reduction_pct'] if r['sample'] != s_low else 100.0 for r in rows if r['arm'] == 'cluster_095']
    expected = sorted(vals[1:] + [100.0])[7]
    neg['A_drop_all_peptides'] = {'arm': 'cluster_095', 'acquisition': s_low, 'median_before': statistics.median(vals), 'median_after_recount': statistics.median(planted),
                                  'expected_after (8th of 15 after shift)': expected, 'n_lt1_after': sum(v > 0 for v in planted), 'detected': statistics.median(planted) != statistics.median(vals) and statistics.median(planted) == expected}
    # B: delete one peptide from S03 cluster_090 and rerun row comparison -> must flag exactly that row
    S, arm = 'S03', 'cluster_090'; p = set(sets[S, arm]); p.pop()
    flagged = int(coll[S, arm]['peptides']) != len(p)
    neg['B_drop_one_peptide'] = {'acquisition': S, 'arm': arm, 'recount': len(p), 'collector': int(coll[S, arm]['peptides']), 'row_flagged': flagged}
    # C: change q threshold to 0.05 on one file: count must change

res = {'rows': rows, 'cfg': cfgrows, 'ko': korows, 'summary': summary, 'diffs': diffs, 'neg': neg, 'missing_from_mapping': missing_from_mapping}
(OUT / f'recount{tag}.json').write_text(json.dumps(res, indent=1, default=str))
for name, data in [('per_acquisition_recount', rows), ('config_and_entries', cfgrows), ('KO_recount', korows), ('summary', summary)]:
    with open(OUT / f'{name}{tag}.tsv', 'w') as f:
        keys = list(dict.fromkeys(k for r in data for k in r))
        w = csv.DictWriter(f, fieldnames=keys, delimiter='\t'); w.writeheader(); w.writerows(data)
print(json.dumps({'summary': summary, 'n_diffs': len(diffs), 'diffs': diffs[:50], 'neg': neg, 'missing_from_mapping': missing_from_mapping}, indent=1, default=str))
