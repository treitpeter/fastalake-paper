"""D-R04: peptide/precursor-group cross-fitting folds versus the published native-scan folds.

Question: native-scan folds withhold the scored spectrum from selection but may keep another measurement of the
same peptide in the selection evidence. Does the cross-fitted FastaLake (FL) vs matched-metagenome (MG) result
change when every group of repeated measurements is kept inside one fold?

Groups (predeclared, defined without any search/identification outcome): connected components over the MS2
spectra of one acquisition joined by
  (a) the same de novo predicted sequence (column peptide_prediction_detokenized_unmodified, the column the
      fasta_extractor reads; modifications absent, I->L), or
  (b) the same precursor: same charge, precursor m/z within 10 ppm and scan start time within 1.0 min.
Folds: components sorted by size (descending, ties by SHA-256 of seed:min scan) and assigned greedily to the fold
holding the fewest spectra (lowest index on ties), so every fold holds ~10% of MS2 and every out-of-fold selection
uses ~90% of the predictions, matching the native design's prediction depth. MS2 spectra without a prediction (none
expected) use SHA-256(seed:scan) mod 10.
Everything else is the published cross-fit, imported unchanged (final/mg295/crossfit/design: frozen fold mzML
writer, out-of-fold CSV writer, frozen extract()/select() with razor/hash-acc parsimony, frozen Sage config and
binary, pooled Sage-equivalent LDA + picked-peptide q, accepted() rule). Cohort fold caches are rebuilt from the
original reservoirs with the out-of-fold predictions of every cohort acquisition (PREDICT 306, CapScan 347; the
same pools as the native comparator: iteration4 for PREDICT, public347 + mg295 for CapScan).
Arms per fold: FL (group-fold selection) and MG on the identical held-out spectra; plus FL+E and MG+E with the
A-16/D-R03 paired entrapment (matched error control), pooled and rescored identically.
Modes: map <j> | cache <c*10+k> | run <i> | negctrl <i> | collect
"""
import sys, os, json, time, gzip, shutil, hashlib, csv, re
from pathlib import Path
import numpy as np, pandas as pd
from lxml import etree as ET

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
A16 = ROOT.parent / 'crossfit_entrapment/src'
sys.path.insert(0, str(A16))
import ecv  # noqa: E402  (sha-checks and imports the frozen cross-fit design; entrapment functions)
C = ecv.C; P = C.P; sage_lda = ecv.sage_lda; W = ecv.W
assert ecv.sha(A16 / 'ecv.py') == (ROOT / 'ECV_SHA256.txt').read_text().split()[0], 'A-16 ecv.py changed after D-R04 was frozen'

GSEED = 'FastaLake-DR04-groupfold-20261008'
PPM = 10.0; RT_MIN = 1.0; K = 10
COH = ['PREDICT', 'CAPSCAN']
csv.field_size_limit(10_000_000)


def map_tasks():
    return [t for t in C.MAN['tasks'] if t['cohort'] in COH]


def hfold(scan):
    return int.from_bytes(hashlib.sha256(f'{GSEED}:{int(scan)}'.encode()).digest()[:8], 'big') % K


def rt_of_ms2(mzml):
    out = {}; M = C.M
    for _, e in ET.iterparse(mzml, events=['end'], tag=M + 'spectrum', huge_tree=True):
        cv = {x.get('accession'): x for x in e.iter(M + 'cvParam')}
        if int(cv['MS:1000511'].get('value')) == 2:
            scan = int(re.search(r'scan=(\d+)', e.get('id')).group(1))
            st = cv['MS:1000016']; v = float(st.get('value')); u = st.get('unitAccession')
            out[scan] = v if u == 'UO:0000031' else v / 60.0 if u == 'UO:0000010' else None
            assert out[scan] is not None, ('unknown RT unit', u)
        e.clear()
        while e.getprevious() is not None:
            del e.getparent()[0]
    return out


class UF:
    def __init__(s, n): s.p = np.arange(n)
    def f(s, x):
        r = x
        while s.p[r] != r: r = s.p[r]
        while s.p[x] != r: s.p[x], x = r, s.p[x]
        return r
    def u(s, a, b):
        a, b = s.f(a), s.f(b)
        if a != b: s.p[max(a, b)] = min(a, b)


def build_map(j):
    t = map_tasks()[j]
    out = ROOT / f"maps/{t['cohort']}"; out.mkdir(parents=True, exist_ok=True)
    T0 = time.time()
    rows = []
    with open(t['csv']) as f:
        for r in csv.DictReader(f):
            rows.append((int(float(r['spec_idx'])), float(r['prec_mz']), int(float(r['prec_charge'])),
                         r['peptide_prediction_detokenized_unmodified'].strip().upper().replace('I', 'L')))
    rt = rt_of_ms2(t['mzml'])
    scans = np.array([r[0] for r in rows]); assert len(set(scans.tolist())) == len(scans)
    missing_rt = [s for s in scans if s not in rt]; assert not missing_rt, missing_rt[:5]
    n = len(rows); uf = UF(n)
    seq_first = {}; e_seq = 0
    for x, r in enumerate(rows):
        if r[3]:
            if r[3] in seq_first:
                uf.u(seq_first[r[3]], x); e_seq += 1
            else:
                seq_first[r[3]] = x
    mz = np.array([r[1] for r in rows]); ch = np.array([r[2] for r in rows]); rts = np.array([rt[s] for s in scans])
    e_prec = 0
    for c in np.unique(ch):
        idx = np.flatnonzero(ch == c); idx = idx[np.argsort(mz[idx], kind='mergesort')]
        m = mz[idx]; r_ = rts[idx]
        for a in range(len(idx)):
            lim = m[a] * (1 + PPM * 1e-6); b = a + 1
            while b < len(idx) and m[b] <= lim:
                if abs(r_[b] - r_[a]) <= RT_MIN:
                    uf.u(idx[a], idx[b]); e_prec += 1
                b += 1
    root = np.array([uf.f(x) for x in range(n)])
    comps = {}
    for x, rr in enumerate(root):
        comps.setdefault(rr, []).append(x)
    order = sorted(comps.values(), key=lambda m_: (-len(m_), hashlib.sha256(f'{GSEED}:{int(scans[m_].min())}'.encode()).hexdigest()))
    load = [0] * K; fold = np.empty(n, int); comp_id = np.empty(n, int); comp_size = np.empty(n, int)
    for cid, mem in enumerate(order):
        k = min(range(K), key=lambda z: (load[z], z)); load[k] += len(mem)
        fold[mem] = k; comp_id[mem] = cid; comp_size[mem] = len(mem)
    # guarantee: no predicted sequence spans two folds
    sf = {}
    for x, r in enumerate(rows):
        if r[3]:
            assert sf.setdefault(r[3], fold[x]) == fold[x], 'sequence split across folds'
    df = pd.DataFrame(dict(scan=scans, fold=fold, component=comp_id, component_size=comp_size))
    df.to_csv(out / f"{t['sample']}.tsv.gz", sep='\t', index=False, compression='gzip')
    sizes = np.array([len(m_) for m_ in order])
    st = dict(sample=t['sample'], cohort=t['cohort'], ms2_with_prediction=n, components=len(order), sequence_edges=e_seq, precursor_edges=e_prec,
              singletons=int((sizes == 1).sum()), size_quantiles={q: float(np.quantile(sizes, q)) for q in (0.5, 0.9, 0.99)}, largest=int(sizes.max()),
              largest_fraction=float(sizes.max() / n), spectra_in_multi=int(sizes[sizes > 1].sum()), fold_ms2=load,
              fold_fraction_min=min(load) / n, fold_fraction_max=max(load) / n, map_sha256=ecv.sha(out / f"{t['sample']}.tsv.gz"),
              seconds=time.time() - T0, utc=P.now(), ppm=PPM, rt_min=RT_MIN, seed=GSEED)
    (out / f"{t['sample']}.json").write_text(json.dumps(st, indent=1))
    print(json.dumps(st), flush=True)


def load_map(cohort, sample):
    st = json.loads((ROOT / f'maps/{cohort}/{sample}.json').read_text())
    p = ROOT / f'maps/{cohort}/{sample}.tsv.gz'; assert ecv.sha(p) == st['map_sha256']
    df = pd.read_csv(p, sep='\t')
    return dict(zip(df.scan.astype(int), df.fold.astype(int)))


def gfold(fmap):
    return lambda s: fmap[s] if s in fmap else hfold(s)


def cache(index):
    co = COH[index // 10]; k = index % 10
    out = ROOT / f'caches/{co}/fold{k}'; inp = out / 'inputs'
    tasks = [t for t in C.MAN['tasks'] if t['cohort'] == co]
    assert len(tasks) == C.MAN['cohorts'][co]['n']
    shas = {}; n = {}; T0 = time.time()
    for t in tasks:
        fm = gfold(load_map(co, t['sample']))
        p = inp / f"{t['sample']}_predictions.csv"
        cnt, _ = C.write_partition_csv(t['csv'], {'oof': p}, lambda s: ['oof'] if fm(s) != k else [])
        shas[t['sample']] = C.sha(p); n[t['sample']] = cnt['oof']
    split_s = time.time() - T0
    P.extract(inp, C.MAN['cohorts'][co]['reservoir'], out, 16)
    stats = json.loads((out / 'selected.stats.json').read_text())
    rec = dict(cohort=co, fold=k, K=K, group_seed=GSEED, reservoir=C.MAN['cohorts'][co]['reservoir'], n_acquisitions=len(tasks),
               oof_csv_sha256=shas, oof_rows=n, split_seconds=split_s, stats=stats, fasta_sha256=C.sha(out / 'selected.fasta'), utc=P.now())
    shutil.rmtree(inp); (out / 'hits.tsv').unlink()
    C.dump(out / 'COMPLETE.json', rec)
    print(json.dumps({x: v for x, v in rec.items() if x not in ('oof_csv_sha256', 'oof_rows')}), flush=True)


def native_dir(row):
    i = int(row['i'])
    if row['cohort'] == 'PREDICT':
        return ecv.IT4 / f"outputs/PREDICT/{row['sample']}"
    base = ecv.REVIEW / 'final/public347/crossfit_capscan347' if i < 578 else ecv.MG295
    return base / f"outputs/CAPSCAN/{row['sample']}"


def task_for(row):
    t = C.MAN['tasks'][int(row['split_index'])]
    assert t['cohort'] == row['cohort'] and t['sample'] == row['sample'] and t['mzml'] == row['mzml_split'] and t['csv'] == row['csv']
    return t


def run(i, neg=False):
    row = C.roster()[i]; assert int(row['i']) == i
    t = task_for(row)
    assert C.sha(t['config']) == t['config_sha256'], 'config identity'
    natJ = json.loads((native_dir(row) / 'COMPLETE.json').read_text())
    assert natJ['config'] == t['config'] and natJ['mg_fasta'] == row['mg_fasta'] and C.sha(row['mg_fasta']) == natJ['mg_fasta_sha256'], 'native comparator identity'
    fm = gfold(load_map(row['cohort'], row['sample']))
    sel_route = (lambda s: C.fold_of(s)) if neg else fm  # negative control: native-scan routing for the selection input
    out = ROOT / f"outputs/{row['cohort']}/{row['sample']}"; scr = ROOT / f"scratch/{row['cohort']}/{row['sample']}"
    if (out / 'COMPLETE.json').exists():
        print('already complete', i); return
    shutil.rmtree(scr, ignore_errors=True); scr.mkdir(parents=True); out.mkdir(parents=True, exist_ok=True)
    res = dict(mode='negctrl' if neg else 'groupfold', i=i, cohort=row['cohort'], sample=row['sample'], subject_id=row['subject_id'], group=row['group'],
               native_dir=str(native_dir(row)), group_seed=GSEED, sage_sha256=C.sha(P.SAGE), mg_fasta=row['mg_fasta'], mg_fasta_sha256=C.sha(row['mg_fasta']),
               config=t['config'], code_sha256=ecv.sha(__file__), folds={})
    T = time.time()
    rec, fold_scans = C.fold_mzml(t['mzml'], t['csv'], {k: scr / f'fold{k}.mzML' for k in range(K)}, fm, out / 'FOLD_MZML.json')
    oof = {k: scr / f'oof{k}' / f"{row['sample']}_predictions.csv" for k in range(K)}
    cnt, _ = C.write_partition_csv(t['csv'], oof, lambda s: [k for k in range(K) if sel_route(s) != k])
    seq_of = {}
    with open(t['csv']) as f:
        for r in csv.DictReader(f):
            seq_of[int(float(r['spec_idx']))] = r['peptide_prediction_detokenized_unmodified'].strip().upper().replace('I', 'L')
    # leakage guards for all folds before anything is selected or searched
    for k in range(K):
        sel_scans = C.csv_scans(oof[k])
        leaked = sel_scans & set(fold_scans[k])
        assert not leaked, f'LEAKAGE: {len(leaked)} fold-{k} spectra have predictions in the fold-{k} selection input'
        held_seq = {seq_of[s] for s in fold_scans[k] if seq_of.get(s)}
        sel_seq = {seq_of[s] for s in sel_scans if seq_of.get(s)}
        assert not (held_seq & sel_seq), f'GROUP LEAKAGE: {len(held_seq & sel_seq)} predicted sequences of fold {k} also in its selection input'
    res['leakage_guards'] = 'passed: no held-out scan and no held-out predicted sequence in any selection input'
    mg_e = ecv.make_entrapment(row['mg_fasta'], scr / 'MG_E.fasta', 'MG', out)
    width = W.tol_width(t['config'])
    psms = {a: [] for a in ['FL', 'MG', 'FL_E', 'MG_E']}
    for k in range(K):
        fk = {}
        cache_rec = json.loads((ROOT / f"caches/{row['cohort']}/fold{k}/COMPLETE.json").read_text())
        fk['oof_csv_sha256'] = C.sha(oof[k]); fk['oof_rows'] = cnt[k]
        assert fk['oof_csv_sha256'] == cache_rec['oof_csv_sha256'][row['sample']], 'out-of-fold CSV differs from the one in the fold cache'
        fd = out / f'fold{k}'
        fasta = C.select(row['sample'], oof[k], ROOT / f"caches/{row['cohort']}/fold{k}/selected.fasta", fd / 'selection')
        fk.update(W.tidy_selection(fd / 'selection'))
        fk['fl_fasta'] = str(fasta); fk['fl_fasta_sha256'] = C.sha(fasta); fk['fl_fasta_sequences'] = len(P.seqset(fasta))
        fk['fold_ms2'] = rec['outputs'][str(k)]['ms2']
        fk['FL_E_db'] = ecv.make_entrapment(fasta, scr / f'FL{k}_E.fasta', f'FL{k}', out)
        for arm, db in [('FL', fasta), ('MG', row['mg_fasta']), ('FL_E', scr / f'FL{k}_E.fasta'), ('MG_E', scr / 'MG_E.fasta')]:
            tsv = C.sage_search(t['config'], db, scr / f'fold{k}.mzML', fd / arm)
            W.tidy_search(fd / arm)
            d = read_fold_g(Path(str(tsv) + '.gz'), k, fm)
            if arm in ('FL', 'MG'):
                acc = C.accepted_set(d); fk[arm] = dict(psms=len(d), **W.counts(acc))
                W.write_set(fd / arm / 'ACCEPTED.txt', acc)
            psms[arm].append(d)
        (scr / f'fold{k}.mzML').unlink()
        res['folds'][str(k)] = fk
        C.dump(out / 'PROGRESS.json', res)
    for arm in ['FL', 'MG']:
        d, acc, naive, how = W.pool(psms[arm], width)
        res[f'{arm}_pooled'] = dict(psms=len(d), decoy_psms=int((d.label == -1).sum()), **how, **W.counts(acc),
                                    accepted_sha256=W.write_set(out / f'{arm}_POOLED_ACCEPTED.txt', acc))
    # matched error control: paired entrapment on the group-fold arms (A-16/D-R03 estimator)
    comb = scr / 'FL_all_E.fasta'
    with open(comb, 'w') as w:
        for k in range(K):
            w.write(open(scr / f'FL{k}_E.fasta').read())
    mp = out / 'FL_all_entrapment_map.tsv'
    pd.concat([pd.read_csv(out / f'FL{k}_entrapment_map.tsv', sep='\t', dtype=str) for k in range(K)]).drop_duplicates('accession').to_csv(mp, sep='\t', index=False)
    for arm, cb, m in [('FL_E', comb, mp), ('MG_E', scr / 'MG_E.fasta', out / 'MG_entrapment_map.tsv')]:
        r, a, cls, pair, table = ecv.analyse_arm(psms[arm], width, cb, m)
        res[arm] = r; table.to_csv(out / f'{arm}_canonical.tsv.gz', sep='\t', index=False, compression='gzip')
    res['mg_entrapment_db'] = mg_e
    res['total_seconds'] = time.time() - T; res['utc'] = P.now(); res['job'] = os.getenv('SLURM_JOB_ID')
    shutil.rmtree(scr)
    C.dump(out / 'COMPLETE.json', res)
    print(json.dumps({x: v for x, v in res.items() if x in ('i', 'sample', 'FL_pooled', 'MG_pooled', 'total_seconds')}), flush=True)


def read_fold_g(path, k, fm):
    d = sage_lda.read(path)
    extra = pd.read_csv(path, sep='\t', usecols=['proteins', 'posterior_error'], dtype={'proteins': str})
    d['proteins'] = extra['proteins'].values; d['sage_log10_pep'] = extra['posterior_error'].astype(np.float64).values
    sc = d.scannr.str.extract(r'scan=(\d+)')[0].astype(int)
    assert d.label.eq(-1).any() and len(d) > 0, 'missing scores or decoys'
    bad = [s for s in sc.unique() if fm(s) != k]
    assert not bad, f'spectrum from another fold entered fold {k}: {bad[:5]}'
    d['fold'] = k
    return d


B = 10000; BSEED = 20261004


def collect():
    recs = []; missing = []
    for r in C.roster():
        f = ROOT / f"outputs/{r['cohort']}/{r['sample']}/COMPLETE.json"
        if not f.exists():
            missing.append(int(r['i'])); continue
        g = json.loads(f.read_text()); nat = json.loads((Path(g['native_dir']) / 'COMPLETE.json').read_text())
        st = json.loads((ROOT / f"maps/{r['cohort']}/{r['sample']}.json").read_text())
        d = dict(i=int(r['i']), cohort=r['cohort'], sample=r['sample'], group=r['group'],
                 subject_id=r['subject_id'] if r['subject_id'] not in ('', 'nan') else 'NOID_' + r['sample'],
                 components=st['components'], largest_fraction=st['largest_fraction'], spectra_in_multi_fraction=st['spectra_in_multi'] / st['ms2_with_prediction'],
                 fold_fraction_min=st['fold_fraction_min'], fold_fraction_max=st['fold_fraction_max'],
                 oof_rows_group_median=float(np.median([g['folds'][str(k)]['oof_rows'] for k in range(K)])),
                 oof_rows_native_median=float(np.median([nat['folds'][str(k)]['oof_rows'] for k in range(K)])),
                 fl_fasta_seq_group_median=float(np.median([g['folds'][str(k)]['fl_fasta_sequences'] for k in range(K)])),
                 fl_fasta_seq_native_median=float(np.median([nat['folds'][str(k)]['fl_fasta_sequences'] for k in range(K)])),
                 fallback_group=';'.join(a for a in ['FL', 'MG'] if g[f'{a}_pooled']['rescoring'] != 'lda'),
                 fallback_native=';'.join(a for a in ['FL', 'MG'] if nat[f'{a}_pooled']['rescoring'] != 'lda'))
        ndir = Path(g['native_dir']); gdir = f.parent
        for L in ['5', '9']:
            for a in ['FL', 'MG']:
                d[f'{a}_group_{L}'] = g[f'{a}_pooled'][f'peptides_{L}_50']; d[f'{a}_native_{L}'] = nat[f'{a}_pooled'][f'peptides_{L}_50']
            gg = lambda x, y: 100 * (x / y - 1) if y > 0 else np.nan
            d[f'gain_group_{L}'] = gg(d[f'FL_group_{L}'], d[f'MG_group_{L}']); d[f'gain_native_{L}'] = gg(d[f'FL_native_{L}'], d[f'MG_native_{L}'])
            d[f'gain_diff_{L}'] = d[f'gain_group_{L}'] - d[f'gain_native_{L}']
            d[f'fl_group_vs_native_{L}'] = gg(d[f'FL_group_{L}'], d[f'FL_native_{L}'])
            ag = {p for p in (gdir / 'FL_POOLED_ACCEPTED.txt').read_text().split() if len(p) >= int(L)}
            an = {p for p in (ndir / 'FL_POOLED_ACCEPTED.txt').read_text().split() if len(p) >= int(L)}
            d[f'fl_identity_jaccard_{L}'] = len(ag & an) / max(1, len(ag | an)); d[f'fl_native_retained_{L}'] = len(ag & an) / max(1, len(an))
        for arm in ['FL_E', 'MG_E']:
            for gt in ['q', 'q_pep']:
                for L in ['5_50', '5_8', '9_50']:
                    e = g[arm]['gates'][gt][L]
                    d[f'{arm}_{gt}_{L}_paired_estimate'] = e['paired_estimate']; d[f'{arm}_{gt}_{L}_diagnostic_with_unresolved'] = e['diagnostic_with_unresolved']
        recs.append(d)
    a = pd.DataFrame(recs); a.to_csv(ROOT / 'PER_ACQUISITION.tsv', sep='\t', index=False)
    cap = a.cohort == 'CAPSCAN'
    strata = [('PREDICT', a.cohort == 'PREDICT'), ('CAPSCAN_all', cap), ('CAPSCAN_stool', cap & (a.group == 'stool')), ('CAPSCAN_intestinal', cap & (a.group != 'stool'))]
    cols = ['gain_group', 'gain_native', 'gain_diff', 'fl_group_vs_native', 'fl_identity_jaccard', 'fl_native_retained']
    ecols = [c for c in a.columns if c.endswith('_paired_estimate') or c.endswith('_diagnostic_with_unresolved')]
    rng = np.random.default_rng(BSEED); out = []
    for name, m in strata:
        g = a[m].reset_index(drop=True)
        if not len(g):
            continue
        subj = g.subject_id.values; us = np.unique(subj); members = [np.flatnonzero(subj == s) for s in us]
        draws = [np.concatenate([members[k] for k in rng.integers(0, len(us), len(us))]) for _ in range(B)]
        for c in [f'{x}_{L}' for x in cols for L in ['5', '9']] + ecols + ['largest_fraction', 'spectra_in_multi_fraction', 'oof_rows_group_median', 'oof_rows_native_median']:
            v = g[c].values.astype(float); bs = np.array([np.nanmedian(v[ix]) for ix in draws])
            out.append(dict(stratum=name, quantity=c, n=len(g), n_participants=len(us), median=float(np.nanmedian(v)), ci95_low=float(np.nanpercentile(bs, 2.5)),
                            ci95_high=float(np.nanpercentile(bs, 97.5)), n_pos=int((v > 0).sum()), n_zero=int((v == 0).sum()), n_neg=int((v < 0).sum()), n_nan=int(np.isnan(v).sum())))
        for L in ['9']:
            out.append(dict(stratum=name, quantity=f'n_acq_FL_more_than_MG_group_{L}', n=len(g), median=int((g[f'FL_group_{L}'] > g[f'MG_group_{L}']).sum())))
            out.append(dict(stratum=name, quantity=f'n_acq_FL_more_than_MG_native_{L}', n=len(g), median=int((g[f'FL_native_{L}'] > g[f'MG_native_{L}']).sum())))
    s = pd.DataFrame(out); s.to_csv(ROOT / 'SUMMARY.tsv', sep='\t', index=False)
    res = dict(endpoint='D-R04', n=len(a), missing=missing, by_cohort=a.groupby('cohort').size().to_dict(),
               headline=s[s.quantity.isin([f'{x}_{L}' for x in ['gain_group', 'gain_native', 'gain_diff', 'fl_native_retained'] for L in ['5', '9']])].to_dict('records'),
               zeros={'FL_group_9_zero': int((a.FL_group_9 == 0).sum()), 'MG_group_9_zero': int((a.MG_group_9 == 0).sum())},
               fallbacks_group=int((a.fallback_group != '').sum()))
    (ROOT / 'RESULT.json').write_text(json.dumps(res, indent=1, default=str))
    L_ = ['# D-R04 result: peptide/precursor-group folds vs native-scan folds', '', f'Acquisitions: {len(a)} {res["by_cohort"]}; missing {missing}', '',
          '| stratum | quantity | median | 95% CI | n |', '|---|---|---|---|---|']
    for r in res['headline']:
        L_.append(f"| {r['stratum']} | {r['quantity']} | {r['median']:.2f} | {r['ci95_low']:.2f} to {r['ci95_high']:.2f} | {r['n']} |")
    L_ += ['', 'Group sizes, fold balance, out-of-fold prediction depth and matched entrapment estimates: SUMMARY.tsv / PER_ACQUISITION.tsv.']
    (ROOT / 'RESULT.md').write_text('\n'.join(L_) + '\n'); print('\n'.join(L_))


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'collect':
        collect()
    else:
        arg = sys.argv[2]
        idx = int(arg) if arg.isdigit() else int(open(ROOT / arg).read().split()[int(os.environ['SLURM_ARRAY_TASK_ID'])])
        {'map': build_map, 'cache': cache, 'run': run, 'negctrl': lambda j: run(j, neg=True)}[mode](idx)
