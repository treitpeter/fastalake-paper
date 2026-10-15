"""A-16 / D-R03: paired entrapment control on the ACTUAL cross-fitted FastaLake (FL) and matched-metagenome (MG)
arms of the clinical benchmark (296 PREDICT + 295 CapScan = 591 matched acquisitions), plus the same-spectrum
all-prediction FL arm (FULL) on the same spectra.

Reused producer code (imported in place, sha256-checked, not modified):
  final/mg295/crossfit/design/{common,worker,sage_lda,sage_q}.py  (= iteration4 cross-fit design; 591-row roster)
  fold rule, fold mzML writer, frozen Sage binary and per-acquisition Sage config, accepted() rule, pooled
  Sage-equivalent LDA rescoring + picked-peptide q, Sage PSM PEP model of cf_pep.py (PEP_MODEL.json).
Inputs reused, not recomputed: the 10 cross-fitted FL fold FASTAs of every acquisition (native cross-fit outputs;
iteration4 for rows 0-577, final/mg295 for rows 578-590), the archived all-prediction FASTA, the matched MG FASTA.
The fold mzML files are regenerated with the frozen writer and must be byte-identical (sha256) to the recorded ones.

What is new: after each database is fixed, a 1:1 paired shuffled entrapment set (src/make_entrapment.cpp, the A91
tryptic-block shuffle applied to every protein) is appended; Sage then searches each fold's held-out spectra
against FL+E (cross-fitted), MG+E and FULL+E with the frozen config (only fasta/mzML/output changed), and each
arm's 10 fold tables are pooled and rescored exactly as in the primary benchmark.

Modes: python ecv.py run <i> | negctrl <i> | collect
"""
import sys, os, json, time, gzip, shutil, re, hashlib, subprocess, csv
from pathlib import Path
import numpy as np, pandas as pd

HERE = Path(__file__).resolve().parent; ROOT = HERE.parent
REVIEW = Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03')
IT4 = REVIEW / 'iteration4/crossfit'; MG295 = REVIEW / 'final/mg295/crossfit'
FROZEN = {MG295 / 'design/common.py': 'd323fe1f6908d67f26c0a58af1a4b4d45516fc254638eee6f9cb500f17c35a33',
          MG295 / 'design/worker.py': '3cd642c5b2d0f0bff1eec298745f3eab98a88ec38ec6ccbdfe7f174900f2d6ba',
          MG295 / 'design/sage_lda.py': '33c5e3d0d278c54cd70e38fe9ed1972f7968272dde93843a885d0f3d489a1aa1',
          MG295 / 'design/sage_q.py': '829a7e59bbd9dc86c8ea8647cccadfd12aa96f442b5db4ad81d108cb59548c10'}


def sha(p):
    h = hashlib.sha256()
    with open(p, 'rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


for _p, _h in FROZEN.items():
    assert sha(_p) == _h, f'frozen producer changed: {_p}'
sys.path.insert(0, str(MG295 / 'design'))
import common as C  # noqa: E402  (mg295 copy: absolute paths, 591-row roster, MANIFEST_347)
import sage_lda, sage_q  # noqa: E402
import worker as W  # noqa: E402  (tidy_search, tol_width, counts)

PEP_MODEL = dict(bins=1000, bw=1.0, monotonic=True)  # pipelines/capscan_crossfit_10fold/PEP_MODEL.json (validated, job 6057983)
SEED = 20261008  # entrapment construction seed (README)
BIN = ROOT / 'bin/make_entrapment'
ARMS = ['FL', 'MG', 'FULL']
LENS = {'5_50': (5, 50), '5_8': (5, 8), '9_50': (9, 50)}
GATES = ['q', 'q_pep']
MAN_PROD = json.loads((C.REC / 'MANIFEST.json').read_text())  # production manifest used by iteration4 (rows 0-577)


def task_for(row):
    i = int(row['i']); man = MAN_PROD if i < 578 else C.MAN
    t = man['tasks'][int(row['split_index'])]
    assert t['cohort'] == row['cohort'] and t['sample'] == row['sample'], (t['sample'], row['sample'])
    assert t['mzml'] == row['mzml_split'] and t['csv'] == row['csv']
    return t


def native_dir(row):
    base = IT4 if int(row['i']) < 578 else MG295
    return base / f"outputs/{row['cohort']}/{row['sample']}"


def make_entrapment(target, combined, tag, out):
    m = out / f'{tag}_entrapment_map.tsv'; st = out / f'{tag}_entrapment_stats.json'
    subprocess.run([str(BIN), str(target), str(combined), str(m), str(st), str(SEED)], check=True)
    s = json.loads(st.read_text()); s.update(target_sha256=sha(target), combined_sha256=sha(combined), target=str(target))
    st.write_text(json.dumps(s, indent=1))
    return s


def read_fold(path, k):
    d = sage_lda.read(path)
    extra = pd.read_csv(path, sep='\t', usecols=['proteins', 'posterior_error'], dtype={'proteins': str})
    assert len(extra) == len(d)
    d['proteins'] = extra['proteins'].values; d['sage_log10_pep'] = extra['posterior_error'].astype(np.float64).values
    sc = d.scannr.str.extract(r'scan=(\d+)')[0].astype(int)
    assert d.label.eq(-1).any() and len(d) > 0, 'missing scores or decoys'
    bad = [s for s in sc.unique() if C.fold_of(s) != k]
    assert not bad, f'spectrum from another fold entered fold {k}: {bad[:5]}'
    d['fold'] = k
    return d


def raw(pep):
    p = re.sub(r'\[[^\]]*\]', '', pep).replace('-', '')
    assert p.isalpha() and p.isupper(), pep
    return p


def log10_pep(scores, decoy):
    b, lo, step = sage_q.estimator(scores.astype(np.float64), decoy, bins=PEP_MODEL['bins'], bw_adjust=PEP_MODEL['bw'], monotonic=PEP_MODEL['monotonic'])
    pe = sage_q.posterior_error(b, lo, step, scores.astype(np.float64))
    with np.errstate(divide='ignore'):
        return np.log10(pe).astype(np.float32)


def load_seqs(fasta, wanted):
    out = {}; a = None; buf = []
    with open(fasta) as f:
        for line in f:
            line = line.rstrip('\n')
            if line.startswith('>'):
                if a in wanted: out[a] = ''.join(buf)
                a = line[1:].split()[0] if len(line) > 1 else ''; buf = []
            else:
                buf.append(line.strip())
    if a in wanted: out[a] = ''.join(buf)
    return out


def classify(d, combined_fasta, map_tsv):
    """Per canonical (mods stripped, I->L) rank-one target-label peptide: class T (target-only or shared with a
    control: reported as target by an ordinary search), E (entrapment-only, paired to exactly one canonical
    target string at the same block coordinates) or U (entrapment-only, pairing not unique)."""
    r1 = d[(d['rank'] == 1) & (d['label'] == 1)]
    accs = r1.proteins.str.split(';')
    is_e = accs.map(lambda xs: all(x.startswith('FLX') for x in xs)).values
    has_e = accs.map(lambda xs: any(x.startswith('FLX') for x in xs)).values
    canon = r1.peptide.map(C.strip).values
    cls = {}; shared = set()
    for c, e, he in zip(canon, is_e, has_e):
        if not e:
            cls[c] = 'T'
            if he: shared.add(c)
        elif c not in cls:
            cls[c] = 'E'
    ent = [(c, raw(p), xs) for c, p, xs, e in zip(canon, r1.peptide.values, accs.values, is_e) if e and cls[c] == 'E']
    ctl = {x for _, _, xs in ent for x in xs}
    m = pd.read_csv(map_tsv, sep='\t', dtype=str).set_index('accession')['source_accession']
    src = {x: m[x] for x in ctl}
    seq = load_seqs(combined_fasta, ctl | set(src.values()))
    pair_sets = {}
    for c, rp, xs in ent:
        s = pair_sets.setdefault(c, set())
        for x in xs:
            cs = seq[x]; ts = seq[src[x]]; assert len(cs) == len(ts)
            pos = cs.find(rp); found = False
            while pos >= 0:
                s.add(C.strip(ts[pos:pos + len(rp)])); found = True
                pos = cs.find(rp, pos + 1)
            assert found, (rp, x)
    pair = {}
    for c, s in pair_sets.items():
        if len(s) == 1:
            pair[c] = next(iter(s))
        else:
            cls[c] = 'U'
    return cls, pair, shared


def best_rows(d, gate_mask):
    """best (pooled q, -pooled score) per canonical peptide among rank-one target-label PSMs; and accepted set
    (any such PSM passing the gate), lengths 5-50 as the frozen rule."""
    r = d[(d['rank'] == 1) & (d['label'] == 1)].copy()
    r['canon'] = r.peptide.map(C.strip)
    r = r[r.canon.str.len().between(5, 50)]
    b = r.sort_values(['pq', 'ps'], ascending=[True, False], kind='mergesort').drop_duplicates('canon')
    best = dict(zip(b.canon, zip(b.pq.astype(float), -b.ps.astype(float))))
    acc = set(r.loc[gate_mask.loc[r.index].values, 'canon'])
    return best, acc, r


def estimate(acc, best, cls, pair, lo, hi, restrict=None):
    acc = {p for p in acc if lo <= len(p) <= hi}
    if restrict is not None:
        acc &= restrict
    T = E = U = E0 = B = J = 0
    for p in acc:
        c = cls.get(p, 'T')
        if c == 'T':
            T += 1
        elif c == 'U':
            U += 1
        else:
            E += 1; t = pair[p]
            if t not in acc or cls.get(t, 'T') != 'T':
                E0 += 1
            else:
                B += best[p] < best[t]; J += best[p] == best[t]
    N = T + E + U; A = E + E0 + 2 * B
    f = lambda x: x / N if N else None
    return dict(n_accepted=N, n_target=T, n_entrapment_paired=E, n_entrapment_unresolved=U, partner_not_accepted=E0,
                entrapment_ranks_higher=B, exact_ties=J, paired_numerator=A, paired_estimate=f(A), paired_ties_upper=f(A + 2 * J),
                diagnostic_with_unresolved=f(A + U), combined_estimate=f(2 * (E + U)), lower_fraction=f(E + U))


def analyse_arm(dfs, width, combined, map_tsv):
    d = pd.concat(dfs, ignore_index=True)
    assert d.scannr.is_unique, 'a spectrum was scored in two folds'
    q, s, how = sage_lda.pooled_q(d, width)
    d['pq'] = q; d['ps'] = s
    dec = (d.label == -1).values
    d['plp'] = log10_pep(s, dec) if how == 'lda' else np.float32(0.0)
    cls, pair, shared = classify(d, combined, map_tsv)
    masks = {'q': d.pq <= 0.01, 'q_pep': (d.pq <= 0.01) & (d.plp <= -2)}
    res = dict(rescoring=how, psms=len(d), decoy_psms=int(dec.sum()), shared_target_control_strings=len(shared),
               classes={k: int(v) for k, v in pd.Series(list(cls.values())).value_counts().items()}, gates={})
    accs = {}
    for g, mk in masks.items():
        best, acc, r1 = best_rows(d, mk)
        accs[g] = (acc, best)
        res['gates'][g] = {}
        for L, (lo, hi) in LENS.items():
            e = estimate(acc, best, cls, pair, lo, hi)
            # observation (PSM) denominator: accepted rank-one target-label PSMs in the length stratum
            ok = mk.loc[r1.index].values & r1.canon.str.len().between(lo, hi).values
            cl = r1.canon.map(lambda p: cls.get(p, 'T')).values
            n_psm = int(ok.sum()); e_psm = int((ok & (cl != 'T')).sum())
            e.update(n_psm_accepted=n_psm, n_psm_entrapment=e_psm, psm_combined_estimate=(2 * e_psm / n_psm) if n_psm else None)
            res['gates'][g][L] = e
    table = pd.DataFrame([dict(canon=p, cls=cls.get(p, 'T'), pair=pair.get(p, ''), best_q=v[0], best_score=-v[1],
                               acc_q=p in accs['q'][0], acc_q_pep=p in accs['q_pep'][0]) for p, v in accs['q'][1].items()])
    return res, accs, cls, pair, table


def run(i, neg=False):
    row = C.roster()[i]; assert int(row['i']) == i
    t = task_for(row); nat = native_dir(row)
    J = json.loads((nat / 'COMPLETE.json').read_text())
    assert J['sample'] == row['sample'] and J['K'] == 10 and J['mode'] == 'crossfit' and J['fold_seed'] == C.FOLD_SEED
    assert C.sha(t['config']) == t['config_sha256'] and J['config'] == t['config'], 'config identity'
    assert C.sha(t['fasta']) == t['fasta_sha256'], 'all-prediction FASTA identity'
    assert row['mg_fasta'] == J['mg_fasta'] and C.sha(row['mg_fasta']) == J['mg_fasta_sha256'], 'MG FASTA identity'
    assert sha(C.P.SAGE) == J['sage_sha256'], 'Sage binary identity'
    # input identity checks for all 10 cross-fitted FL fold FASTAs BEFORE anything is searched
    fl = {}
    for k in range(10):
        src = J['folds'][str((k + 1) % 10 if neg else k)]['fl_fasta']
        assert C.sha(src) == J['folds'][str(k)]['fl_fasta_sha256'], f'FL fold FASTA identity failed for fold {k}'
        fl[k] = src
    tag = 'negctrl' if neg else 'run'
    out = ROOT / f"per_acq/{row['cohort']}/{row['sample']}"; scr = ROOT / f"scratch/{row['cohort']}/{row['sample']}"
    shutil.rmtree(scr, ignore_errors=True); scr.mkdir(parents=True)
    if out.exists() and (out / 'RESULT.json').exists():
        print('already complete', i); return
    out.mkdir(parents=True, exist_ok=True)
    res = dict(i=i, cohort=row['cohort'], sample=row['sample'], subject_id=row['subject_id'], group=row['group'], native_dir=str(nat),
               entrapment_seed=SEED, code_sha256=sha(__file__), make_entrapment_sha256=sha(BIN), pep_model=PEP_MODEL, folds={})
    T0 = time.time()
    rec, _ = C.fold_mzml(t['mzml'], t['csv'], {k: scr / f'fold{k}.mzML' for k in range(10)}, C.fold_of, out / 'FOLD_MZML.json')
    for k in range(10):
        assert rec['outputs'][str(k)]['sha256'] == J['folds'][str(k)]['fold_mzml_sha256'], f'fold {k} mzML differs from the cross-fit run'
    res['fold_mzml_identical'] = True
    dbs = {'MG': make_entrapment(row['mg_fasta'], scr / 'MG_E.fasta', 'MG', out),
           'FULL': make_entrapment(t['fasta'], scr / 'FULL_E.fasta', 'FULL', out)}
    width = W.tol_width(t['config'])
    psms = {a: [] for a in ARMS}
    for k in range(10):
        fk = {}
        fk['FL_db'] = make_entrapment(fl[k], scr / f'FL{k}_E.fasta', f'FL{k}', out)
        for arm, db in [('FL', scr / f'FL{k}_E.fasta'), ('MG', scr / 'MG_E.fasta'), ('FULL', scr / 'FULL_E.fasta')]:
            T1 = time.time()
            tsv = C.sage_search(t['config'], db, scr / f'fold{k}.mzML', out / f'fold{k}/{arm}')
            fk[f'{arm}_search_seconds'] = time.time() - T1
            W.tidy_search(out / f'fold{k}/{arm}')
            psms[arm].append(read_fold(Path(str(tsv) + '.gz'), k))
        (scr / f'fold{k}.mzML').unlink()
        res['folds'][str(k)] = fk
        json.dump(res, open(out / 'PROGRESS.json', 'w'), indent=1, default=str)
    res['databases'] = dbs
    accs = {}
    for arm in ARMS:
        comb = scr / ('MG_E.fasta' if arm == 'MG' else 'FULL_E.fasta')
        mp = out / f'{arm}_entrapment_map.tsv'
        if arm == 'FL':
            # pairing per fold database: concatenate the 10 fold maps and resolve sequences from all fold databases
            mp = out / 'FL_all_entrapment_map.tsv'
            pd.concat([pd.read_csv(out / f'FL{k}_entrapment_map.tsv', sep='\t', dtype=str) for k in range(10)]).drop_duplicates('accession').to_csv(mp, sep='\t', index=False)
            comb = scr / 'FL_all_E.fasta'
            with open(comb, 'w') as w:
                for k in range(10):
                    w.write(open(scr / f'FL{k}_E.fasta').read())
        r, a, cls, pair, table = analyse_arm(psms[arm], width, comb, mp)
        r['native_no_entrapment_pooled'] = J[f'{arm}_pooled']
        res[arm] = r; accs[arm] = (a, cls, pair)
        table.to_csv(out / f'{arm}_canonical.tsv.gz', sep='\t', index=False, compression='gzip')
    # outcomes among gained identities: accepted with FL+E but not with MG+E (same gate), FL scores and pairs
    res['gained_FL_vs_MG'] = {}
    for g in GATES:
        aFL, bFL = accs['FL'][0][g]; aMG, _ = accs['MG'][0][g]
        gained = aFL - aMG
        res['gained_FL_vs_MG'][g] = {L: estimate(aFL, bFL, accs['FL'][1], accs['FL'][2], lo, hi, restrict=gained) for L, (lo, hi) in LENS.items()}
    res['total_seconds'] = time.time() - T0; res['utc'] = C.P.now(); res['job'] = os.getenv('SLURM_JOB_ID'); res['task'] = os.getenv('SLURM_ARRAY_TASK_ID')
    shutil.rmtree(scr)
    json.dump(res, open(out / 'RESULT.json', 'w'), indent=1, default=str)
    print(json.dumps({k: v for k, v in res.items() if k in ('i', 'sample', 'total_seconds')}), flush=True)


# ------------------------------------------------------------------ collect
B = 10000; BSEED = 20261004


def collect():
    rows = []
    ro = C.roster()
    missing = []
    for r in ro:
        f = ROOT / f"per_acq/{r['cohort']}/{r['sample']}/RESULT.json"
        if not f.exists():
            missing.append(r['i']); continue
        j = json.loads(f.read_text())
        base = dict(i=int(r['i']), cohort=r['cohort'], sample=r['sample'], group=r['group'],
                    subject_id=r['subject_id'] if r['subject_id'] not in ('', 'nan') else 'NOID_' + r['sample'])
        for arm in ARMS:
            base[f'{arm}_rescoring'] = j[arm]['rescoring']
            for g in GATES:
                for L in LENS:
                    e = j[arm]['gates'][g][L]
                    for k in ['n_accepted', 'n_target', 'n_entrapment_paired', 'n_entrapment_unresolved', 'partner_not_accepted', 'entrapment_ranks_higher',
                              'exact_ties', 'paired_numerator', 'paired_estimate', 'diagnostic_with_unresolved', 'combined_estimate', 'n_psm_accepted',
                              'n_psm_entrapment', 'psm_combined_estimate']:
                        base[f'{arm}_{g}_{L}_{k}'] = e[k]
            for L in ['5', '9']:
                base[f'{arm}_native_pooled_{L}'] = j[arm]['native_no_entrapment_pooled'][f'peptides_{L}_50']
        for g in GATES:
            for L in LENS:
                e = j['gained_FL_vs_MG'][g][L]
                for k in ['n_accepted', 'n_target', 'n_entrapment_paired', 'n_entrapment_unresolved', 'paired_numerator', 'paired_estimate', 'diagnostic_with_unresolved']:
                    base[f'GAINED_{g}_{L}_{k}'] = e[k]
        rows.append(base)
    a = pd.DataFrame(rows)
    a.to_csv(ROOT / 'PER_ACQUISITION.tsv', sep='\t', index=False)
    cap = a.cohort == 'CAPSCAN'
    strata = [('PREDICT', a.cohort == 'PREDICT'), ('CAPSCAN_all', cap), ('CAPSCAN_stool', cap & (a.group == 'stool')), ('CAPSCAN_intestinal', cap & (a.group != 'stool'))]
    rng = np.random.default_rng(BSEED); out = []
    for name, m in strata:
        g = a[m].reset_index(drop=True)
        if not len(g):
            continue
        subj = g.subject_id.values; us = np.unique(subj); members = [np.flatnonzero(subj == s) for s in us]
        draws = [np.concatenate([members[k] for k in rng.integers(0, len(us), len(us))]) for _ in range(B)]
        for arm in ARMS + ['GAINED']:
            for gt in GATES:
                for L in LENS:
                    pre = f'{arm}_{gt}_{L}'
                    num_p = g[f'{pre}_paired_numerator'].astype(float); n = g[f'{pre}_n_accepted'].astype(float)
                    u = g[f'{pre}_n_entrapment_unresolved'].astype(float)
                    for est in ['paired_estimate', 'diagnostic_with_unresolved'] + ([] if arm == 'GAINED' else ['combined_estimate', 'psm_combined_estimate']):
                        v = g[f'{pre}_{est}'].astype(float).values
                        bs = np.array([np.nanmedian(v[ix]) for ix in draws])
                        rec = dict(stratum=name, arm=arm, gate=gt, length=L, estimator=est, n_acquisitions=len(g), n_participants=len(us),
                                   median_pct=100 * float(np.nanmedian(v)), ci95_low_pct=100 * float(np.nanpercentile(bs, 2.5)),
                                   ci95_high_pct=100 * float(np.nanpercentile(bs, 97.5)), n_nan=int(np.isnan(v).sum()))
                        if est == 'paired_estimate':
                            rec['pooled_sum_pct'] = 100 * float(num_p.sum() / n.sum()) if n.sum() else None
                        if est == 'diagnostic_with_unresolved':
                            rec['pooled_sum_pct'] = 100 * float((num_p + u).sum() / n.sum()) if n.sum() else None
                        out.append(rec)
            # FL minus MG paired difference of the per-acquisition paired estimate
        for gt in GATES:
            for L in LENS:
                dv = (g[f'FL_{gt}_{L}_paired_estimate'].astype(float) - g[f'MG_{gt}_{L}_paired_estimate'].astype(float)).values
                bs = np.array([np.nanmedian(dv[ix]) for ix in draws])
                out.append(dict(stratum=name, arm='FL_minus_MG', gate=gt, length=L, estimator='paired_estimate_difference', n_acquisitions=len(g),
                                n_participants=len(us), median_pct=100 * float(np.nanmedian(dv)), ci95_low_pct=100 * float(np.nanpercentile(bs, 2.5)),
                                ci95_high_pct=100 * float(np.nanpercentile(bs, 97.5)), n_nan=int(np.isnan(dv).sum())))
                gain = 100 * (g[f'FL_{gt}_{L}_n_target'].astype(float) / g[f'MG_{gt}_{L}_n_target'].replace(0, np.nan).astype(float) - 1)
                bs = np.array([np.nanmedian(gain.values[ix]) for ix in draws])
                out.append(dict(stratum=name, arm='FL_vs_MG_target_gain_with_entrapment', gate=gt, length=L, estimator='gain_pct', n_acquisitions=len(g),
                                n_participants=len(us), median_pct=float(np.nanmedian(gain)), ci95_low_pct=float(np.nanpercentile(bs, 2.5)),
                                ci95_high_pct=float(np.nanpercentile(bs, 97.5)), n_nan=int(np.isnan(gain).sum())))
    s = pd.DataFrame(out); s.to_csv(ROOT / 'SUMMARY.tsv', sep='\t', index=False)
    A92 = REVIEW / 'final/reruns_347/A8_fusion_A92/values/fusion_pooled_current.tsv'
    A91 = REVIEW / 'final/reruns_347/A9_full_route_A91/values'
    res = dict(endpoint='A-16 / D-R03', n_acquisitions=len(a), missing_rows=missing, by_cohort=a.groupby('cohort').size().to_dict(),
               lda_fallback={arm: int((a[f'{arm}_rescoring'] != 'lda').sum()) for arm in ARMS},
               headline_q001=s[(s.gate == 'q') & (s.estimator.isin(['paired_estimate', 'diagnostic_with_unresolved', 'paired_estimate_difference']))].to_dict('records'),
               fusion_A92_values_present=A92.exists(), full_route_A91_values_present=A91.exists() and any(A91.iterdir()),
               note='Entrapment appended after each database was fixed. For the cross-fitted FL and the MG arms the database does not depend on the '
                    'searched spectra, so this is an ordinary-search error control of the actual primary comparison. For the same-spectrum FULL arm '
                    'it does not capture selection from the searched spectra (see A92 linked fusion and A91 full-route controls).')
    (ROOT / 'RESULT.json').write_text(json.dumps(res, indent=1, default=str))
    lines = ['# A-16 / D-R03 result', '', f'Acquisitions collected: {len(a)} ({res["by_cohort"]}); missing: {missing}', '',
             'Peptide q <= 0.01 (pooled cross-fit rescoring), paired estimator (E + E0 + 2B)/N and with the unresolved term, median over acquisitions [participant bootstrap 95% CI], pooled sum:', '',
             '| stratum | arm | length | estimator | median % | 95% CI | pooled % |', '|---|---|---|---|---|---|---|']
    for r in res['headline_q001']:
        lines.append(f"| {r['stratum']} | {r['arm']} | {r['length']} | {r['estimator']} | {r['median_pct']:.3f} | {r['ci95_low_pct']:.3f}–{r['ci95_high_pct']:.3f} | {r.get('pooled_sum_pct', '')} |")
    lines += ['', 'PEP<=0.01 gate, observation (PSM) denominators, gained-identity strata and the FL vs MG gain under entrapment are in SUMMARY.tsv.',
              f'A92 fusion values present at collect time: {res["fusion_A92_values_present"]}; A91 full-route values present: {res["full_route_A91_values_present"]}.']
    (ROOT / 'RESULT.md').write_text('\n'.join(lines) + '\n')
    print('\n'.join(lines))


def preflight():
    """All 591 rows: native cross-fit record complete, every input path present; output dir writable;
    make_entrapment self-test on a toy FASTA (pairing, validation, ambiguous/short exclusion)."""
    ro = C.roster(); assert len(ro) == 591, len(ro)
    bad = []; n = {'PREDICT': 0, 'CAPSCAN': 0}
    for r in ro:
        try:
            t = task_for(r); nat = native_dir(r); J = json.loads((nat / 'COMPLETE.json').read_text())
            assert J['sample'] == r['sample'] and J['K'] == 10 and len(J['folds']) == 10
            for a in ARMS:
                assert J[f'{a}_pooled']['peptides_5_50'] >= 0
            for k in range(10):
                assert Path(J['folds'][str(k)]['fl_fasta']).exists()
            for p in [t['config'], t['fasta'], t['mzml'], t['csv'], r['mg_fasta']]:
                assert Path(p).exists(), p
            n[r['cohort']] += 1
        except Exception as e:  # noqa
            bad.append((r['i'], repr(e)))
    toy = ROOT / 'scratch/preflight'; toy.mkdir(parents=True, exist_ok=True)
    (toy / 't.fa').write_text('>sp|P1|A desc\nMKWVTFISLLFLFSSAYSRGVFRRDAHKSEVAHRFK\n>p2\nACDEK\n>p3\nAXCDEFGHIKLMN\n')
    st = make_entrapment(toy / 't.fa', toy / 'c.fa', 'toy', toy)
    assert st['controls_written'] == 1 and st['ambiguous_no_control'] == 1 and st['short_no_control'] == 1
    rec = dict(rows=len(ro), ok_by_cohort=n, failures=bad, toy=st, utc=C.P.now(), job=os.getenv('SLURM_JOB_ID'))
    (ROOT / 'PREFLIGHT.json').write_text(json.dumps(rec, indent=1))
    print(json.dumps(rec, indent=1))
    assert not bad and n == {'PREDICT': 296, 'CAPSCAN': 295}


if __name__ == '__main__':
    mode = sys.argv[1]
    if mode == 'collect':
        collect()
    elif mode == 'preflight':
        preflight()
    else:
        arg = sys.argv[2]
        i = int(arg) if arg.isdigit() else int(open(ROOT / arg).read().split()[int(os.environ['SLURM_ARRAY_TASK_ID'])])
        run(i, neg=(mode == 'negctrl'))
