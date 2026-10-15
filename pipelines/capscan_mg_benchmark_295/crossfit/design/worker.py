"""D11 10-fold spectrum cross-fit worker.

Modes (python worker.py <mode> <index>; index = roster row unless stated):
  cache <c*10+k>  fold-k cohort cache (c: 0 PREDICT, 1 CAPSCAN): out-of-fold prediction CSVs of ALL
                  cohort acquisitions (as the frozen parity caches used all 306/331) -> frozen extract()
                  from the original reservoir, 16 threads. Records per-acquisition CSV hashes.
  crossfit <i>    10 folds: fold-k mzML (fold-k MS2 + all MS1); selection from out-of-fold predictions
                  via the fold-k cache (frozen extractor + razor/hash-acc parsimony); Sage FL (cross-fit
                  selection), Sage FULL (archived all-prediction FASTA) and Sage MG on the fold-k spectra with the frozen per-acquisition config (only fasta/mzML/output
                  changed). Then pooled peptide q (sage_q) per arm and union of per-fold accepted sets.
  k1 <i>          no-fold control: K=1 (all spectra in one fold, selection from all predictions via the
                  production cache) -> must reproduce the archived FASTA and the frozen full-spectra
                  FL and MG accepted sets.
  negctrl <i>     cross-fit with the selection CSV replaced by the full CSV (in-fold predictions leak):
                  must fail the leakage assertion before any search.
  cachecheck <i>  fold-0 selection from the 61 GB reservoir directly vs from the fold-0 cache: equal.
  validate <i>    (a) parity routing reproduces the archived odd-scan mzML byte-for-byte and the
                  archived even/odd half CSVs byte-for-byte; (b) sage_q reproduces Sage peptide_q on
                  archived single (unsplit) Sage searches.
"""
import sys, os, json, time, gzip, shutil, traceback
from pathlib import Path
import numpy as np, pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common as C
import sage_q
import sage_lda

P = C.P
COH = ['PREDICT', 'CAPSCAN']


def gz(path):
    with open(path, 'rb') as a, gzip.open(str(path) + '.gz', 'wb', compresslevel=4) as b:
        shutil.copyfileobj(a, b)
    Path(path).unlink()


def tidy_search(d):
    (d / 'results.sage.pin').unlink(missing_ok=True)
    for n in ['results.sage.tsv', 'lfq.tsv']:
        if (d / n).exists(): gz(d / n)


def tidy_selection(d):
    st = json.loads((d / 'retrieval/selected.stats.json').read_text())
    rec = dict(retrieval_proteins=st['total_output_proteins'], retrieval_peptides=st['total_peptides'],
               retrieval_fasta_sha256=C.sha(d / 'retrieval/selected.fasta'))
    (d / 'retrieval/selected.fasta').unlink(); (d / 'retrieval/hits.tsv').unlink(missing_ok=True)
    return rec


def tol_width(config):
    tol = json.loads(Path(config).read_text())['precursor_tol']
    assert 'ppm' in tol, 'Da tolerance not implemented (all frozen configs use ppm)'
    return float(tol['ppm'][1] - tol['ppm'][0])


def pool(dfs, width):
    """PRIMARY pooled FDR: one Sage-equivalent LDA rescoring over all fold PSMs (sage_lda), then Sage's
    picked-peptide q (sage_q); frozen accepted() rule. Also returns the invalid naive pooling (raw
    per-fold scores) as a diagnostic only."""
    d = pd.concat(dfs, ignore_index=True)
    assert d.scannr.is_unique, 'a spectrum was scored in two folds'
    d['pooled_q'], _, how = sage_lda.pooled_q(d, width)
    d['naive_q'] = sage_q.peptide_q(d)
    sens = []
    for seed in (1, 2):
        d['pert_q'], _, _ = sage_lda.pooled_q(d, width, perturb_seed=seed)
        sens.append(counts(C.accepted_set(d, 'pert_q')))
    return d, C.accepted_set(d, 'pooled_q'), C.accepted_set(d, 'naive_q'), dict(rescoring=how, solver_perturbation=sens)


def fold_psms(path, k, kk):
    d = sage_lda.read(path)
    sc = d.scannr.str.extract(r'scan=(\d+)')[0].astype(int)
    assert d.label.eq(-1).any() and len(d) > 0, 'missing scores or decoys'
    bad = [s for s in sc.unique() if C.fold_of(s, kk) != k]
    assert not bad, f'spectrum from another fold entered fold {k}: {bad[:5]}'
    d['fold'] = k
    return d


def counts(acc):
    return dict(peptides_5_50=len(acc), peptides_9_50=sum(len(p) >= 9 for p in acc))


def write_set(path, acc):
    Path(path).write_text('\n'.join(sorted(acc)) + '\n')
    return C.sha(path)


def cache(index):
    co = COH[index // 10]; k = index % 10
    out = C.CF / f'caches/{co}/fold{k}'; inp = out / 'inputs'
    tasks = [t for t in C.MAN['tasks'] if t['cohort'] == co]
    assert len(tasks) == C.MAN['cohorts'][co]['n']
    shas = {}; n = {}
    t0 = time.time()
    for t in tasks:
        p = inp / f"{t['sample']}_predictions.csv"
        cnt, _ = C.write_partition_csv(t['csv'], {'oof': p}, lambda s: ['oof'] if C.fold_of(s) != k else [])
        shas[t['sample']] = C.sha(p); n[t['sample']] = cnt['oof']
    split_s = time.time() - t0
    P.extract(inp, C.MAN['cohorts'][co]['reservoir'], out, 16)
    stats = json.loads((out / 'selected.stats.json').read_text())
    rec = dict(cohort=co, fold=k, K=C.K, fold_seed=C.FOLD_SEED, reservoir=C.MAN['cohorts'][co]['reservoir'], n_acquisitions=len(tasks),
               oof_csv_sha256=shas, oof_rows=n, split_seconds=split_s, stats=stats, fasta_sha256=C.sha(out / 'selected.fasta'),
               hits_bytes=(out / 'hits.tsv').stat().st_size, utc=P.now())
    shutil.rmtree(inp); (out / 'hits.tsv').unlink()
    C.dump(out / 'COMPLETE.json', rec)
    print(json.dumps({k2: v for k2, v in rec.items() if k2 not in ('oof_csv_sha256', 'oof_rows')}), flush=True)


def crossfit(i, leak=False, kk=C.K, tag='crossfit'):
    row = C.roster()[i]; t = C.task_for(row)
    assert C.sha(t['config']) == t['config_sha256'] and C.sha(t['fasta']) == t['fasta_sha256']
    out = C.CF / f"outputs/{row['cohort']}/{row['sample']}"
    if tag != 'crossfit': out = C.CF / f"{tag}_outputs/{row['cohort']}/{row['sample']}"
    scr = C.CF / f"scratch/{tag}/{row['cohort']}/{row['sample']}"
    shutil.rmtree(scr, ignore_errors=True); scr.mkdir(parents=True); out.mkdir(parents=True, exist_ok=True)
    res = dict(mode=tag, i=i, cohort=row['cohort'], sample=row['sample'], subject_id=row['subject_id'], group=row['group'], K=kk,
               fold_seed=C.FOLD_SEED, sage_sha256=C.sha(P.SAGE), pipeline_sha256=C.PIPE_SHA, mg_fasta=row['mg_fasta'],
               mg_fasta_sha256=C.sha(row['mg_fasta']), config=t['config'], folds={})
    T = time.time()
    rec, fold_scans = C.fold_mzml(t['mzml'], t['csv'], {k: scr / f'fold{k}.mzML' for k in range(kk)}, lambda s: C.fold_of(s, kk), out / 'FOLD_MZML.json')
    res['fold_mzml_seconds'] = time.time() - T
    oof = {k: scr / f'oof{k}' / f"{row['sample']}_predictions.csv" for k in range(kk)}
    cnt, all_pred = C.write_partition_csv(t['csv'], oof, lambda s: [k for k in range(kk) if C.fold_of(s, kk) != k])
    ARMS = ['FL', 'FULL', 'MG']
    psms = {a: [] for a in ARMS}; union = {a: set() for a in ARMS}
    for k in range(kk):
        fk = {}; T0 = time.time()
        sel_csv = oof[k] if not leak else Path(t['csv'])
        # leakage guard: no prediction of a fold-k spectrum may enter the fold-k selection
        sel_scans = C.csv_scans(sel_csv)
        leaked = sel_scans & set(fold_scans[k])
        assert not leaked, f'LEAKAGE: {len(leaked)} fold-{k} spectra have predictions in the fold-{k} selection input'
        assert all(C.fold_of(s, kk) != k for s in sel_scans)
        cache_rec = json.loads((C.CF / f"caches/{row['cohort']}/fold{k}/COMPLETE.json").read_text())
        fk['oof_csv_sha256'] = C.sha(oof[k]); fk['oof_rows'] = cnt[k]
        assert fk['oof_csv_sha256'] == cache_rec['oof_csv_sha256'][row['sample']], 'out-of-fold CSV differs from the one in the fold cache'
        fd = out / f'fold{k}'
        fasta = C.select(row['sample'], oof[k], C.CF / f"caches/{row['cohort']}/fold{k}/selected.fasta", fd / 'selection')
        fk.update(tidy_selection(fd / 'selection'))
        fk['fl_fasta'] = str(fasta); fk['fl_fasta_sha256'] = C.sha(fasta); fk['fl_fasta_sequences'] = len(P.seqset(fasta))
        fk['selection_seconds'] = time.time() - T0
        fk['fold_ms2'] = rec['outputs'][str(k)]['ms2']; fk['fold_mzml_sha256'] = rec['outputs'][str(k)]['sha256']
        for arm, db in [('FL', fasta), ('FULL', t['fasta']), ('MG', row['mg_fasta'])]:
            T1 = time.time()
            tsv = C.sage_search(t['config'], db, scr / f'fold{k}.mzML', fd / arm)
            fk[f'{arm}_search_seconds'] = time.time() - T1
            d = fold_psms(tsv, k, kk)
            acc = C.accepted_set(d)
            fk[arm] = dict(psms=len(d), decoy_psms=int((d.label == -1).sum()), **counts(acc))
            union[arm] |= acc; psms[arm].append(d)
            write_set(fd / arm / 'ACCEPTED.txt', acc)
            tidy_search(fd / arm)
        (scr / f'fold{k}.mzML').unlink()
        res['folds'][str(k)] = fk
        C.dump(out / 'PROGRESS.json', res)
    for arm in ARMS:
        d, acc, naive, how = pool(psms[arm], tol_width(t['config']))
        res[f'{arm}_pooled'] = dict(psms=len(d), decoy_psms=int((d.label == -1).sum()), **how, **counts(acc),
                                    accepted_sha256=write_set(out / f'{arm}_POOLED_ACCEPTED.txt', acc))
        res[f'{arm}_pooled_naive_diagnostic'] = counts(naive)
        res[f'{arm}_union'] = dict(**counts(union[arm]), accepted_sha256=write_set(out / f'{arm}_UNION_ACCEPTED.txt', union[arm]))
    res['total_seconds'] = time.time() - T; res['utc'] = P.now()
    shutil.rmtree(scr)
    C.dump(out / 'COMPLETE.json', res)
    print(json.dumps({k: v for k, v in res.items() if k != 'folds'}), flush=True)


def repool(i):
    """Recompute the pooled sets of a finished cross-fit acquisition from its saved fold tables
    (used for the pilot, which was first pooled with the invalid raw-score method)."""
    row = C.roster()[i]; t = C.task_for(row)
    out = C.CF / f"outputs/{row['cohort']}/{row['sample']}"
    res = json.loads((out / 'COMPLETE.json').read_text()); kk = res['K']
    for arm in ['FL', 'FULL', 'MG']:
        d, acc, naive, how = pool([fold_psms(out / f'fold{k}/{arm}/results.sage.tsv.gz', k, kk) for k in range(kk)], tol_width(t['config']))
        res[f'{arm}_pooled'] = dict(psms=len(d), decoy_psms=int((d.label == -1).sum()), **how, **counts(acc),
                                    accepted_sha256=write_set(out / f'{arm}_POOLED_ACCEPTED.txt', acc))
        res[f'{arm}_pooled_naive_diagnostic'] = counts(naive)
    res['repooled_utc'] = P.now()
    C.dump(out / 'COMPLETE.json', res); print(json.dumps({k: v for k, v in res.items() if k != 'folds'}), flush=True)


def k1pool(i):
    """K=1 through the pooled-FDR path: the single all-spectra search of each arm, rescored and q-valued
    by sage_lda + sage_q, must give the frozen full-spectra sets."""
    row = C.roster()[i]; t = C.task_for(row)
    out = C.CF / f"k1_outputs/{row['cohort']}/{row['sample']}"
    res = json.loads((out / 'COMPLETE.json').read_text())
    man = json.loads(C.CONF_MAN.read_text())
    o = next(x for x in man['ordinary'] if x['cohort'] == row['cohort'] and x['sample'] == row['sample'])
    for arm, key in [('FL', 'FL_public'), ('MG', 'MG')]:
        d = sage_lda.read(out / f'{arm}/results.sage.tsv.gz')
        d2, acc, _, how = pool([d], tol_width(t['config']))
        r = dict(**how, max_abs_score_diff=float(np.abs(sage_lda.rescore(d, tol_width(t['config']))[0] - d.sage_discriminant_score.to_numpy(np.float32)).max()))
        for L in ['5', '9']:
            exp = set(Path(o['expected_sets'][key][L]['path']).read_text().split())
            mine = {p for p in acc if len(p) >= int(L)}
            r[f'equals_frozen_{L}_50'] = mine == exp; r[f'n_{L}'] = len(mine); r[f'frozen_{L}'] = len(exp)
        res[f'{arm}_pooled_path'] = r
    C.dump(out / 'COMPLETE.json', res); print(json.dumps(res), flush=True)
    assert all(res[f'{a}_pooled_path'][f'equals_frozen_{L}_50'] for a in ['FL', 'MG'] for L in ['5', '9'])


def k1(i):
    row = C.roster()[i]; t = C.task_for(row)
    assert C.sha(t['config']) == t['config_sha256'] and C.sha(t['fasta']) == t['fasta_sha256']
    out = C.CF / f"k1_outputs/{row['cohort']}/{row['sample']}"; scr = C.CF / f"scratch/k1/{row['cohort']}/{row['sample']}"
    shutil.rmtree(scr, ignore_errors=True); scr.mkdir(parents=True); out.mkdir(parents=True, exist_ok=True)
    res = dict(mode='k1', i=i, cohort=row['cohort'], sample=row['sample'])
    T = time.time()
    rec, _ = C.fold_mzml(t['mzml'], t['csv'], {0: scr / 'all.mzML'}, lambda s: 0, out / 'FOLD_MZML.json')
    res['all_spectra_rewritten'] = rec['outputs']['0']['ms2'] == rec['original_ms2']
    fasta = C.select(row['sample'], Path(t['csv']), t['production_cache'], out / 'selection')
    res.update(tidy_selection(out / 'selection'))
    a = P.seqset(fasta); b = P.seqset(t['fasta'])
    res['selection_equals_archived_fasta'] = a == b; res['selected'] = len(a); res['archived'] = len(b)
    man = json.loads(C.CONF_MAN.read_text())
    o = next(x for x in man['ordinary'] if x['cohort'] == row['cohort'] and x['sample'] == row['sample'])
    for arm, db, key in [('FL', fasta, 'FL_public'), ('MG', row['mg_fasta'], 'MG')]:
        T1 = time.time()
        tsv = C.sage_search(t['config'], db, scr / 'all.mzML', out / arm)
        res[f'{arm}_search_seconds'] = time.time() - T1
        d = sage_q.read_psms(tsv)
        d['re_q'] = sage_q.peptide_q(d)
        acc = C.accepted_set(d); acc_re = C.accepted_set(d, 're_q')
        rel = np.abs(d.re_q - d.peptide_q) / np.maximum(d.peptide_q, 1e-12)
        res[arm] = dict(**counts(acc), reimpl_equal_sets=acc == acc_re, reimpl_threshold_disagreements=int(((d.re_q <= .01) != (d.peptide_q <= .01)).sum()),
                        reimpl_max_rel_diff=float(rel.max()))
        for L in ['5', '9']:
            exp = set(Path(o['expected_sets'][key][L]['path']).read_text().split())
            mine = {p for p in acc if len(p) >= int(L)}
            res[arm][f'equals_frozen_{L}_50'] = mine == exp; res[arm][f'frozen_{L}_50'] = len(exp)
            res[arm][f'only_k1_{L}'] = len(mine - exp); res[arm][f'only_frozen_{L}'] = len(exp - mine)
        write_set(out / arm / 'ACCEPTED.txt', acc); tidy_search(out / arm)
    res['total_seconds'] = time.time() - T; res['utc'] = P.now()
    shutil.rmtree(scr)
    C.dump(out / 'COMPLETE.json', res); print(json.dumps(res), flush=True)


def cachecheck(i):
    row = C.roster()[i]; t = C.task_for(row)
    out = C.CF / f"cachecheck_outputs/{row['cohort']}/{row['sample']}"; out.mkdir(parents=True, exist_ok=True)
    csv0 = out / f"oof0/{row['sample']}_predictions.csv"
    C.write_partition_csv(t['csv'], {0: csv0}, lambda s: [0] if C.fold_of(s) != 0 else [])
    cache_rec = json.loads((C.CF / f"caches/{row['cohort']}/fold0/COMPLETE.json").read_text())
    assert C.sha(csv0) == cache_rec['oof_csv_sha256'][row['sample']]
    res = dict(mode='cachecheck', i=i, cohort=row['cohort'], sample=row['sample'])
    fa = C.select(row['sample'], csv0, C.CF / f"caches/{row['cohort']}/fold0/selected.fasta", out / 'via_cache')
    fb = C.select(row['sample'], csv0, t['reservoir'], out / 'via_reservoir')
    ra = P.seqset(out / 'via_cache/retrieval/selected.fasta'); rb = P.seqset(out / 'via_reservoir/retrieval/selected.fasta')
    a = P.seqset(fa); b = P.seqset(fb)
    res.update(retrieval_equal=ra == rb, retrieval_n=(len(ra), len(rb)), selection_equal=a == b, selection_n=(len(a), len(b)),
               selection_sha=(C.sha(fa), C.sha(fb)), utc=P.now())
    for d in ['via_cache', 'via_reservoir']:
        tidy_selection(out / d)
    shutil.rmtree(out / 'oof0')
    C.dump(out / 'COMPLETE.json', res); print(json.dumps(res), flush=True)
    assert res['retrieval_equal'] and res['selection_equal']


def validate(i):
    row = C.roster()[i]; t = C.task_for(row)
    out = C.CF / f"validation/{row['cohort']}/{row['sample']}"; shutil.rmtree(out, ignore_errors=True); out.mkdir(parents=True)
    res = dict(mode='validate', i=i, cohort=row['cohort'], sample=row['sample'])
    # (a1) parity routing through the generalised writer reproduces the archived odd-scan mzML
    rec, _ = C.fold_mzml(t['mzml'], t['csv'], {1: out / 'odd.mzML'}, lambda s: s % 2, out / 'ODD.json')
    res['odd_mzml_sha256'] = rec['outputs']['1']['sha256']; res['odd_mzml_identical'] = res['odd_mzml_sha256'] == row['odd_mzml_sha256']
    (out / 'odd.mzML').unlink()
    # (a2) parity routing through the CSV writer reproduces the archived half CSVs
    halves = {p: out / f'{p}.csv' for p in ['even', 'odd']}
    C.write_partition_csv(t['csv'], halves, lambda s: ['even' if s % 2 == 0 else 'odd'])
    inp = json.loads((C.SH / f"inputs/{row['cohort']}/COMPLETE.json").read_text())
    exp = next(x for x in inp['samples'] if x['sample'] == row['sample'])['half_sha256']
    res['half_csv_identical'] = {p: C.sha(f) == exp[p] for p, f in halves.items()}
    for f in halves.values(): f.unlink()
    # (b) Sage peptide q re-implementation on unsplit single Sage searches
    man = json.loads(C.CONF_MAN.read_text())
    o = next(x for x in man['ordinary'] if x['cohort'] == row['cohort'] and x['sample'] == row['sample'])
    srcs = {f'split_{a}': C.REC / f"outputs/{row['cohort']}/{row['sample']}/{a}/results.sage.tsv" for a in ['FULL', 'HALF_SELF', 'HALF_DISJOINT']}
    srcs['frozen_full_FL'] = Path(o['searches']['FL_public']['path']); srcs['frozen_full_MG'] = Path(o['searches']['MG']['path'])
    srcs['round3_MG_ODD'] = C.REVIEW / f"iteration3/disjoint/outputs/{row['cohort']}/{row['sample']}/MG_ODD/results.sage.tsv.gz"
    res['q'] = {}
    for name, p in srcs.items():
        if not p.exists():
            res['q'][name] = 'missing'; continue
        d = sage_lda.read(p)
        T = time.time(); q = sage_q.peptide_q(d); dt = time.time() - T
        sc, how = sage_lda.rescore(d, tol_width(t['config']))
        e = d[['peptide', 'label']].copy(); e['sage_discriminant_score'] = sc
        d['lda_q'] = sage_q.peptide_q(e)
        s = d.peptide_q.to_numpy(np.float32)
        rel = np.abs(q - s) / np.maximum(s, 1e-12)
        d['re_q'] = q
        a, b = C.accepted_set(d), C.accepted_set(d, 're_q')
        res['q'][name] = dict(psms=len(d), max_abs_diff=float(np.abs(q - s).max()), max_rel_diff=float(rel.max()), n_rel_diff_gt_1e5=int((rel > 1e-5).sum()),
                              threshold_disagreements=int(((q <= .01) != (s <= .01)).sum()), accepted_equal=a == b, accepted_sage=len(a), accepted_reimpl=len(b), seconds=dt,
                              lda_rescoring=how, lda_max_abs_score_diff=float(np.abs(sc - d.sage_discriminant_score.to_numpy(np.float32)).max()),
                              lda_q_threshold_disagreements=int(((d.lda_q <= .01) != (s <= .01)).sum()), lda_accepted_equal=C.accepted_set(d, 'lda_q') == a)
    res['utc'] = P.now()
    C.dump(out / 'COMPLETE.json', res); print(json.dumps(res), flush=True)
    assert res['odd_mzml_identical'] and all(res['half_csv_identical'].values())
    assert all(v == 'missing' or (v['accepted_equal'] and v['threshold_disagreements'] == 0 and v['lda_accepted_equal']) for v in res['q'].values())


if __name__ == '__main__':
    mode = sys.argv[1]; ids = sys.argv[2]
    idx = int(ids) if ids.isdigit() else int(open(C.DES / ids).read().split()[int(os.environ['SLURM_ARRAY_TASK_ID'])])
    {'cache': cache, 'crossfit': crossfit, 'k1': k1, 'cachecheck': cachecheck, 'validate': validate, 'repool': repool, 'k1pool': k1pool,
     'negctrl': lambda j: crossfit(j, leak=True, tag='negctrl')}[mode](idx)
