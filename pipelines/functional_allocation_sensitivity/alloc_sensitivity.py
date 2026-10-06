"""GH-029: functional-allocation sensitivity on frozen annotated-origin consensus outputs.

Read-only on every upstream file. Writes only below this script's directory.

Variants (per peptide record, per annotation level):
  equal            production rule (profiles.py): intensity split equally over distinct eligible terms
  unique           only records with exactly one distinct eligible term are allocated; multi-term -> unassigned
  complete         only records whose every compatible parent carries an annotation at this level
                   (function['missing'] is False) are allocated, equal split; partial -> unassigned
  unique_complete  both restrictions
Human-assigned withholding is identical in every variant (annotated-origin consensus, profiles.human).

Modes:
  pilot      CAMPI_REF only; asserts that the 'equal' variant reproduces the frozen profile matrices exactly
  profiles   HSTOOL_noclust, PREDICT_noclust, CAPSCAN_noclust; same equality assertion for 'equal'
  endpoints  HSTOOL KO pair precision (Fig. 3f), CapScan own-stool OG similarity (Fig. 5c),
             PREDICT functional cross-omics (Supp. Fig. 11 / Fig. 4i functional rows), each per variant,
             plus negative controls (label permutation).
"""
from pathlib import Path
import ast, gzip, hashlib, json, math, os, re, sys
from itertools import combinations
import numpy as np
import pandas as pd
from scipy import sparse

HERE = Path(__file__).resolve().parent
C = next(p for p in HERE.parents if p.name == 'manuscript_completion_2026-10-01')
B = C / 'biology_refit'
J = C.parents[1]
R = next(p for p in HERE.parents if p.name == 'FastaLake')
A = C.parent / 'github_workflow_reconciliation_2026-09-30/resumed_audit/annotation_repair'
CONS = A / 'annotated_origin_consensus_2026-10-01/full'
FROZEN = C / 'profiles'
OUT = HERE / 'out'
LEVELS = ['KO', 'OG', 'Pathway']
VARIANTS = ['equal', 'unique', 'complete', 'unique_complete']
UNASSIGNED = 'Unassigned_at_this_resolution'
HUMAN = 'Human_assigned_signal_withheld'


def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''):
            h.update(b)
    return h.hexdigest()


def load(p):
    return json.loads(Path(p).read_text())


def dump(p, x):
    Path(p).write_text(json.dumps(x, indent=2, default=str) + '\n')


def human(t):  # identical to profiles.py
    return t['rank'] == 'species' and t['supported_lineage'][-1:] == ['s__Homo sapiens']


def eligible(f, level):  # identical to profiles.py
    if level == 'KO':
        return sorted({t for t in f['terms'] if re.fullmatch(r'K\d{5}', t)})
    if level == 'Pathway':
        return sorted({t for t in f['terms'] if re.fullmatch(r'map\d{5}', t) and not t.startswith('map011')})
    return sorted({s[0] for s in f['annotation_sets'] if len(s) == 1 and re.fullmatch(r'[^@]+@(?:2\|Bacteria|2157\|Archaea)', s[0])})


def variant_terms(f, level, variant):
    ts = eligible(f, level)
    if variant in ('unique', 'unique_complete') and len(ts) != 1:
        ts = []
    if variant in ('complete', 'unique_complete') and f['missing']:
        ts = []
    return ts


def projector(recs, level, exclude_human, variant):
    assignments = []
    for r in recs:
        ts = variant_terms(r['function'], level, variant)
        assignments.append([HUMAN] if exclude_human and human(r['taxonomy']) else (ts or [UNASSIGNED]))
    labels = sorted({t for ts in assignments for t in ts} - {UNASSIGNED, HUMAN}) + [UNASSIGNED] + ([HUMAN] if exclude_human else [])
    ix = {t: i for i, t in enumerate(labels)}
    rr, cc, vv = [], [], []
    for j, ts in enumerate(assignments):
        for t in ts:
            rr.append(ix[t]); cc.append(j); vv.append(1 / len(ts))
    p = sparse.coo_matrix((vv, (rr, cc)), shape=(len(labels), len(recs))).tocsr()
    assert np.allclose(p.sum(0).A.ravel(), 1, rtol=0, atol=1e-12)
    return labels, p


def unit_tests():
    tax = dict(rank='unresolved', supported_lineage=[])
    def r(ts, missing=False):
        return dict(function=dict(terms=ts, missing=missing, annotation_sets=[[t] for t in ts]), taxonomy=tax)
    recs = [r(['K00001', 'K00002']), r(['K00001'], True), r(['K00003']), r([])]
    x = np.array([4., 6., 5., 7.])
    expect = {
        'equal': {'K00001': 8., 'K00002': 2., 'K00003': 5., UNASSIGNED: 7.},
        'unique': {'K00001': 6., 'K00003': 5., UNASSIGNED: 11.},
        'complete': {'K00001': 2., 'K00002': 2., 'K00003': 5., UNASSIGNED: 13.},
        'unique_complete': {'K00003': 5., UNASSIGNED: 17.},
    }
    for v, e in expect.items():
        names, p = projector(recs, 'KO', False, v)
        got = dict(zip(names, p @ x))
        assert got == e, (v, got)
        assert math.isclose(sum(got.values()), x.sum())
    recs[2]['taxonomy'] = dict(rank='species', supported_lineage=['s__Homo sapiens'])
    names, p = projector(recs, 'KO', True, 'unique')
    assert dict(zip(names, p @ x)) == {'K00001': 6., UNASSIGNED: 11., HUMAN: 5.}


def build_arm(arm):
    d = CONS / arm
    src = load(d / 'COMPLETE.json'); assert src['status'] == 'PASS'
    used = ['sample_totals.npz', 'quantity_definitions.jsonl.gz'] + [f'{l}_{s}' for l in LEVELS for s in ['quantity_ids.json', 'evidence_signal.npz']]
    for name in used:
        assert sha(d / name) == src['outputs'][name]['sha256'], name
    z = np.load(d / 'sample_totals.npz'); total = z['intensity']; samples = list(z['samples'])
    defs = {}
    with gzip.open(d / 'quantity_definitions.jsonl.gz', 'rt') as f:
        for r in map(json.loads, f):
            defs[r['quantity_id']] = r
    dest = OUT / 'profiles' / arm; dest.mkdir(parents=True, exist_ok=True)
    cov = []; checks = {}
    for lev in LEVELS:
        ids = load(d / f'{lev}_quantity_ids.json'); recs = [defs[k] for k in ids]
        mat = sparse.load_npz(d / f'{lev}_evidence_signal.npz')
        assert np.allclose(mat.sum(0).A.ravel(), total, rtol=1e-9, atol=1e-5)
        for exclude, pname in [(False, 'all'), (True, 'other_unresolved')]:
            for v in VARIANTS:
                labels, p = projector(recs, lev, exclude, v)
                res = (p @ mat).tocsr()
                assert np.allclose(res.sum(0).A.ravel(), total, rtol=1e-9, atol=1e-5)
                mask = np.array([t not in (UNASSIGNED, HUMAN) for t in labels])
                named = res[mask]; feats = list(np.array(labels)[mask])
                un = res[labels.index(UNASSIGNED)].toarray().ravel()
                hw = res[labels.index(HUMAN)].toarray().ravel() if exclude else np.zeros_like(total)
                stem = f'{lev}_{pname}__{v}'
                sparse.save_npz(dest / (stem + '.npz'), named.T.tocsr())
                dump(dest / (stem + '_axes.json'), dict(samples=samples, features=feats, variant=v, human_exclusion=exclude))
                alloc = named.sum(0).A.ravel()
                for i, s in enumerate(samples):
                    cov.append(dict(arm=arm, level=lev, profile=pname, variant=v, sample=s, total_LFQ=float(total[i]),
                                    allocated=float(alloc[i]), unassigned=float(un[i]), human_withheld=float(hw[i])))
                if v == 'equal':
                    fa = load(FROZEN / arm / f'{lev}_{pname}_axes.json')
                    fm = sparse.load_npz(FROZEN / arm / f'{lev}_{pname}.npz')
                    assert fa['samples'] == samples and fa['features'] == feats, ('axes differ', arm, stem)
                    diff = abs(fm - named.T.tocsr())
                    mx = float(diff.max()) if diff.nnz else 0.0
                    assert mx <= max(1e-6, float(total.max()) * 1e-12), ('equal variant does not reproduce frozen', arm, stem, mx)
                    checks[stem] = dict(reproduces_frozen=True, max_abs_diff=mx, frozen_sha256=sha(FROZEN / arm / f'{lev}_{pname}.npz'))
        print('BUILT', arm, lev, flush=True)
    pd.DataFrame(cov).to_csv(dest / 'coverage.tsv', sep='\t', index=False)
    dump(dest / 'COMPLETE.json', dict(status='PASS', arm=arm, positive_control=checks, source_receipt_sha256=sha(d / 'COMPLETE.json'),
                                      code_sha256=sha(__file__), job=os.environ.get('SLURM_JOB_ID')))


def frame(arm, key):
    a = load(OUT / 'profiles' / arm / (key + '_axes.json')); v = sparse.load_npz(OUT / 'profiles' / arm / (key + '.npz'))
    return pd.DataFrame(v.toarray(), index=a['samples'], columns=a['features'])


def import_functions(path, names, env):
    tree = ast.parse(Path(path).read_text())
    defs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    assert {n.name for n in defs} == set(names), (path, names)
    exec(compile(ast.Module(body=defs, type_ignores=[]), str(path), 'exec'), env)
    return env


def hstool_precision():
    hp = import_functions(C / 'hstool_precision.py', ['measure', 'pairs', 'summarize'], dict(np=np, pd=pd, METRICS=['r_per_pair', 'cv_per_pair_pct', 'n_features_per_pair']))
    design = pd.read_csv(C / 'hstool_precision_refit/design.tsv', sep='\t', dtype=str)
    tiers = ['TechnicalReplicates', 'WorkflowReplicates', 'OnePiece', 'Intra', 'Interindividual']
    groups = {t: design.loc[design.group.eq(t), 'sample_id'].tolist() for t in tiers}
    pairlist = [(t, a, b) for t in tiers for a, b in combinations(groups[t], 2)]; assert len(pairlist) == 316
    samples = design.sample_id.tolist(); rows = []
    for pname in ['all', 'other_unresolved']:
        for v in VARIANTS:
            f = frame('HSTOOL_noclust', f'KO_{pname}__{v}').T.reindex(columns=samples)
            s = hp['summarize'](hp['pairs'](f, pairlist), ['design']); s['profile'] = pname; s['variant'] = v; rows.append(s)
    out = pd.concat(rows, ignore_index=True)
    ref = pd.read_csv(C / 'hstool_precision_refit/native_feature_summary.tsv', sep='\t')
    for pname in ['all', 'other_unresolved']:
        a = out[(out.profile == pname) & (out.variant == 'equal')].set_index('design')
        b = ref[ref.profile == 'new_' + pname].set_index('design')
        for m in ['r_per_pair', 'cv_per_pair_pct', 'n_features_per_pair']:
            assert np.allclose(a.loc[b.index, m], b[m], rtol=1e-9, atol=1e-9), ('HSTOOL positive control', pname, m)
    # Negative control: shuffle KO labels independently per acquisition (destroys feature identity).
    rng = np.random.default_rng(20261003)
    f = frame('HSTOOL_noclust', 'KO_other_unresolved__equal').T.reindex(columns=samples)
    sh = f.apply(lambda c: pd.Series(rng.permutation(c.to_numpy()), index=c.index))
    neg = hp['summarize'](hp['pairs'](sh, pairlist), ['design']); neg['profile'] = 'other_unresolved'; neg['variant'] = 'NEGATIVE_CONTROL_feature_shuffle'
    out = pd.concat([out, neg], ignore_index=True)
    out.to_csv(OUT / 'hstool_KO_precision.tsv', sep='\t', index=False)
    return out


def capscan_own_stool():
    from scipy.spatial.distance import cdist
    SD = R / 'MM_review_analysis/A33_4_6_peptide_centric_downstream_2026-09-19/results/CAPSCAN__simplified/source_data'
    REG = ['duodenum', 'jejunum', 'ileum', 'colon', 'stool']
    design = pd.read_csv(SD / 'Figure5/design.tsv', sep='\t').set_index('sample_id')
    ix = pd.MultiIndex.from_product([range(1, 15), REG], names=['participant', 'region'])
    old = pd.read_csv(SD / 'CAPSCAN_noclust/microbial_OG.tsv', sep='\t', index_col=0).loc[design.index]
    roster = old.index[old.sum(1) > 0]; assert len(roster) == 330
    perms = np.argsort(np.random.default_rng(20260910).random((99999, 14)), axis=1)
    def cells(raw):
        d = design.loc[roster]; v = raw.loc[roster]; ok = v.sum(1) > 0
        rel = v.div(v.sum(1).where(ok), axis=0)
        spec = rel.groupby([d.participant, d.region, d.specimen_id]).mean()
        return spec.groupby(level=[0, 1]).mean().loc[ix], int((~ok).sum())
    def stats(x, prev):
        keep = (x > 0).mean(0) >= prev; v = x[:, keep]; h = np.sqrt(v / v.sum(1, keepdims=True)); cube = h.reshape(14, 5, -1)
        dd = np.stack([cdist(cube[i, :4], cube[:, 4]).mean(0) for i in range(14)]); own = np.diag(dd); other = (dd.sum(1) - own) / 13
        obs = float((other - own).mean()); pdg = dd[np.arange(14)[None, :], perms].mean(1); null = (dd.sum(1).mean() - 14 * pdg) / 13
        p = (1 + int((abs(null) >= abs(obs) - 1e-14).sum())) / 100000
        return dict(prevalence=prev, features=int(keep.sum()), own_closer=int((own < other).sum()), mean_own=float(own.mean()),
                    mean_other=float(other.mean()), mean_advantage=obs, p=p)
    rows = []
    for v in VARIANTS:
        cc, zero = cells(frame('CAPSCAN_noclust', f'OG_other_unresolved__{v}').reindex(design.index))
        x = cc.to_numpy()
        for prev in [0., .1, .2]:
            rows.append(dict(variant=v, acquisitions_without_allocated_signal=zero, **stats(x, prev)))
    # Negative control: stool profiles assigned to the wrong participants (fixed derangement).
    cc, _ = cells(frame('CAPSCAN_noclust', 'OG_other_unresolved__equal').reindex(design.index))
    cube = cc.to_numpy().reshape(14, 5, -1).copy(); cube[:, 4] = cube[np.roll(np.arange(14), 1), 4]
    rows.append(dict(variant='NEGATIVE_CONTROL_stool_rotated', acquisitions_without_allocated_signal=0, **stats(cube.reshape(70, -1), 0.)))
    out = pd.DataFrame(rows)
    ref = pd.read_csv(B / 'capscan_refit/similarity_summary.tsv', sep='\t')
    e = out[out.variant == 'equal'].set_index('prevalence'); r = ref.set_index('prevalence')
    for m in ['features', 'own_closer', 'mean_own', 'mean_other', 'mean_advantage', 'p']:
        assert np.allclose(e.loc[r.index, m], r[m], rtol=1e-9, atol=1e-12), ('CapScan positive control', m)
    out.to_csv(OUT / 'capscan_own_stool.tsv', sep='\t', index=False)
    return out


def predict_crossomics():
    from scipy.spatial.distance import pdist, squareform
    from scipy.stats import rankdata
    X = J / 'source_data/current/Crossomics'
    OLDE = R / 'analyses/enterotypes_2026-09-07'
    env = dict(np=np, pd=pd, pdist=pdist, squareform=squareform, rankdata=rankdata, combinations=combinations,
               SEED=20260907, NPERM=1999, results=[], coverage=[], rosters=[], matrices={})
    import_functions(X / 'provenance/check_prevalence.py', ['relative', 'logmet', 'transform', 'evaluate'], env)
    def rd(p, **kw): return pd.read_csv(p, sep='\t', **kw)
    def bframe(folder, key):
        a = load(folder / (key + '_axes.json')); v = sparse.load_npz(folder / (key + '.npz'))
        return pd.DataFrame(v.toarray(), index=a['samples'], columns=a['features'])
    def collapse(d):
        ids = d.index.to_series().str.extract(r'(MicrobPred_\d+)', expand=False); assert ids.notna().all()
        return d.groupby(ids).mean()
    first = rd(OLDE / 'PREDICT_fullstack_primary_samples.tsv', index_col=0); assert len(first) == 77
    host = collapse(bframe(B / 'groups/PREDICT', 'supported_human_groups'))
    met = env['logmet'](rd(OLDE / 'PREDICT_metabolite_intensity.tsv', index_col=0))
    rows = []
    variants = VARIANTS + ['NEGATIVE_CONTROL_protein_patients_permuted']
    for v in variants:
        for lev in LEVELS:
            vv = 'equal' if v.startswith('NEGATIVE') else v
            prot = collapse(frame('PREDICT_noclust', f'{lev}_other_unresolved__{vv}'))
            if v.startswith('NEGATIVE'):
                rng = np.random.default_rng(20261003); sub = prot.loc[first.index]
                prot = pd.DataFrame(sub.to_numpy()[rng.permutation(len(sub))], index=sub.index, columns=sub.columns)
            fs = dict(metaG=rd(OLDE / ('metaG_OG_TPM.tsv' if lev == 'OG' else 'PREDICT_metaG_' + lev + '.tsv'), index_col=0),
                      metaT=rd(OLDE / ('metaT_OG_TPM.tsv' if lev == 'OG' else 'PREDICT_metaT_' + lev + '.tsv'), index_col=0),
                      microbial_metaP=prot, host_metaP=host, metabolites=met)
            before = len(env['results']); cbefore = len(env['coverage'])
            env['evaluate']('PREDICT', 'displayed77', lev, fs, first.index, 0.)
            for r in env['results'][before:]:
                rows.append(dict(variant=v, **r))
            for r in env['coverage'][cbefore:]:
                r['variant'] = v
        print('CROSS', v, flush=True)
    out = pd.DataFrame(rows)
    out = out[(out.layer_a == 'microbial_metaP') | (out.layer_b == 'microbial_metaP')].copy()
    ref = rd(B / 'cross_refit/functional_correlations.tsv')
    ref = ref[(ref.cohort == 'PREDICT') & (ref.prevalence == 0) & ((ref.layer_a == 'microbial_metaP') | (ref.layer_b == 'microbial_metaP'))]
    m = out[out.variant == 'equal'].merge(ref, on=['level', 'layer_a', 'layer_b'], suffixes=('', '_ref'))
    assert len(m) == len(ref) == 12 and np.allclose(m.rho, m.rho_ref, atol=1e-12) and np.allclose(m.p, m.p_ref, atol=1e-12), 'PREDICT positive control'
    cov = pd.DataFrame(env['coverage']); cov = cov[cov.layer == 'microbial_metaP']
    out.to_csv(OUT / 'predict_crossomics_protein_rows.tsv', sep='\t', index=False)
    cov.to_csv(OUT / 'predict_crossomics_protein_coverage.tsv', sep='\t', index=False)
    return out


def coverage_summary():
    rows = []
    for arm in ['HSTOOL_noclust', 'PREDICT_noclust', 'CAPSCAN_noclust']:
        c = pd.read_csv(OUT / 'profiles' / arm / 'coverage.tsv', sep='\t')
        c = c[c.profile == 'other_unresolved']
        c['allocated_fraction'] = c.allocated / c.total_LFQ.where(c.total_LFQ > 0)
        rows.append(c.groupby(['arm', 'level', 'variant']).allocated_fraction.median().reset_index())
    out = pd.concat(rows); out.to_csv(OUT / 'allocated_fraction_medians.tsv', sep='\t', index=False); return out


def main(mode):
    unit_tests()
    OUT.mkdir(exist_ok=True)
    if mode == 'pilot':
        build_arm('CAMPI_REF')
    elif mode == 'profiles':
        for arm in ['HSTOOL_noclust', 'CAPSCAN_noclust', 'PREDICT_noclust']:
            build_arm(arm)
    elif mode == 'endpoints':
        for arm in ['HSTOOL_noclust', 'CAPSCAN_noclust', 'PREDICT_noclust']:
            assert load(OUT / 'profiles' / arm / 'COMPLETE.json')['status'] == 'PASS'
        print(hstool_precision().to_string(index=False), flush=True)
        print(capscan_own_stool().to_string(index=False), flush=True)
        print(predict_crossomics()[['variant', 'level', 'layer_a', 'layer_b', 'rho', 'p']].to_string(index=False), flush=True)
        print(coverage_summary().to_string(index=False), flush=True)
        dump(OUT / 'ENDPOINTS_COMPLETE.json', dict(status='PASS', code_sha256=sha(__file__), job=os.environ.get('SLURM_JOB_ID'),
             outputs={p.name: sha(p) for p in OUT.iterdir() if p.is_file() and p.suffix == '.tsv'}))
    else:
        raise SystemExit(mode)
    print('DONE', mode, flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
