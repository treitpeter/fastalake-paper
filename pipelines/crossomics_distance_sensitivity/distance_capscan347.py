"""A4 copy of biology_refit/distance_sensitivity.py (refit mode), CapScan section only.

Edits relative to production (loop body, variant/metric/adjustment/outer grid, gram/distance/correlate functions
and the 1,999-permutation design are unchanged; functions are exec'd from byte-identical copies in code/):
 E1 inputs are a JSON config (groups, design, specimen roster, donors, DNA native cells, metabolites, output).
 E2 PREDICT omitted. correlate() applies BH within each cohort x features x distance x adjustment x outer table
    (six CapScan tests), so CapScan p/q do not depend on PREDICT, and PREDICT rows do not involve participant 15.
 E3 hard-coded common.DONORS (12 participants) replaced by the config list (347: the 13 DNA-complete participants
    recorded in C2 fig5g/crossomics_347/COMPLETE.json dna_participants).
 E4 DNA distances are computed from the DNA native cells with the production gram_for/gram_distance. Production
    read them from the cached method_comparison distance_matrices.npz, which exists only for the 12-participant
    population. Control mode asserts recomputed == cached (|diff| < 1e-7, production tolerance).
 E5 Metabolite profiles: log10 metabolite table restricted to the specimen roster, mean per participant-region
    (rule of native_pyramids.other_layers / C2 crossomics). Control mode asserts this equals the production
    prepared Metabolites_log_profiles.npy and then uses the npy, as production did.
 E6 the 'previous' table merged one-to-one is the production distance_refit CAPSCAN rows (production merged the
    older method_comparison table); no equality is asserted in 347 mode.
 E7 347 mode cross-checks region-adjusted Hellinger DNA and metabolite distances against the independently
    built cohort_reconciliation current331_plus16_all native distances.
usage: distance_capscan347.py <config.json>"""
import json, sys, math, hashlib, ast, os
from pathlib import Path
from itertools import combinations
from types import SimpleNamespace
import numpy as np, pandas as pd
from scipy import sparse
from scipy.stats import rankdata, spearmanr, pearsonr
from scipy.spatial.distance import squareform, pdist

cfg = json.loads(Path(sys.argv[1]).read_text()); mode = cfg['mode']; assert mode in ['control331', 'fixed347']
CODE = Path(__file__).resolve().parent
def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''): h.update(b)
    return h.hexdigest()
def rd(p, **kw): return pd.read_csv(p, sep='\t', **kw)
def wr(d, p): d.to_csv(p, sep='\t', index=False)
def frame(folder, key):
    a = json.loads((Path(folder)/(key + '_axes.json')).read_text()); v = sparse.load_npz(Path(folder)/(key + '.npz'))
    return pd.DataFrame(v.toarray(), index=a['samples'], columns=a['features'])
def functions(path, names, env):
    tree = ast.parse(Path(path).read_text()); defs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    assert {n.name for n in defs} == set(names)
    exec(compile(ast.Module(body=defs, type_ignores=[]), str(path), 'exec'), env); return env

FROZEN = {'distance_sensitivity.py': 'f23a5468d394757c4ba73eee691bde3d19f0896a1716dfe6b37f1d7205992d65',
          'experiment.py': '382d8f70e3352f493e8b71b328d82bd9fa91086d390b7d945ecc26c78eb62d62',
          'native_pyramids.py': '64ba668b266b9be42bb0a3f273878a2b93daba53022cc3c5a70ab04d4d56c8da'}
for f, h in FROZEN.items(): assert sha(CODE/f) == h, ('frozen code changed', f)
for p, h in cfg['sha256'].items(): assert sha(p) == h, ('input changed', p)
out = Path(cfg['out']); out.mkdir(parents=True, exist_ok=False)

REG = ['duodenum', 'jejunum', 'ileum', 'colon', 'stool']
env = dict(np=np, pd=pd, sparse=sparse, math=math, rankdata=rankdata, spearmanr=spearmanr, pearsonr=pearsonr, pdist=pdist, squareform=squareform, combinations=combinations)
functions(CODE/'native_pyramids.py', ['bh'], env); env['native'] = SimpleNamespace(SEED=20260907, bh=env['bh'])
functions(CODE/'experiment.py', ['projection', 'relative', 'rank_gram', 'gram_distance', 'met_matrix', 'gram_for', 'rank_checks', 'correlate'], env); env['rank_checks']()

DONORS = [int(d) for d in cfg['donors']]
des = rd(cfg['design'], index_col=0); des.specimen = des.specimen.astype(str)
dspec = rd(cfg['specimen_roster'], index_col=0); dspec.index = dspec.index.astype(str)
ci = pd.MultiIndex.from_product([DONORS, REG], names=['subject_id', 'group'])
layers = {}
for key, folderkey in [('Non-host protein', 'other_unresolved_groups'), ('Host protein', 'supported_human_groups')]:
    raw = frame(cfg['groups'], folderkey)
    x = raw.loc[des.index]; rel = x.div(x.sum(1).replace(0, np.nan), axis=0); spec = rel.groupby(des.specimen).median(); si = spec.index.intersection(dspec.index)
    raw = spec.loc[si].groupby([dspec.loc[si, 'subject_id'], dspec.loc[si, 'group']]).mean().loc[ci]
    layers[key] = sparse.csr_matrix(raw.to_numpy())
met = rd(cfg['metabolites'], index_col=0); met.index = met.index.astype(str); assert set(dspec.index) <= set(met.index)
metcells = np.log10(met.loc[dspec.index].where(met.loc[dspec.index] > 0)).groupby([dspec.subject_id, dspec.group]).mean().loc[ci].to_numpy()
checks = {}
if mode == 'control331':
    npy = np.load(cfg['metabolites_npy']); assert npy.shape == metcells.shape and np.array_equal(np.isnan(npy), np.isnan(metcells))
    err = float(np.nanmax(abs(npy - metcells))); assert err < 1e-12, err; checks['metabolite_builder_vs_production_npy'] = err
    layers['Metabolites'] = npy
else: layers['Metabolites'] = metcells
dna = sparse.load_npz(cfg['dna_cells']).tocsr(); assert dna.shape[0] == len(ci), (dna.shape, len(ci))
cached = np.load(cfg['dna_cached']) if cfg.get('dna_cached') else None
recon = np.load(cfg['recon_native']) if cfg.get('recon_native') else None
prev = rd(cfg['previous_rows']); prev = prev[prev.cohort == 'CAPSCAN']
rows = []; archive = {}; quality = {}; errors = {}
co = 'CAPSCAN'
for variant in ['all_observed', 'prevalence10', 'wuyts_recipe']:
    for metric in ['hellinger', 'rank']:
        if variant == 'wuyts_recipe' and metric == 'hellinger': continue
        grams = {}; quality[variant + '__' + metric] = {}
        grams['DNA'], quality[variant + '__' + metric]['DNA'] = env['gram_for'](dna, 'DNA', variant, metric)                    # E4
        for key, x in layers.items(): grams[key], quality[variant + '__' + metric][key] = env['gram_for'](x, key, variant, metric)
        for adj in ['unadjusted', 'region_adjusted']:
            prefix = '__'.join([variant, metric, adj]); ds = {}
            for key in ['DNA', 'Non-host protein', 'Host protein', 'Metabolites']:
                ds[key] = env['gram_distance'](grams[key], adj == 'region_adjusted', metric == 'rank')
            if cached is not None:
                for key in ['DNA', 'Metabolites']:
                    e = float(abs(ds[key] - cached[prefix + '__' + key]).max()); assert e < 1e-7, (prefix, key, e); errors['cached__' + prefix + '__' + key] = e
            if recon is not None and metric == 'hellinger' and adj == 'region_adjusted' and variant in ('all_observed', 'prevalence10'):   # E7
                label = 'sensitivity_all_features' if variant == 'all_observed' else 'primary_prevalence10'
                for key in ['DNA', 'Metabolites']:
                    e = float(abs(ds[key] - recon[label + '__' + key]).max()); assert e < 1e-7, (label, key, e); errors['recon__' + label + '__' + key] = e
            archive.update({prefix + '__' + k: d for k, d in ds.items()})
            for outer in ['pearson', 'spearman']:
                t = env['correlate'](ds, co, outer); t.insert(1, 'features', variant); t.insert(2, 'distance', metric); t.insert(3, 'adjustment', adj); t.insert(4, 'outer_correlation', outer); rows.append(t)
table = pd.concat(rows, ignore_index=True); keys = ['cohort', 'features', 'distance', 'adjustment', 'outer_correlation', 'layer_a', 'layer_b']
comp = prev.merge(table, on=keys, validate='one_to_one', suffixes=('_previous', '_revised')); assert len(comp) == len(table) == len(prev) == 120
if mode == 'control331':
    for c in ['coefficient', 'p', 'p_positive_tail', 'q']: assert np.allclose(comp[c + '_previous'], comp[c + '_revised'], atol=1e-10, rtol=0), c
    checks['production_distance_refit_CAPSCAN_rows_reproduced'] = True
wr(comp, out/'CAPSCAN_method_comparison.tsv'); np.savez_compressed(out/'CAPSCAN_distances.npz', **archive)
(out/'CAPSCAN_quality.json').write_text(json.dumps(quality, indent=2, default=float) + '\n')
wr(table, out/'comparisons.tsv'); controlled = table[table.features != 'wuyts_recipe']; assert len(controlled) == 96; wr(controlled, out/'controlled_CAPSCAN_96.tsv')
(out/'COMPLETE.json').write_text(json.dumps(dict(status='PASS', mode=mode, job=os.environ.get('SLURM_JOB_ID'), donors=DONORS, n_profiles=len(ci),
    capscan_rows=len(table), capscan_controlled=len(controlled), checks=checks, distance_errors=errors, config=cfg,
    script_sha256=sha(__file__), outputs={p.name: sha(p) for p in out.iterdir() if p.is_file()}), indent=2, default=str) + '\n')
print(mode, 'PASS', len(table), 'rows,', len(DONORS), 'participants', flush=True)
