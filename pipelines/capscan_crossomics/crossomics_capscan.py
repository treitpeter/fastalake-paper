"""C2 copy of biology_refit/crossomics.py (refit mode), CapScan section only.

Edits relative to production (everything else is copied line for line):
 E1 inputs are a JSON config (profiles, groups, design, metabolites, positivity roster, native DNA distances, output);
 E2 PREDICT section omitted (BH families are per cohort/scope/prevalence, so CapScan p/q are unaffected);
 E3 hard-coded DONORS=[1,2,3,4,6,8,9,10,11,12,13,14] and assert len(donors14)==14 replaced by the rule they
    encoded: donors14 = participants with protein and metabolite cells in all five groups; DONORS = those of
    donors14 that also have DNA cells in all five groups;
 E4 contributing-specimen roster: production uses positivity of the old A33 OG profile; config 'positivity'
    selects that file (control) or the refit OG profile (347; the two positive sets are identical on the 331);
 E5 native DNA distances read from config (production: Crossomics/native_CAPSCAN_distance_matrices.npz);
 E6 refit-vs-cross_baseline equality checks dropped (they compare to the 331 population by design).
usage: crossomics_capscan.py <config.json>"""
import ast, json, math, sys, hashlib
from itertools import combinations
from pathlib import Path
import numpy as np, pandas as pd
from scipy import sparse
from scipy.spatial.distance import pdist, squareform
from scipy.stats import rankdata, spearmanr
from statsmodels.stats.multitest import multipletests

cfg = json.loads(Path(sys.argv[1]).read_text()); out = Path(cfg['out']); out.mkdir(parents=True, exist_ok=False)
R = Path('LAB_ROOT')
X = R/'review/v15_integration/jess_integration_2026-09-29/source_data/current/Crossomics'
REG = ['duodenum', 'jejunum', 'ileum', 'colon', 'stool']
def rd(p, **kw): return pd.read_csv(p, sep='\t', **kw)
def wr(d, p): d.to_csv(p, sep='\t', index=False)
def frame(folder, key):
    a = json.loads((Path(folder)/(key + '_axes.json')).read_text()); v = sparse.load_npz(Path(folder)/(key + '.npz'))
    return pd.DataFrame(v.toarray(), index=a['samples'], columns=a['features'])
def functions(path, names, env):
    tree = ast.parse(Path(path).read_text()); defs = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name in names]
    assert {n.name for n in defs} == set(names)
    exec(compile(ast.Module(body=defs, type_ignores=[]), str(path), 'exec'), env); return env

SEED = 20260907; NPERM = 1999; coverage = []; results = []; matrices = {}; rosters = []
functions(X/'provenance/check_prevalence.py', ['relative', 'logmet', 'transform', 'evaluate'], globals())
_original_evaluate = evaluate
def evaluate(co, scope, lev, frames, ci, prev, pairset=None):
    _original_evaluate(co, scope, lev, frames, ci, prev, pairset)
    prefix = '__'.join([co, scope, lev, str(prev)])
    for layer in frames:
        matrices[prefix + '__n' + str(len(ci)) + '__' + layer] = matrices.pop(prefix + '__' + layer)
ns = dict(np=np, pd=pd, sparse=sparse, math=math, SEED=SEED, combinations=combinations, rankdata=rankdata, squareform=squareform, spearmanr=spearmanr)
functions(X/'native/native_pyramids.py', ['profile_distance', 'bh', 'correlations'], ns)
def fresh(co, lev):
    return frame(cfg['groups'], 'supported_human_groups') if lev == 'host' else frame(cfg['profiles'], lev + '_other_unresolved')

des = rd(cfg['design'], index_col=0); des.specimen = des.specimen.astype(str)
met = rd(cfg['metabolites'], index_col=0); met.index = met.index.astype(str)
if cfg['positivity'] == 'refit_OG': posog = fresh('CAPSCAN', 'OG').reindex(des.index)
else: posog = rd(cfg['positivity'], index_col=0).reindex(des.index)
newog = fresh('CAPSCAN', 'OG').reindex(des.index)
wr(pd.DataFrame(dict(acquisition=des.index, previous_positive=posog.sum(1).gt(0).to_numpy(), current_positive=newog.sum(1).gt(0).to_numpy())), out/'eligibility_audit.tsv')
valid = set(des.loc[posog.sum(1) > 0, 'specimen']); met = met.loc[met.index.isin(valid)]
dspec = des.groupby('specimen').first().loc[met.index]; host = relative(fresh('CAPSCAN', 'host').loc[des.index]).groupby(des.specimen).median(); cells = {}
for lev in ['OG', 'KO', 'Pathway']:
    micro = relative(fresh('CAPSCAN', lev).loc[des.index]).groupby(des.specimen).median(); cells[lev] = {}
    for k, f in dict(microbial_metaP=micro, host_metaP=host, metabolites=logmet(met)).items():
        si = f.index.intersection(dspec.index); cells[lev][k] = f.loc[si].groupby([dspec.loc[si, 'subject_id'], dspec.loc[si, 'group']]).mean()
common = set.intersection(*(set(x.index) for fs in cells.values() for x in fs.values()))
donors14 = sorted(d for d in {c[0] for c in common} if all((d, r) in common for r in REG))          # E3
runs = rd(R/'analyses/capscan_metag_analysis_2026-09-10/manifests/runs.tsv'); runs.specimen = runs.specimen.astype(str); dna = {}
for lev in ['OG', 'KO', 'Pathway']:
    d = rd(R/'analyses/capscan_metag_analysis_2026-09-10/aggregate'/('metaG_' + lev + '.tsv'), index_col=0); d = d.loc[d.sum(1) > 0]; d.index = runs.set_index('run_accession').loc[d.index, 'specimen'].values
    d = d.groupby(level=0).mean(); si = d.index.intersection(dspec.index); dna[lev] = relative(d.loc[si]).groupby([dspec.loc[si, 'subject_id'], dspec.loc[si, 'group']]).mean()
have = set.intersection(*(set(v.index) for v in dna.values()))
DONORS = [d for d in donors14 if all((d, r) in have for r in REG)]                                   # E3
print('protein/metabolite complete participants', len(donors14), donors14, '| DNA complete', len(DONORS), DONORS, flush=True)
ci = lambda ds: pd.MultiIndex.from_product([ds, REG], names=['subject_id', 'group'])
assert all(set(ci(DONORS)) <= set(v.index) for v in dna.values())
for lev in ['OG', 'KO', 'Pathway']:
    for prev in [0., .1]:
        evaluate('CAPSCAN', 'displayed_mixed12_14', lev, cells[lev], ci(donors14), prev)
        evaluate('CAPSCAN', 'displayed_mixed12_14', lev, dict(cells[lev], metaG=dna[lev]), ci(DONORS), prev, [('metaG', k) for k in ['microbial_metaP', 'host_metaP', 'metabolites']])
        evaluate('CAPSCAN', 'matched12', lev, dict(metaG=dna[lev], **cells[lev]), ci(DONORS), prev)
    print('functional CAPSCAN', lev, flush=True)
res = pd.DataFrame(results); res['q'] = np.nan
for key, g in res.groupby(['cohort', 'scope', 'prevalence']):
    assert len(g) == 16; res.loc[g.index, 'q'] = multipletests(g.p, method='fdr_bh')[1]
wr(res, out/'functional_correlations.tsv'); wr(pd.DataFrame(coverage).drop_duplicates(), out/'feature_coverage.tsv'); wr(pd.DataFrame(rosters).drop_duplicates(), out/'signal_retention.tsv'); np.savez_compressed(out/'functional_distances.npz', **matrices)
wr(dspec.reset_index(), out/'capscan_specimen_roster.tsv')
native = np.load(cfg['native_dna']); nativerows = []; nativemats = {}; qualities = {}; controlerrs = {}
hg = frame(cfg['groups'], 'supported_human_groups'); ng = frame(cfg['groups'], 'other_unresolved_groups'); fs = {}
for k, f in [('Non-host protein', ng), ('Host protein', hg)]:
    byspec = relative(f.loc[des.index]).groupby(des.specimen).median(); si = byspec.index.intersection(dspec.index)
    fs[k] = byspec.loc[si].groupby([dspec.loc[si, 'subject_id'], dspec.loc[si, 'group']]).mean().loc[ci(DONORS)]
fs['Metabolites'] = logmet(met).groupby([dspec.subject_id, dspec.group]).mean().loc[ci(DONORS)]
for label, prev in [('primary_prevalence10', .1), ('sensitivity_all_features', 0.)]:
    ds = {'DNA': native[label + '__DNA']}; assert ds['DNA'].shape == (5*len(DONORS),)*2, ('native DNA population differs', ds['DNA'].shape, len(DONORS)); qualities[label] = {}
    for k, f in fs.items():
        ds[k], q = ns['profile_distance'](f.to_numpy() if k == 'Metabolites' else sparse.csr_matrix(f.to_numpy()), k, prev, capscan=True); qualities[label][k] = q
        if label + '__' + k in native.files: controlerrs[label + '__' + k] = float(np.max(abs(ds[k] - native[label + '__' + k])))
    tab = ns['correlations'](ds, 'CAPSCAN'); tab.insert(1, 'analysis', label); nativerows.append(tab)
    nativemats.update({'CAPSCAN__' + label + '__' + k: v for k, v in ds.items()})
nat = pd.concat(nativerows, ignore_index=True); wr(nat, out/'native_correlations.tsv'); np.savez_compressed(out/'native_distances.npz', **nativemats)
(out/'native_quality.json').write_text(json.dumps(qualities, indent=2, default=float))
(out/'COMPLETE.json').write_text(json.dumps(dict(status='PASS', config=cfg, protein_metabolite_participants=donors14, dna_participants=DONORS,
    functional_rows=len(res), native_rows=len(nat), native_distance_max_abs_diff_vs_supplied=controlerrs,
    script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()), indent=2, default=str) + '\n')
print(res[res.scope.eq('matched12') & res.prevalence.eq(0)].to_string(index=False)); print(nat.to_string(index=False))
