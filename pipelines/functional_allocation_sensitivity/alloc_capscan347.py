"""A3: CapScan part of the GH-029 functional-allocation sensitivity (iteration1/proposals_methods/alloc_sensitivity,
job 6048208) on the 331 (control) or 347 (participant 15 added) roster.

Imports the unmodified code/alloc_sensitivity.py (sha256 f62592df..., = ENDPOINTS_COMPLETE.json code_sha256) and reuses
unit_tests / projector / build_arm / frame unchanged; only module paths are redirected (CONS = annotated-origin
consensus folder, FROZEN = profiles folder for the 'equal' reproduction assertion, OUT = this run's output).
Edits (in this wrapper only):
 E1 capscan_own_stool: hard-coded participants range(1,15), 14 and 13 generalized to the participants observed in all
    five groups among positive acquisitions (C2 capscan_fig5bcd.py rule); permutation draws keep seed 20260910 with
    shape (99999, n), so n=14 reproduces production draws exactly; the rotated-stool negative control uses roll by 1.
 E2 design / positivity from config: control = A33 Figure5/design.tsv + A33 microbial_OG.tsv (production);
    347 = C2 roster/design_347.tsv + C2 refit OG_other_unresolved profile (C2 fig5bcd rule; identical positive set on 331).
 E3 positive-control reference from config: control = production capscan_own_stool.tsv and CAPSCAN rows of
    allocated_fraction_medians.tsv (must be equal); 347 = C2 fig5/fig5bcd_347/similarity_summary.tsv for the 'equal' variant.
 E4 HSTOOL and PREDICT endpoints are not run (no CapScan input; their printed values are unchanged).
usage: alloc_capscan347.py <config.json>"""
import importlib.util, json, sys, os, hashlib
from pathlib import Path
import numpy as np, pandas as pd
from scipy.spatial.distance import cdist

cfg = json.loads(Path(sys.argv[1]).read_text()); mode = cfg['mode']; assert mode in ('control331', 'fixed347')
CODE = Path(__file__).resolve().parent
def sha(p):
    h = hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda: f.read(8 << 20), b''): h.update(b)
    return h.hexdigest()
assert sha(CODE/'alloc_sensitivity.py') == 'f62592df0f5200a6c7fa0d048a68e3ab110e11e1153123d4056d519a2855c490'
for p, h in cfg['sha256'].items(): assert sha(p) == h, ('input changed', p)
spec = importlib.util.spec_from_file_location('alloc_sensitivity', CODE/'alloc_sensitivity.py'); a = importlib.util.module_from_spec(spec); spec.loader.exec_module(a)
out = Path(cfg['out']); out.mkdir(parents=True, exist_ok=False)
a.CONS = Path(cfg['consensus_root']); a.FROZEN = Path(cfg['frozen_profiles_root']); a.OUT = out
ARM = 'CAPSCAN_noclust'; REG = ['duodenum', 'jejunum', 'ileum', 'colon', 'stool']
a.unit_tests(); a.build_arm(ARM)

design = pd.read_csv(cfg['design'], sep='\t').set_index('sample_id')
if cfg['positivity'] == 'refit_OG':
    pos = a.frame(ARM, 'OG_other_unresolved__equal').reindex(design.index)          # == C2 refit OG profile (asserted equal in build_arm)
else:
    pos = pd.read_csv(cfg['positivity'], sep='\t', index_col=0).loc[design.index]
roster = design.index[pos.sum(1) > 0]
obsr = design.loc[roster].groupby('participant').region.nunique(); parts = sorted(obsr.index[obsr.eq(5)]); n = len(parts)
if mode == 'control331': assert len(roster) == 330 and parts == list(range(1, 15))
ix = pd.MultiIndex.from_product([parts, REG], names=['participant', 'region'])
perms = np.argsort(np.random.default_rng(20260910).random((99999, n)), axis=1)
def cells(raw):
    d = design.loc[roster]; v = raw.loc[roster]; ok = v.sum(1) > 0
    rel = v.div(v.sum(1).where(ok), axis=0)
    spec = rel.groupby([d.participant, d.region, d.specimen_id]).mean()
    return spec.groupby(level=[0, 1]).mean().loc[ix], int((~ok).sum())
def stats(x, prev):
    keep = (x > 0).mean(0) >= prev; v = x[:, keep]; h = np.sqrt(v / v.sum(1, keepdims=True)); cube = h.reshape(n, 5, -1)
    dd = np.stack([cdist(cube[i, :4], cube[:, 4]).mean(0) for i in range(n)]); own = np.diag(dd); other = (dd.sum(1) - own) / (n - 1)
    obs = float((other - own).mean()); pdg = dd[np.arange(n)[None, :], perms].mean(1); null = (dd.sum(1).mean() - n * pdg) / (n - 1)
    p = (1 + int((abs(null) >= abs(obs) - 1e-14).sum())) / 100000
    return dict(prevalence=prev, features=int(keep.sum()), own_closer=int((own < other).sum()), mean_own=float(own.mean()),
                mean_other=float(other.mean()), mean_advantage=obs, p=p)
rows = []
for v in a.VARIANTS:
    cc, zero = cells(a.frame(ARM, f'OG_other_unresolved__{v}').reindex(design.index)); x = cc.to_numpy()
    for prev in [0., .1, .2]: rows.append(dict(variant=v, acquisitions_without_allocated_signal=zero, participants=n, **stats(x, prev)))
cc, _ = cells(a.frame(ARM, 'OG_other_unresolved__equal').reindex(design.index))
cube = cc.to_numpy().reshape(n, 5, -1).copy(); cube[:, 4] = cube[np.roll(np.arange(n), 1), 4]
rows.append(dict(variant='NEGATIVE_CONTROL_stool_rotated', acquisitions_without_allocated_signal=0, participants=n, **stats(cube.reshape(5 * n, -1), 0.)))
own = pd.DataFrame(rows); own.to_csv(out/'capscan_own_stool.tsv', sep='\t', index=False)
c = pd.read_csv(out/'profiles'/ARM/'coverage.tsv', sep='\t'); c = c[c.profile == 'other_unresolved']
c['allocated_fraction'] = c.allocated / c.total_LFQ.where(c.total_LFQ > 0)
med = c.groupby(['arm', 'level', 'variant']).allocated_fraction.median().reset_index(); med.to_csv(out/'allocated_fraction_medians_CAPSCAN.tsv', sep='\t', index=False)
checks = {}
cols = ['features', 'own_closer', 'mean_own', 'mean_other', 'mean_advantage', 'p']
if mode == 'control331':
    ref = pd.read_csv(cfg['ref_own_stool'], sep='\t'); mine = own.drop(columns='participants')
    assert list(ref.variant) == list(mine.variant) and np.allclose(ref[cols + ['prevalence', 'acquisitions_without_allocated_signal']], mine[cols + ['prevalence', 'acquisitions_without_allocated_signal']], rtol=1e-12, atol=0), 'own-stool not reproduced'
    rm = pd.read_csv(cfg['ref_alloc_medians'], sep='\t'); rm = rm[rm.arm == ARM].reset_index(drop=True)
    assert rm[['arm', 'level', 'variant']].equals(med[['arm', 'level', 'variant']]) and np.allclose(rm.allocated_fraction, med.allocated_fraction, rtol=1e-12, atol=0), 'medians not reproduced'
    checks['production_capscan_own_stool_reproduced'] = checks['production_CAPSCAN_allocated_medians_reproduced'] = True
else:
    ref = pd.read_csv(cfg['ref_similarity'], sep='\t').set_index('prevalence'); e = own[own.variant == 'equal'].set_index('prevalence')
    assert np.allclose(e.loc[ref.index, cols], ref[cols], rtol=1e-9, atol=1e-12), 'equal variant does not reproduce C2 fig5bcd_347'
    assert int(ref.participants.iloc[0]) == n == 15
    checks['C2_fig5bcd_347_equal_reproduced'] = True
neg = own[own.variant.str.startswith('NEGATIVE')].iloc[0]; assert neg.p > 0.05, 'negative control unexpectedly significant'
(out/'COMPLETE.json').write_text(json.dumps(dict(status='PASS', mode=mode, job=os.environ.get('SLURM_JOB_ID'), participants=parts,
    roster_acquisitions=int(len(roster)), design_rows=int(len(design)), checks=checks, negative_control=dict(own_closer=int(neg.own_closer), p=float(neg.p)),
    config=cfg, wrapper_sha256=sha(__file__), outputs={p.name: sha(p) for p in out.iterdir() if p.is_file()}), indent=2, default=str) + '\n')
print(own.to_string(index=False)); print(med.to_string(index=False)); print(mode, 'PASS', flush=True)
