"""C1 cross-fit under PEP<=0.01, two options, read-only on iteration4/crossfit/outputs.
 baseline: pooled LDA rescoring + picked-peptide q over the 10 fold PSM tables (iteration4 design/sage_lda.pooled_q);
           must reproduce <ARM>_POOLED_ACCEPTED.txt (identity recorded per arm).
 opt_a:    baseline pooled q<=0.01 AND the per-fold Sage posterior_error (log10 PEP)<=-2 on the SAME PSM.
           Sentinel posterior_error==1.0 (unset) excluded.
 opt_b:    baseline pooled q<=0.01 AND log10(PEP_pooled)<=-2, PEP_pooled re-estimated from the pooled LDA scores
           with the Sage PSM PEP model (variant fixed by the 'validate' mode against Sage's own outputs).
           If the pooled rescoring falls back (no LDA), Sage would leave PEP unset -> zero; recorded.
Modes: validate <i>   reconstruct Sage posterior_error from Sage's own discriminant scores on the fold searches of roster row i
       run <chunk>    rows chunk*CH .. chunk*CH+CH-1
"""
import sys, json, hashlib, os
from pathlib import Path
import numpy as np, pandas as pd
IT4 = Path('LAB_ROOT/review/v15_integration/jess_integration_2026-09-29/analysis_followup/manuscript_completion_2026-10-01/biology_refit/claude_review_2026-10-03/iteration4/crossfit')
sys.path.insert(0, str(Path(__file__).resolve().parent))  # sage_q, sage_lda: identical copies of iteration4/crossfit/design
import sage_q, sage_lda
H = Path(os.environ.get('FASTALAKE_OUT_DIR', Path(__file__).resolve().parent))  # run directory: PEP_MODEL.json, per_acq/
ROSTER = IT4.parents[1] / 'iteration3/disjoint/roster.tsv'
ARMS = ['FL', 'FULL', 'MG']; CH = 29
VARIANT = json.loads((H / 'PEP_MODEL.json').read_text()) if (H / 'PEP_MODEL.json').exists() else None

def sha(p):
    """sha256 hex digest of a file."""
    h = hashlib.sha256(); h.update(Path(p).read_bytes()); return h.hexdigest()

def strip(pep):
    """Canonical peptide: modifications removed, I mapped to L."""
    import re
    p = re.sub(r'\[[^\]]*\]', '', pep).replace('I', 'L')
    assert p.isalpha() and p.isupper(), pep
    return p

def roster():
    """The 578-acquisition matched roster (iteration3/disjoint/roster.tsv) as a list of dicts."""
    return list(pd.read_csv(ROSTER, sep='\t', dtype=str).to_dict('records'))

def fold_table(path):
    """Read one fold's Sage PSM table (results.sage.tsv.gz)."""
    d = sage_lda.read(path)
    d['posterior_error'] = pd.read_csv(path, sep='\t', usecols=['posterior_error'])['posterior_error'].astype(np.float64).values
    return d

def pep_model(scores, decoy, v):
    """Sage PSM PEP model (KDE over discriminant scores, settings from PEP_MODEL.json): log10 PEP per PSM."""
    sel = np.ones(len(scores), bool)
    b, lo, step = sage_q.estimator(scores[sel].astype(np.float64), decoy[sel], bins=v['bins'], bw_adjust=v['bw'], monotonic=v['monotonic'])
    pe = sage_q.posterior_error(b, lo, step, scores.astype(np.float64))
    with np.errstate(divide='ignore'):
        return np.log10(pe).astype(np.float32)

def accepted(d, mask):
    """Canonical peptides of 5-50 and 9-50 residues among rank-one targets passing the mask."""
    m = (d['rank'] == 1) & (d['label'] == 1) & mask
    acc = {strip(p) for p in d.loc[m, 'peptide']}
    acc = {p for p in acc if 5 <= len(p) <= 50}
    return dict(peptides_5_50=len(acc), peptides_9_50=sum(len(p) >= 9 for p in acc)), acc

def validate(i):
    """Rebuild Sage posterior_error from Sage's own scores for the fold searches of roster row i and compare."""
    r = roster()[i]; out = IT4 / f"outputs/{r['cohort']}/{r['sample']}"
    res = []
    for arm in ARMS:
        for k in range(10):
            d = fold_table(out / f'fold{k}/{arm}/results.sage.tsv.gz')
            ok = d.posterior_error.values != 1.0
            s = d.sage_discriminant_score.values.astype(np.float32); dec = (d.label == -1).values
            for mono in (True, False):
                for bw in (1.0, 2.0):
                    v = dict(bins=1000, bw=bw, monotonic=mono)
                    p = pep_model(s, dec, v)
                    diff = np.abs(p[ok] - d.posterior_error.values[ok])
                    res.append(dict(arm=arm, fold=k, **v, n=int(ok.sum()), max_abs=float(np.nanmax(diff)), p99=float(np.nanpercentile(diff, 99)),
                                    gate_disagree=int(((p[ok] <= -2) != (d.posterior_error.values[ok] <= -2)).sum())))
    df = pd.DataFrame(res); df.to_csv(H / f'validate_{i}.tsv', sep='\t', index=False)
    print(df.groupby(['monotonic', 'bw'])[['max_abs', 'p99', 'gate_disagree']].max().to_string(), flush=True)

def run_one(i, r):
    """Baseline, option a and option b accepted-peptide counts for the three arms of roster row i."""
    out = IT4 / f"outputs/{r['cohort']}/{r['sample']}"
    cj = json.loads((out / 'COMPLETE.json').read_text())
    width = float(json.loads(Path(cj['config']).read_text())['precursor_tol']['ppm'][1] - json.loads(Path(cj['config']).read_text())['precursor_tol']['ppm'][0])
    rec = dict(i=i, cohort=r['cohort'], sample=r['sample'], subject_id=r['subject_id'], group=r['group'])
    for arm in ARMS:
        d = pd.concat([fold_table(out / f'fold{k}/{arm}/results.sage.tsv.gz') for k in range(10)], ignore_index=True)
        assert d.scannr.is_unique
        q, s, how = sage_lda.pooled_q(d, width)
        base_c, base = accepted(d, q <= 0.01)
        frozen = set(x for x in (out / f'{arm}_POOLED_ACCEPTED.txt').read_text().split('\n') if x)
        a_c, a = accepted(d, (q <= 0.01) & (d.posterior_error.values <= -2))
        dec = (d.label == -1).values
        if how == 'lda':
            lp = pep_model(s, dec, VARIANT)
            b_c, b = accepted(d, (q <= 0.01) & (lp <= -2))
        else:
            b_c = dict(peptides_5_50=0, peptides_9_50=0); b = set()
        assert a <= base and b <= base
        rec[arm] = dict(rescoring=how, baseline=base_c, baseline_equals_frozen=(base == frozen), baseline_frozen_symdiff=len(base ^ frozen),
                        opt_a=a_c, opt_b=b_c, unset_pep_psms=int((d.posterior_error.values == 1.0).sum()),
                        unset_pep_psms_q01_rank1_target=int(((d.posterior_error.values == 1.0) & (q <= .01) & (d['rank'] == 1) & (d.label == 1)).sum()))
    return rec

def run(chunk):
    """Run rows chunk*CH to chunk*CH+CH-1 and write per_acq/<i>.json."""
    rows = roster(); assert VARIANT is not None, 'PEP_MODEL.json missing'
    for i in range(chunk * CH, min(len(rows), chunk * CH + CH)):
        f = H / 'per_acq' / f'{i:03d}.json'
        if f.exists(): continue
        rec = run_one(i, rows[i]); rec['job'] = os.getenv('SLURM_JOB_ID'); rec['code_sha256'] = sha(__file__); rec['pep_model'] = VARIANT
        f.write_text(json.dumps(rec, indent=1)); print('done', i, {a: rec[a]['baseline_equals_frozen'] for a in ARMS}, flush=True)

if __name__ == '__main__':
    {'validate': validate, 'run': run}[sys.argv[1]](int(sys.argv[2]))
