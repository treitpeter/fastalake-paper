"""C3 (JAG30): repeated-measures test of CapScan FastaLake-over-metagenome gains by sampling group.
Input: iteration4/crossfit/run/PER_ACQUISITION.tsv (CAPSCAN rows, 282; labels = frozen roster, provisional; participant 15 stool-only).
Outcomes per acquisition: rel = 100*(FL/MG-1); logr = ln(FL/MG).
 cross-fitted (primary): FL_pooled_L / MG_pooled_L ; same-spectrum (v25): frozen_fl_L / frozen_mg_L  (L=5: 5-50, L=9: 9-50)
Model: y ~ C(group, ref=stool) + (1|participant), ML fit (statsmodels MixedLM). Omnibus: LRT vs y ~ 1 + (1|participant), 4 df,
plus Wald 4 df. Contrasts: pooled fluid (mean of 4 region coefs) vs stool; 6 pairwise region contrasts, Holm; region omnibus (3 df Wald).
Robustness (participant level): per participant x group median; paired Wilcoxon fluid-vs-stool (participant median over all fluid
acquisitions); exact sign-flip permutation of mean paired difference; pairwise region Wilcoxon + sign-flip (participants with both), Holm;
Friedman across 4 regions (participants with all 4).
"""
import sys, itertools, numpy as np, pandas as pd, scipy.stats as st, statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests
IN, OUT = sys.argv[1], sys.argv[2]
a = pd.read_csv(IN, sep='\t', dtype={'subject_id': str}); a = a[a.cohort == 'CAPSCAN'].copy()
assert len(a) == 282, len(a)
RMAP = {'duodenum': 'R1', 'jejunum': 'R2', 'ileum': 'R3', 'colon': 'R4', 'stool': 'stool'}
a['g'] = a.group.map(RMAP); assert a.g.notna().all()
print('n per group', a.g.value_counts().to_dict(), 'participants', a.subject_id.nunique(), flush=True)
SPEC = {'crossfit': ('FL_pooled_{L}', 'MG_pooled_{L}'), 'same_spectrum': ('frozen_fl_{L}', 'frozen_mg_{L}')}
REG = ['R1', 'R2', 'R3', 'R4']
rows, prow = [], []

import warnings; warnings.filterwarnings('ignore')
def fit(f, d):
    # statsmodels MixedLM can report convergence at a bad point (lbfgs left the intercept at 0 on the rel scale);
    # fit with four optimisers and keep the highest ML log-likelihood among fits with a positive-definite fixed-effect covariance.
    best, last = None, None
    for meth in ['lbfgs', 'bfgs', 'powell', 'nm']:
        try:
            m = smf.mixedlm(f, d, groups=d.subject_id).fit(reml=False, method=meth, maxiter=5000)
            cp = m.cov_params().loc[m.fe_params.index, m.fe_params.index].values
            if m.converged and np.isfinite(m.llf) and np.all(np.isfinite(cp)) and np.all(np.linalg.eigvalsh(cp) > 0):
                if best is None or m.llf > best.llf + 1e-9:
                    m.method_used = meth; best = m
        except Exception as e:
            last = e
    if best is None: raise RuntimeError(f'no method converged for {f}: {last}')
    return best

def signflip(d):
    d = np.asarray(d, float); n = len(d)
    if n == 0: return np.nan
    obs = abs(d.mean()); signs = np.array(list(itertools.product([1, -1], repeat=n)))
    return float(np.mean(np.abs((signs * d).mean(1)) >= obs - 1e-12))

for view, (fl, mg) in SPEC.items():
    for L in ['5', '9']:
        F = a[fl.format(L=L)].astype(float); M = a[mg.format(L=L)].astype(float)
        for scale in ['rel', 'logr']:
            d = a[['subject_id', 'g']].copy()
            d['y'] = 100 * (F / M - 1) if scale == 'rel' else np.log(F / M)
            d = d[np.isfinite(d.y)]
            base = dict(view=view, length=f'{L}-50', scale=scale, n_acq=len(d), n_part=d.subject_id.nunique())
            sd = d.y.std(); d['ys'] = d.y / sd   # fit on SD-scaled y for numerical stability; estimates rescaled by sd below
            m1 = fit("ys ~ C(g, Treatment('stool'))", d); m0 = fit("ys ~ 1", d)
            lr = 2 * (m1.llf - m0.llf); base['fit_method'] = m1.method_used; base['participant_var_frac'] = float(m1.cov_re.iloc[0, 0] / (m1.cov_re.iloc[0, 0] + m1.scale)); p_lrt = st.chi2.sf(lr, 4)
            names = [f"C(g, Treatment('stool'))[T.{r}]" for r in REG]
            b = m1.fe_params[names].values * sd; V = m1.cov_params().loc[names, names].values * sd ** 2
            w4 = float(b @ np.linalg.solve(V, b)); p_w4 = st.chi2.sf(w4, 4)
            # region-only omnibus (R2-R1, R3-R1, R4-R1)
            C3 = np.array([[-1, 1, 0, 0], [-1, 0, 1, 0], [-1, 0, 0, 1]], float)
            cb = C3 @ b; w3 = float(cb @ np.linalg.solve(C3 @ V @ C3.T, cb)); p_w3 = st.chi2.sf(w3, 3)
            def con(c, label, fam):
                c = np.asarray(c, float); est = float(c @ b); se = float(np.sqrt(c @ V @ c)); z = est / se
                rows.append(dict(base, test=label, family=fam, estimate=est, se=se, ci_low=est - 1.96 * se, ci_high=est + 1.96 * se,
                                 stat=z, df=1, p=2 * st.norm.sf(abs(z))))
            rows.append(dict(base, test='omnibus_5groups_LRT', family='omnibus', estimate=np.nan, se=np.nan, ci_low=np.nan, ci_high=np.nan, stat=lr, df=4, p=p_lrt))
            rows.append(dict(base, test='omnibus_5groups_Wald', family='omnibus', estimate=np.nan, se=np.nan, ci_low=np.nan, ci_high=np.nan, stat=w4, df=4, p=p_w4))
            rows.append(dict(base, test='omnibus_regions_only_Wald', family='omnibus', estimate=np.nan, se=np.nan, ci_low=np.nan, ci_high=np.nan, stat=w3, df=3, p=p_w3))
            rows.append(dict(base, test='stool_intercept', family='level', estimate=float(m1.fe_params['Intercept']) * sd, se=float(m1.bse['Intercept']) * sd,
                             ci_low=np.nan, ci_high=np.nan, stat=np.nan, df=np.nan, p=np.nan, participant_var=float(m1.cov_re.iloc[0, 0]) * sd ** 2, resid_var=float(m1.scale) * sd ** 2))
            con([.25] * 4, 'fluid_pooled_minus_stool', 'fluid_vs_stool')
            for i, r in enumerate(REG):
                c = np.zeros(4); c[i] = 1; con(c, f'{r}_minus_stool', 'region_vs_stool')
            for i, j in itertools.combinations(range(4), 2):
                c = np.zeros(4); c[j] = 1; c[i] = -1; con(c, f'{REG[j]}_minus_{REG[i]}', 'region_pairwise')
            # participant-level robustness
            pm = d.groupby(['subject_id', 'g']).y.median().unstack()
            fm = d[d.g != 'stool'].groupby('subject_id').y.median()
            both = pm.stool.dropna().index.intersection(fm.index); dd = (fm[both] - pm.stool[both]).values
            prow.append(dict(base, test='fluid_minus_stool', family='fluid_vs_stool', n_pairs=len(dd), mean_diff=dd.mean(), median_diff=np.median(dd),
                             n_pos=int((dd > 0).sum()), p_wilcoxon=st.wilcoxon(dd).pvalue, p_signflip=signflip(dd)))
            for r in REG:
                x = pm[[r, 'stool']].dropna(); dd = (x[r] - x.stool).values
                prow.append(dict(base, test=f'{r}_minus_stool', family='region_vs_stool', n_pairs=len(dd), mean_diff=dd.mean(), median_diff=np.median(dd),
                                 n_pos=int((dd > 0).sum()), p_wilcoxon=st.wilcoxon(dd).pvalue, p_signflip=signflip(dd)))
            for i, j in itertools.combinations(REG, 2):
                x = pm[[i, j]].dropna(); dd = (x[j] - x[i]).values
                prow.append(dict(base, test=f'{j}_minus_{i}', family='region_pairwise', n_pairs=len(dd), mean_diff=dd.mean(), median_diff=np.median(dd),
                                 n_pos=int((dd > 0).sum()), p_wilcoxon=st.wilcoxon(dd).pvalue, p_signflip=signflip(dd)))
            x = pm[REG].dropna(); fr = st.friedmanchisquare(*[x[r] for r in REG])
            prow.append(dict(base, test='friedman_4regions', family='omnibus', n_pairs=len(x), mean_diff=np.nan, median_diff=np.nan, n_pos=np.nan,
                             p_wilcoxon=fr.pvalue, p_signflip=np.nan))
R = pd.DataFrame(rows); P = pd.DataFrame(prow)
key = ['view', 'length', 'scale', 'family']
for df_, cols in [(R, ['p']), (P, ['p_wilcoxon', 'p_signflip'])]:
    for c in cols:
        df_[c + '_holm'] = np.nan
        for k, g in df_.groupby(key):
            if k[3] in ('region_vs_stool', 'region_pairwise'):
                ok = g[c].notna()
                df_.loc[g.index[ok], c + '_holm'] = multipletests(g.loc[ok, c], method='holm')[1]
R.to_csv(f'{OUT}/MIXED_MODEL.tsv', sep='\t', index=False); P.to_csv(f'{OUT}/PARTICIPANT_LEVEL.tsv', sep='\t', index=False)
pd.set_option('display.width', 250); pd.set_option('display.max_rows', 500)
print(R.drop(columns=['participant_var', 'resid_var'], errors='ignore').round(4).to_string())
print(P.round(4).to_string())
