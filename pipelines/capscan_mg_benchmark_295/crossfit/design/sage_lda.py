"""Re-implementation of Sage 0.14.6 PSM rescoring (sage_core::ml::linear_discriminant::score_psms),
so that the PSMs of the 10 fold searches of one acquisition and arm can be rescored with ONE model,
as a single Sage run over all spectra would do, before the picked-peptide q (sage_q.peptide_q).

Why: each fold search fits its own LDA, and the raw discriminant scores of different folds are on
different scales. Pooling them unchanged gave invalid q values in the pilot (archived FULL FASTA pooled
over folds: 5,631 accepted peptides against 10,657 in the frozen single search; union 10,616).

Source read: github.com/lazear/sage v0.14.6 crates/sage/src/ml/{linear_discriminant,gauss,matrix,kde}.rs,
sage-cli/src/{main,output}.rs. Steps as in Sage:
 - precursor mass model: KDE posterior error of delta_mass (tsv precursor_ppm), bandwidth x2,
   bins = max(precursor tolerance width, 100) = 100, not monotone.
 - 20 features: rank, charge, ln1p(hyperscore), ln1p(delta_next), ln1p(delta_best), mass PEP,
   isotope_error, average_ppm (tsv fragment_ppm), ln1p(-poisson) (3.5 if not finite),
   ln1p(matched_intensity_pct), matched_peaks, ln1p(longest_b), ln1p(longest_y), longest_y/peptide_len,
   ln1p(peptide_len), missed_cleavages, aligned_rt, ims (tsv ion_mobility),
   sqrt(clamp(delta_rt_model, .001, .999)), sqrt(clamp(delta_ims_model (tsv delta_mobility), .001, .999)).
 - LDA: Sw = cov(decoy) + cov(target) (each /n_class), Sb = sum_c (m_c - xbar)(m_c - xbar)^T;
   M = Gauss-Jordan solve of (Sw + eps I) M = Sb ported operation for operation from gauss.rs, eps from
   1e-8 upward x10 until Sage's own left_solved() test passes; eigenvector by the power method
   from xbar/|xbar| (50 iterations, stop when |norm change| < 1e-8); sign so that target mean score >
   decoy mean score; score = features . eigenvector, stored as f32.
 - If the model fails (no solution or non-finite), Sage's fallback: ln1p(-poisson) as f32 + longest_y_pct/3.
f32 fields are parsed as float32 (as Sage stores them) and widened to f64; f64 fields as float64.
Validation is against Sage outputs (validate and k1pool modes). Sage's solve is numerically unstable on
these data (signed-pivot Gauss-Jordan, scatter condition number ~1e9, three constant features):
perturbing the inputs by one ulp can move scores by up to ~16 units, so bitwise agreement is not always
attainable; pool() records a two-perturbation diagnostic of this sensitivity per acquisition and arm.
"""
import numpy as np, pandas as pd
import sage_q

F64 = ['hyperscore', 'delta_next', 'delta_best', 'poisson']
F32 = ['precursor_ppm', 'isotope_error', 'fragment_ppm', 'matched_intensity_pct', 'aligned_rt', 'ion_mobility',
       'delta_rt_model', 'delta_mobility', 'longest_y_pct']
INT = ['rank', 'charge', 'matched_peaks', 'longest_b', 'longest_y', 'peptide_len', 'missed_cleavages', 'label']
COLS = ['scannr', 'peptide', 'sage_discriminant_score', 'peptide_q'] + F64 + F32 + INT


def read(path):
    d = pd.read_csv(path, sep='\t', usecols=COLS, dtype={'scannr': str, 'peptide': str, **{c: str for c in F32 + F64}})
    for c in F32:
        d[c] = d[c].astype(np.float32).astype(np.float64)
    for c in F64:
        d[c] = d[c].astype(np.float64)
    return d


def mass_pep(dm, decoy, width=40.0):
    b, lo, step = sage_q.estimator(dm, decoy, bins=int(np.ceil(abs(max(width, 100.0)))), bw_adjust=2.0, monotonic=False)
    return sage_q.posterior_error(b, lo, step, dm)


def features(d, width=40.0):
    decoy = (d['label'] == -1).to_numpy()
    with np.errstate(divide='ignore', invalid='ignore'):
        pois = np.log1p(-d['poisson'].to_numpy())
    pois = np.where(np.isfinite(pois), pois, 3.5)
    X = np.column_stack([
        d['rank'], d['charge'], np.log1p(d['hyperscore']), np.log1p(d['delta_next']), np.log1p(d['delta_best']),
        mass_pep(d['precursor_ppm'].to_numpy(), decoy, width), d['isotope_error'], d['fragment_ppm'], pois,
        np.log1p(d['matched_intensity_pct']), d['matched_peaks'], np.log1p(d['longest_b']), np.log1p(d['longest_y']),
        d['longest_y'].to_numpy(np.float64) / d['peptide_len'].to_numpy(np.float64), np.log1p(d['peptide_len']), d['missed_cleavages'],
        d['aligned_rt'], d['ion_mobility'], np.sqrt(np.clip(d['delta_rt_model'], .001, .999)), np.sqrt(np.clip(d['delta_mobility'], .001, .999))]).astype(np.float64)
    return X, decoy


def _gauss_inner(left, right, eps):
    """Sage ml/gauss.rs Gauss::solve_inner, operation for operation."""
    L = left.copy(); R = right.copy(); m, n = L.shape
    for i in range(n):
        L[i, i] += eps
    h = 0; k = 0
    while h < m and k < n:                       # echelon (partial pivoting on the largest signed value, last wins)
        mi, mv = 0, -np.finfo(np.float64).max
        for i in range(h, m):
            if L[i, k] >= mv:
                mi, mv = i, L[i, k]
        if L[mi, k] == 0.0:
            k += 1; continue
        if h != mi:
            L[[h, mi]] = L[[mi, h]]; R[[h, mi]] = R[[mi, h]]
        for i in range(h + 1, m):
            f = L[i, k] / L[h, k]
            L[i, k] = 0.0
            L[i, k + 1:] -= L[h, k + 1:] * f
            R[i, :] -= R[h, :] * f
        h += 1; k += 1
    for i in range(m - 1, -1, -1):               # reduce
        for j in range(n):
            x = L[i, j]
            if x == 0.0:
                continue
            L[i, j:] /= x; R[i, :] /= x
            break
    for i in range(m - 1, -1, -1):               # backfill
        for j in range(n):
            if L[i, j] == 0.0:
                continue
            for kk in range(i):
                f = L[kk, j] / L[i, j]
                L[kk, :] -= L[i, :] * f; R[kk, :] -= R[i, :] * f
            break
    for i in range(n):                           # left_solved
        for j in range(n):
            x = L[i, j]
            if i == j:
                if x != 1.0 and x != 0.0:
                    return None
            elif x > 1e-8:
                return None
    return R


def gauss_solve(left, right):
    eps = 1e-8
    while eps <= 1.0:
        r = _gauss_inner(left, right, eps)
        if r is not None:
            return r, eps
        eps *= 10.0
    return None, None


def lda(X, decoy):
    xbar = X.mean(axis=0)
    Sw = np.zeros((X.shape[1],) * 2); Sb = np.zeros_like(Sw); means = []
    for cls in [True, False]:
        Z = X[decoy == cls]; m = Z.mean(axis=0); Zc = Z - m
        Sw += Zc.T @ Zc / len(Z)
        diff = (m - xbar)[:, None]; Sb += diff @ diff.T; means.append(m)
    Msol, eps = gauss_solve(Sw, Sb)
    if Msol is None:
        return None
    v = xbar / np.linalg.norm(xbar); last = 0.0
    for _ in range(50):
        v1 = Msol @ v; n = np.linalg.norm(v1)
        if abs(n - last) < 1e-8:
            break
        last = n; v = v1 / n
    coef = np.array(means) @ v
    if coef[1] < coef[0]:
        v = -v
    return v if np.all(np.isfinite(v)) else None


def rescore(d, width=40.0, perturb_seed=None):
    """Returns (float32 discriminant scores, 'lda' or 'fallback').
    perturb_seed: multiply every feature by (1 + N(0, 2e-16)) before the fit (diagnostic of the
    ulp-level sensitivity of Sage's solver; the scores are still computed on the unperturbed features)."""
    X, decoy = features(d, width)
    Xf = X if perturb_seed is None else X * (1 + np.random.default_rng(perturb_seed).normal(0, 2e-16, X.shape))
    v = lda(Xf, decoy)
    if v is None:
        with np.errstate(divide='ignore'):
            s = np.log1p((-d['poisson'].to_numpy()).astype(np.float32)) + d['longest_y_pct'].to_numpy(np.float32) / np.float32(3.0)
        return s.astype(np.float32), 'fallback'
    return (X @ v).astype(np.float32), 'lda'


def pooled_q(d, width=40.0, perturb_seed=None):
    """One Sage-equivalent rescoring + picked-peptide q over a pooled PSM table."""
    s, how = rescore(d, width, perturb_seed)
    e = d[['peptide', 'label']].copy(); e['sage_discriminant_score'] = s
    return sage_q.peptide_q(e), s, how
