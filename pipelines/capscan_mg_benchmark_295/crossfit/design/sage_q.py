"""Re-implementation of Sage 0.14.6 peptide-level q values (sage_core::fdr::picked_peptide).

Source read: github.com/lazear/sage tag v0.14.6, crates/sage/src/fdr.rs, ml/kde.rs, ml/mod.rs,
peptide.rs (reverse, Display). Steps, as in Sage:
 1. Competition key per PSM: the peptide string (with modifications) for targets; for decoys the
    decoy string reversed back (residues 1..n-1 reversed together with their modifications, first and
    last residue fixed), i.e. Peptide::reverse().to_string(). generate_decoys is true in every config.
 2. Per key: forward = max target discriminant_score, reverse = max decoy discriminant_score
    (f32::MIN when absent). Competition score = max(forward, reverse); is_decoy = reverse >= forward.
 3. KDE posterior error model (Gaussian kernel, Scott-like bandwidth sigma*(4/3/n)^(1/5) with
    population sigma, separate target and decoy densities weighted by pi = decoy fraction of
    competitions, 1000 bins over [min, max] of competition scores, made monotone from the right,
    linear interpolation between bins).
 4. Rows: one for every present forward entry (target) and every present reverse entry (decoy),
    each with its own score (not only the winner). Sorted by score descending. Running
    decoy = 1 + cumulative PEP (f32), target = cumulative count of target rows (f32); q = decoy/target;
    then reverse cumulative minimum starting at 1.0.
 5. Each PSM gets the q of its own row (its peptide index): targets the forward row of its key,
    decoys the reverse row.
Ties: Sage sorts an unordered hash map with a stable sort, so tie order is arbitrary; here ties are
ordered by key string. Validation against Sage's own peptide_q is in validate.py.
"""
import re
import numpy as np
import pandas as pd

F32_MIN = np.float32(np.finfo(np.float32).min)
TOK = re.compile(r'[A-Z](?:\[[^\]]*\])?')


def reverse_peptide(s):
    nterm = ''; cterm = ''
    m = re.match(r'^(\[[^\]]*\]-)', s)
    if m:
        nterm = m.group(1); s = s[len(nterm):]
    m = re.search(r'(-\[[^\]]*\])$', s)
    if m:
        cterm = m.group(1); s = s[:-len(cterm)]
    toks = TOK.findall(s)
    assert ''.join(toks) == s, s
    n = len(toks) - 1
    if n > 1:
        toks = [toks[0]] + toks[1:n][::-1] + [toks[n]]
    return nterm + ''.join(toks) + cterm


def kde_pdf(sample, grid, bw_adjust=1.0):
    sample = np.asarray(sample, dtype=np.float64)
    n = len(sample)
    sigma = np.sqrt(np.sum((sample - sample.mean()) ** 2) / n)
    h = bw_adjust * sigma * (4.0 / 3.0 / n) ** 0.2
    const = np.sqrt(2 * np.pi) * h * n
    out = np.zeros(len(grid))
    for i in range(0, n, 20000):
        x = sample[i:i + 20000]
        out += np.exp(-0.5 * ((grid[:, None] - x[None, :]) / h) ** 2).sum(axis=1)
    return out / const


def estimator(scores, decoys, bins=1000, bw_adjust=1.0, monotonic=True):
    scores = np.asarray(scores, dtype=np.float64); decoys = np.asarray(decoys, bool)
    pi = decoys.sum() / len(scores)
    lo, hi = scores.min(), scores.max()
    step = (hi - lo) / (bins - 1)
    grid = np.arange(bins) * step + lo
    d = kde_pdf(scores[decoys], grid, bw_adjust) * pi
    t = kde_pdf(scores[~decoys], grid, bw_adjust) * (1 - pi)
    b = d / (t + d)
    if monotonic:
        b = np.maximum.accumulate(b[::-1])[::-1]
    return b, lo, step


def posterior_error(b, lo, step, s):
    s = np.asarray(s, dtype=np.float64)
    raw = np.floor((s - lo) / step)
    lo_i = np.minimum(len(b) - 1, np.where(raw > 0, raw, 0)).astype(np.int64)
    hi_i = np.minimum(len(b) - 1, lo_i + 1)
    lin = (s - (lo_i * step + lo)) / step
    return b[lo_i] + (b[hi_i] - b[lo_i]) * lin


def peptide_q(psms):
    """psms: DataFrame with peptide, label (1/-1), sage_discriminant_score. Returns float32 q per row."""
    lab = psms['label'].to_numpy()
    score = psms['sage_discriminant_score'].to_numpy(dtype=np.float64).astype(np.float32)
    pep = psms['peptide'].to_numpy()
    key = np.where(lab == -1, [reverse_peptide(p) if l == -1 else p for p, l in zip(pep, lab)], pep)
    df = pd.DataFrame(dict(key=key, decoy=lab == -1, score=score))
    g = df.groupby(['key', 'decoy'], sort=True)['score'].max().unstack('decoy')
    fwd = g.get(False, pd.Series(index=g.index, dtype=np.float32)).astype(np.float32)
    rev = g.get(True, pd.Series(index=g.index, dtype=np.float32)).astype(np.float32)
    fwd_v = fwd.fillna(F32_MIN).to_numpy(np.float32); rev_v = rev.fillna(F32_MIN).to_numpy(np.float32)
    comp_score = np.maximum(fwd_v, rev_v); comp_decoy = rev_v >= fwd_v
    b, lo, step = estimator(comp_score.astype(np.float64), comp_decoy)
    keys = g.index.to_numpy()
    rows = []
    has_f = fwd.notna().to_numpy(); has_r = rev.notna().to_numpy()
    rk = np.concatenate([keys[has_f], keys[has_r]])
    rd = np.concatenate([np.zeros(has_f.sum(), bool), np.ones(has_r.sum(), bool)])
    rs = np.concatenate([fwd_v[has_f], rev_v[has_r]]).astype(np.float32)
    order = np.lexsort((rd, rk, -rs.astype(np.float64)))
    rk, rd, rs = rk[order], rd[order], rs[order]
    pe = posterior_error(b, lo, step, rs.astype(np.float64)).astype(np.float32)
    dec = np.add.accumulate(np.concatenate([[np.float32(1.0)], pe]).astype(np.float32), dtype=np.float32)[1:]
    tgt = np.add.accumulate((~rd).astype(np.float32), dtype=np.float32)
    with np.errstate(divide='ignore', invalid='ignore'):
        q = (dec / tgt).astype(np.float32)
    q = np.minimum.accumulate(np.concatenate([q[::-1], []]).astype(np.float32))[::-1]
    q = np.minimum(q, np.float32(1.0))
    qmap = pd.Series(q, index=pd.MultiIndex.from_arrays([rk, rd]))
    out = qmap.reindex(pd.MultiIndex.from_arrays([key, lab == -1])).to_numpy(np.float32)
    assert not np.isnan(out).any()
    return out


def read_psms(path):
    cols = ['scannr', 'peptide', 'label', 'rank', 'sage_discriminant_score', 'peptide_q']
    return pd.read_csv(path, sep='\t', usecols=cols, dtype={'scannr': str, 'peptide': str})
